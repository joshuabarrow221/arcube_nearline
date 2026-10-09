#!/usr/bin/env python3
"""Queue completed files, optionally tracking lifetime input generations.

Ordinary actions retain their path-based discovery. The opt-in lifetime profile
also watches native packet companions and slow-control snapshots. A late packet
or replaced input gets a new revision of the same acquisition; the action's
persistent shards make these retries idempotent. Dry runs never contact MongoDB.
"""
import argparse
import glob
from pathlib import Path
import shlex
import sys
import time
from typing import Callable, Optional


class Watcher:
    def __init__(self, prog: Path, paths: list[Path], exts: list[str],
                 cond: Callable[[Path], bool], sleep_between_scans: int,
                 min_file_age: int, max_file_age: Optional[int], priority: int,
                 dry_run: bool, lifetime_config=None, launchpad=None):
        self.prog = Path(prog)
        self.paths, self.exts, self.cond = paths, exts, cond
        self.sleep_between_scans = sleep_between_scans
        self.min_file_age, self.max_file_age = min_file_age, max_file_age
        self.priority, self.dry_run = priority, dry_run
        self.seen_files = set()
        self.config_path = Path(lifetime_config).resolve() if lifetime_config else None
        self.config = None
        if self.config_path:
            # Discovery deliberately imports only a standard-library module;
            # a read-only dry run needs neither numerical libraries nor MongoDB.
            sys.path.insert(0, str(Path(__file__).parent/'actions'/'lifetime'))
            from monitor_config import load_config
            self.config = load_config(self.config_path)
        self.lpad, self.db = None, None
        if not dry_run:
            if launchpad is None:
                from fireworks import LaunchPad
                launchpad = LaunchPad.auto_load()
            self.lpad = launchpad
            self.db = launchpad.connection[launchpad.name]

    def candidate(self, path):
        """Resolve canonical input and ensure *all* companion files have settled."""
        p = Path(path)
        if self.config:
            from monitor_config import discovery
            item = discovery(p, self.config)
            if item is None:
                return None
            ages = [time.time()-d['mtime_ns']/1e9 for d in item['dependencies']]
            # The newest dependency determines freshness. Thus an old FLOW can
            # be retried when its native file arrives today. Every dependency
            # must pass the minimum age; no half-transferred companion is read.
            if min(ages) < self.min_file_age:
                return None
            if self.max_file_age is not None and min(ages) > self.max_file_age:
                return None
            item['script'] = shlex.join([str(self.prog), item['input'],
                                        '--monitor-config', str(self.config_path)])
            item['seen_key'] = (item['key'], item['revision'])
            return item
        if not self.cond(p):
            return None
        try:
            if not p.is_file():
                return None
            age = time.time()-p.stat().st_mtime
        except FileNotFoundError:
            return None  # A scan may race a producer's rename or removal.
        if age < self.min_file_age or (self.max_file_age is not None and age > self.max_file_age):
            return None
        return dict(input=str(p), script=shlex.join([str(self.prog), str(p)]), seen_key=p)

    def find_firework(self, p, item=None):
        """Task existence is scheduling state, never proof of a successful fit.

        Failed tasks are not automatically requeued every scan. Operators use
        FireWorks' explicit rerun command after repairing an extraction/input
        problem. New source/config revisions can be scheduled independently.
        """
        if self.db is None:
            return None
        item = item or self.candidate(p)
        if item is None:
            return None
        if self.config:
            query = {'spec._purity_key': item['key'], 'spec._purity_revision': item['revision']}
        else:
            # Recognize pre-upgrade unquoted tasks as well as new safely quoted
            # commands, preserving existing generic watcher deduplication.
            query = {'spec._tasks.script': {'$in': [item['script'], f'{self.prog} {p}']}}
        return self.db['fireworks'].find_one(query)

    def maybe_make_firework(self, p: Path):
        item = self.candidate(p)
        if item is None or item['seen_key'] in self.seen_files:
            return None
        if self.find_firework(p, item):
            self.seen_files.add(item['seen_key'])
            return None
        if self.dry_run:
            print('Would queue: '+item['script'], flush=True)
            self.seen_files.add(item['seen_key'])
            return item
        from fireworks import Firework, ScriptTask
        spec = {'_priority': self.priority}
        if self.config:
            spec.update(_purity_key=item['key'], _purity_revision=item['revision'])
        fw = Firework(ScriptTask.from_str(item['script']), name=self.prog.name, spec=spec)
        self.lpad.add_wf(fw)
        self.seen_files.add(item['seen_key'])
        print(f'FW {fw.fw_id:05}: {item["input"]}', flush=True)
        return fw

    def snarf(self):
        if self.config_path:
            from monitor_config import load_config
            self.config = load_config(self.config_path)  # Notice deliberate config edits.
            roots = [('flow_root', '*.FLOW.hdf5'), ('packet_root', 'packet-*.h5'),
                     ('snapshot_root', self.config['snapshot_pattern'])]
            for key, pattern in roots:
                root = Path(self.config[key])
                # Snapshots have a flat, explicit export contract. Acquisition
                # files keep their complete mirrored run/subrun hierarchy.
                paths = root.glob(pattern) if key == 'snapshot_root' else root.rglob(pattern)
                for path in sorted(paths):
                    self.maybe_make_firework(path)
            return
        for path in self.paths:
            for actual_path in glob.glob(str(path)):
                for ext in self.exts:
                    for p in sorted(Path(actual_path).rglob(f'*.{ext}')):
                        self.maybe_make_firework(p)

    def watch_inotify(self):
        """Optional acceleration; periodic scans remain the portable default.

        Close events alone can precede min_file_age. The periodic scan below
        supplies the delayed retry and also sees files created on another host.
        """
        from watchdog.observers import Observer
        from watchdog.events import FileSystemEventHandler
        watcher = self
        class Handler(FileSystemEventHandler):
            def on_closed(self, event):
                if not event.is_directory:
                    watcher.maybe_make_firework(Path(event.src_path))
        observer = Observer()
        roots = [self.config[k] for k in ('flow_root', 'packet_root', 'snapshot_root')] if self.config else self.paths
        for root in roots:
            observer.schedule(Handler(), str(root), recursive=True)
        observer.start()
        try:
            self.watch_dumb()
        finally:
            observer.stop()
            observer.join()

    def watch_dumb(self):
        while True:
            self.snarf()
            print(f'NAPPING FOR {self.sleep_between_scans} SECONDS', flush=True)
            time.sleep(self.sleep_between_scans)


def is_hdf5(p: Path):
    return p.suffix.lower() in ['.h5', '.hdf5']


def is_raw_binary(p: Path):
    return is_hdf5(p) and 'binary-' in p.name


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('prog', type=Path)
    locations = ap.add_mutually_exclusive_group(required=True)
    locations.add_argument('--path', type=Path)
    locations.add_argument('--path-file')
    locations.add_argument('--lifetime-config', type=Path,
                           help='Watch FLOW, native packets and snapshots from one deployment config')
    ap.add_argument('--require-binary', action='store_true')
    ap.add_argument('--ext', help='File extension (if not h5 / hdf5)')
    ap.add_argument('--base-dir', type=Path, default=Path('/'))
    ap.add_argument('--sleep-between-scans', type=int, default=30)
    ap.add_argument('--min-file-age', type=int, default=30,
                    help='Minimum seconds since last write; producers must publish immutable files by rename')
    ap.add_argument('--max-file-age', type=int, default=None)
    ap.add_argument('--priority', type=int, default=0)
    ap.add_argument('--dry-run', action='store_true', help='Print commands without contacting FireWorks/MongoDB')
    ap.add_argument('--once', action='store_true', help='Scan once and exit (useful for validation/backfill)')
    args = ap.parse_args()
    if args.min_file_age < 0 or args.sleep_between_scans <= 0 or (args.max_file_age is not None and args.max_file_age < 0):
        ap.error('ages must be nonnegative and scan interval must be positive')
    if args.lifetime_config and (args.require_binary or args.ext):
        ap.error('the lifetime profile defines its own input types')
    paths = [args.path] if args.path else []
    if args.path_file:
        paths = [args.base_dir/line.strip() for line in Path(args.path_file).read_text().splitlines()
                 if line.strip() and not line.lstrip().startswith('#')]
    exts, cond = ['h5', 'hdf5'], is_hdf5
    if args.require_binary:
        cond = is_raw_binary
    elif args.ext:
        exts = [args.ext]
        cond = lambda p: p.suffix.lower() == f'.{args.ext.lower()}'
    watcher = Watcher(args.prog.absolute(), paths, exts, cond, args.sleep_between_scans,
                      args.min_file_age, args.max_file_age, args.priority, args.dry_run,
                      lifetime_config=args.lifetime_config)
    watcher.snarf() if args.once else watcher.watch_dumb()


if __name__ == '__main__':
    main()
