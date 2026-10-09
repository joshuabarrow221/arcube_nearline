# arcube_nearline

Nearline prompt processing workflows for DUNE's ArgonCube 2x2 demonstrator.

## Installation

``` bash
admin/install.sh
```

The installation will go into `_install`. After installing, edit `config/my_launchpad.yaml` to configure the MongoDB connection.

## Environment setup

``` bash
source admin/load.sh
```

If your `my_launchpad.yaml` points to a new Mongo database, run `lpad reset` to initialize it.

## Monitoring the filesystem and loading the DB

``` bash
./watcher.py actions/packetize.sh --path /global/cfs/cdirs/dune/www/data/2x2/CRS/commission/April2024
```

This will continue to run in the foreground.
Use `--dry-run --once` to inspect discovery without connecting to MongoDB.
Generic actions deduplicate by action and file path. The opt-in lifetime profile
tracks FLOW/native input revisions and slow-control snapshots; see below.

## Launching a worker loop

``` bash
cd $(mktemp -d)
rlaunch rapidfire
```

The temporary directory avoids littering the filesystem. Use `rlaunch rapidfire --nlaunches infinite` to prevent the worker from exiting when there's no available work.

## Purity monitoring

The [lifetime README](actions/lifetime/README.md) describes the recent-history
monitoring window: daily PRM means, six-hour pooled track/packet fits and six-hour
gas equivalents. The interactive plot supports zoom and links to fit diagnostics.
All tracks includes the cosmic-enriched subset; those two curves share data.

Use `watcher.py --lifetime-config ...` with the lifetime action to process incoming
FLOW, native packets and received slow-control snapshots. The watcher does not
connect to DAQ05 itself: companion export/rsync scripts supply those snapshots.
See the [operator guide](docs/lifetime-integration.md) for configuration and
scheduling, and [methods reference](docs/lifetime-methods.md) for inner workings.
