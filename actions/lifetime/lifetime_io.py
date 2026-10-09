"""Output operations for shared NERSC storage or explicitly local WSL drives.

The default remains the hard-link-based flufl lock used on shared filesystems.
Set ARCUBE_NEARLINE_LOCAL_OUTPUT=1 only when every writer runs on this machine:
WSL external drives may reject hard links and POSIX mode changes. Kernel flock
then supplies crash-safe mutual exclusion without needing either operation.
Do not mix lock backends for the same output directory.
"""
from contextlib import contextmanager
import errno
import os
import time


def local_output():
    """Make the local-only concurrency assumption explicit, never automatic."""
    return os.environ.get('ARCUBE_NEARLINE_LOCAL_OUTPUT') == '1'


def public_permissions(path):
    """Make measurement artifacts readable; honor a local drive's fixed modes.

    Only the explicit local mode tolerates unsupported chmod. File creation,
    content writes, and atomic renames must still succeed normally. Credentials
    must never be placed in these public measurement artifacts.
    """
    try:
        os.chmod(path, 0o644)
    except OSError as exc:
        if not local_output() or exc.errno not in (errno.EPERM, errno.ENOTSUP, errno.EOPNOTSUPP):
            raise


@contextmanager
def output_lock(path, default_timeout, lifetime):
    """Lock one output namespace, retaining the shared-storage default.

    The local backend never unlinks its lock file: unlinking would let a second
    inode be locked while another process still holds the first. The OS releases
    flock when the descriptor closes, including process death. Unlike a leased
    distributed lock, it has no expiry while the process is alive.
    """
    if not local_output():
        from flufl.lock import Lock
        with Lock(path, default_timeout=default_timeout, lifetime=lifetime):
            yield
        return
    import fcntl
    from flufl.lock import TimeOutError
    timeout = default_timeout.total_seconds()
    deadline = time.monotonic() + timeout
    with open(path, 'a') as stream:
        while True:
            try:
                fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if time.monotonic() >= deadline:
                    raise TimeOutError('Timed out waiting for local output lock: ' + path)
                time.sleep(min(.1, max(0., deadline-time.monotonic())))
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)
