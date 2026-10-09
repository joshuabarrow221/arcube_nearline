"""Exercise external-drive portability without weakening shared-file safety."""
import errno
import os
from pathlib import Path
import subprocess
import sys
from datetime import timedelta

import pytest

SOURCE=Path(__file__).resolve().parents[1]/'actions/lifetime'
sys.path.insert(0,str(SOURCE))
from lifetime_io import output_lock, public_permissions


def test_fixed_drive_modes_require_explicit_local_opt_in(tmp_path,monkeypatch):
    target=tmp_path/'public.json';target.write_text('{}')
    def unsupported(*args):raise PermissionError(errno.EPERM,'drive has fixed modes')
    monkeypatch.setattr(os,'chmod',unsupported)
    monkeypatch.delenv('ARCUBE_NEARLINE_LOCAL_OUTPUT',raising=False)
    with pytest.raises(PermissionError):public_permissions(target)
    monkeypatch.setenv('ARCUBE_NEARLINE_LOCAL_OUTPUT','1')
    public_permissions(target)
    def denied(*args):raise PermissionError(errno.EACCES,'real access denial')
    monkeypatch.setattr(os,'chmod',denied)
    with pytest.raises(PermissionError):public_permissions(target)


def test_local_lock_excludes_another_process_and_releases(tmp_path,monkeypatch):
    monkeypatch.setenv('ARCUBE_NEARLINE_LOCAL_OUTPUT','1')
    lock=str(tmp_path/'output.lock')
    code='''from lifetime_io import output_lock
from datetime import timedelta
import sys
with output_lock(sys.argv[1],timedelta(seconds=.15),timedelta(seconds=1)):
    print('acquired')
'''
    env=dict(os.environ,PYTHONPATH=str(SOURCE))
    with output_lock(lock,timedelta(seconds=1),timedelta(seconds=1)):
        child=subprocess.run([sys.executable,'-c',code,lock],env=env,capture_output=True,text=True)
        assert child.returncode!=0 and 'TimeOutError' in child.stderr
    child=subprocess.run([sys.executable,'-c',code,lock],env=env,capture_output=True,text=True)
    assert child.returncode==0 and 'acquired' in child.stdout
    assert Path(lock).exists() # Stable inode across successive lock owners.
