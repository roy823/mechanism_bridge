import shutil
import subprocess

import pytest

from mechbridge.provenance import runtime_provenance


@pytest.mark.skipif(shutil.which('git') is None, reason='git is not installed')
def test_runtime_provenance_tracks_commit_and_dirty_state(tmp_path):
    def git(*args):
        subprocess.run(['git', '-C', str(tmp_path), *args], check=True, capture_output=True)
    git('init', '-q')
    (tmp_path/'src').mkdir()
    (tmp_path/'src'/'module.py').write_text('x = 1\n')
    (tmp_path/'notes.md').write_text('outside the code paths\n')
    git('add', 'src/module.py', 'notes.md')
    git('-c', 'user.name=test', '-c', 'user.email=test@example.invalid', '-c', 'commit.gpgsign=false',
        'commit', '-qm', 'initial')
    clean = runtime_provenance(tmp_path, argv=['run_network_exploration.py', '--seed', '17'])
    assert len(clean['git']['commit']) == 40
    assert clean['git']['dirty'] is False and clean['git']['diff_sha256'] is None
    assert clean['argv'] == ['run_network_exploration.py', '--seed', '17']
    assert 'error' not in clean['git']
    (tmp_path/'notes.md').write_text('edited report text\n')
    assert runtime_provenance(tmp_path)['git']['dirty'] is False
    (tmp_path/'src'/'module.py').write_text('x = 2\n')
    dirty = runtime_provenance(tmp_path)
    assert dirty['git']['commit'] == clean['git']['commit']
    assert dirty['git']['dirty'] is True and len(dirty['git']['diff_sha256']) == 64
    (tmp_path/'nested').mkdir()
    nested = runtime_provenance(tmp_path/'nested')['git']
    assert nested['commit'] is None and 'top level' in nested['error']


def test_runtime_provenance_outside_repository(tmp_path):
    info = runtime_provenance(tmp_path/'missing')
    assert info['git']['commit'] is None and info['git']['dirty'] is None
    assert info['git']['error']
    assert info['packages']['numpy']
    assert 'SLURM_JOB_ID' in info['scheduler']
