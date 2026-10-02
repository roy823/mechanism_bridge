"""Runtime provenance for run manifests: code version, interpreter, host, packages, job."""
import hashlib
import importlib.metadata
import os
from pathlib import Path
import platform
import subprocess
import sys

PACKAGES = ('numpy', 'scipy', 'ase', 'rdkit', 'networkx', 'torch', 'aimnet', 'sella', 'pyscf')
# Code and configuration that determine a run; tracked report artifacts are excluded
# so that a dirty reports/ tree neither marks the code dirty nor slows the diff.
CODE_PATHS = ('src', 'scripts', 'configs', 'tests', 'pyproject.toml')
SCHEDULER_VARIABLES = ('SLURM_JOB_ID', 'SLURM_ARRAY_JOB_ID', 'SLURM_ARRAY_TASK_ID',
                       'SLURM_JOB_PARTITION', 'SLURM_JOB_NODELIST', 'SLURM_CPUS_PER_TASK',
                       'CUDA_VISIBLE_DEVICES')


def _git(root, *args):
    """(stdout bytes, None) on success, (None, reason) otherwise."""
    try:
        done = subprocess.run(['git', '-C', str(root), *args], capture_output=True, timeout=60)
    except (OSError, subprocess.SubprocessError) as exc:
        return None, f'{type(exc).__name__}: {exc}'
    if done.returncode:
        return None, done.stderr.decode(errors='replace').strip() or f'git exit {done.returncode}'
    return done.stdout, None


def _version(name):
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def git_provenance(root):
    """Commit and tracked-file state of the repository whose top level is `root`.

    `dirty` covers tracked files under CODE_PATHS only; `diff_sha256` hashes
    `git diff HEAD` over the same paths, so a dirty run can be matched to a saved
    patch. Untracked files are not covered; the per-run source snapshot remains
    the authoritative copy of src/. A root
    nested inside another repository, an rsync copy without .git, or git's
    safe.directory refusal all give null fields plus `error`.
    """
    empty = dict(commit=None, dirty=None, diff_sha256=None)
    toplevel, error = _git(root, 'rev-parse', '--show-toplevel')
    if error:
        return dict(empty, error=error)
    if Path(toplevel.decode().strip()).resolve() != Path(root).resolve():
        return dict(empty, error=f'run root is not the repository top level ({toplevel.decode().strip()})')
    commit, error = _git(root, 'rev-parse', 'HEAD')
    status, status_error = _git(root, 'status', '--porcelain', '--untracked-files=no', '--', *CODE_PATHS)
    diff, diff_error = _git(root, 'diff', 'HEAD', '--binary', '--', *CODE_PATHS)
    error = error or status_error or diff_error
    info = dict(commit=commit.decode().strip() if commit else None,
                dirty=None if status is None else bool(status.strip()),
                diff_sha256=hashlib.sha256(diff).hexdigest() if diff else None,
                scope=list(CODE_PATHS))
    return dict(info, error=error) if error else info


def runtime_provenance(root, argv=None):
    """Identify exactly which code and environment produced a run."""
    return dict(
        git=git_provenance(root),
        argv=list(sys.argv if argv is None else argv),
        python=dict(version=sys.version, executable=sys.executable,
                    implementation=platform.python_implementation()),
        host=dict(node=platform.node(), platform=platform.platform(), machine=platform.machine(),
                  cpu_count=os.cpu_count()),
        packages={name: _version(name) for name in PACKAGES},
        scheduler={name: os.environ.get(name) for name in SCHEDULER_VARIABLES})
