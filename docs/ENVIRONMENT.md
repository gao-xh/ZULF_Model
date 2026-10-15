# Environment

## conda environment `zulf` (local machine)

ZULF_Model has its own conda environment, `zulf` (Python 3.12, conda-forge, Miniforge). Its definition lives in
the owner's environment repository (`~/code/ml-env/zulf.yml`, with the exact lock `zulf-explicit-lock.txt`);
the packages are:

| Package | Used by |
|---|---|
| python 3.12, numpy, scipy | everything (zulf_core, zulf_hypothesis, zulf_processing, zulf_studio) |
| pytorch | zulf_model (models, training) |
| h5py | training shards |
| matplotlib | figures, scripts, Studio plots |
| rdkit | structure insets of scripts/paper_figure.py (optional) |
| pyside6 | ZULF Studio window (zulf_studio/app.py) |
| anthropic, openai | ZULF Studio AI assistant (zulf_studio/assistant.py) |
| pip, ipykernel | editable install, Jupyter kernel |

Create and install:

    mamba env create -f ~/code/ml-env/zulf.yml          # or: conda create -n zulf -c conda-forge python=3.12 ...
    conda activate zulf
    pip install --no-deps -e ~/code/ZULF_Model           # the project itself; dependencies come from conda

Update after a package is added to `zulf.yml` (no `--prune` needed; check that nothing else changes version):

    mamba env update -y -n zulf -f ~/code/ml-env/zulf.yml
    conda list -n zulf --explicit > ~/code/ml-env/zulf-explicit-lock.txt

After a new top-level package is added to this repository (for example `zulf_studio`), run the editable
install again so it is importable outside the repository directory.

In VS Code and Jupyter choose the interpreter / kernel "zulf (Python 3.12)". Project dependencies are not
upgraded casually (reproducibility); the environment changes only with an entry in `zulf.yml`.

Calling the environment's Python without activating it (as agents do):
`/opt/homebrew/Caskroom/miniforge/base/envs/zulf/bin/python`.

## Without conda (other machines, cloud sessions)

    python -m venv .venv && . .venv/bin/activate
    pip install -e ".[dev]"                 # numpy, scipy, torch, h5py, matplotlib
    pip install -e ".[studio,ai]"           # PySide6 and matplotlib for Studio; anthropic and openai for its assistant

The extras are defined in `pyproject.toml`.

## Environment variables

| Variable | Meaning |
|---|---|
| `ZULF_DATA_DIR` | folder of averaged FIDs for the regression and processing scripts. The `zulf` environment sets it to `~/research/zulf/data/raw` on activation; FIDs are found by file name in its measurement subfolders (`zulf_core.io.find_data_file`). Not set in `~/.zshrc`: without the activated environment, prefix commands with it. |
| `ZULF_MODEL_WORKSPACE` | output folder of the agent tools (default `runs/`, not in git) |
| `OMP_NUM_THREADS` | set to 1 for parallel fits (`--workers`); fit_joint_series subprocesses started by Studio set it |
| `ANTHROPIC_API_KEY` / `ANTHROPIC_AUTH_TOKEN` | Claude API credentials for the Studio assistant (or `ant auth login`) |
| `OPENAI_API_KEY`, `OPENAI_MODEL` | OpenAI credentials and default model for the Studio assistant |

Keys are never written into the repository or its configs.

## Data layout

Experimental data is not in git; on the local machine it is under `~/research/zulf/data/` (`original/`, `raw/`,
`processed/`, one folder per measurement; README "Data" and that folder's README.md).

## Checks

    python scripts/check_ascii.py                         # code and docs are ASCII English
    python -m unittest tests.test_studio                  # Studio, its API and assistant (about 3 s; Qt offscreen)
    python -m unittest tests.test_processing tests.test_j_tuner     # quick check, about 1.5 min
