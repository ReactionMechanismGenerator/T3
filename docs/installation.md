# Installation

T3 was tested on **Linux** (Ubuntu 22.04 LTS and later).


## Prerequisites

### 1. Install system packages

You need `git`, a C/C++ compiler, and either `curl` or `wget`.

- **Ubuntu / Debian**:

    ```bash
    sudo apt install git gcc g++ make curl
    ```

- **Fedora / Red Hat**:

    ```bash
    sudo dnf install git gcc gcc-c++ make curl
    ```

- **macOS** (with [Homebrew](https://brew.sh)):

    ```bash
    xcode-select --install
    ```


### 2. Install a conda package manager

We recommend [Miniforge](https://github.com/conda-forge/miniforge)
(which includes `mamba` for faster environment solving):

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

Follow the prompts and restart your terminal when done.


## Install T3 and dependencies

T3 requires four co-installed components:

| Component | Conda env | Why T3 needs it |
|---|---|---|
| **T3** itself | `t3_env` | the orchestrator |
| **ARC** | `t3_env` | QM job submission and post-processing |
| **RMG-Py** | `rmg_env` | mechanism generation and Arkane PDep SA — T3 spawns these as subprocesses |
| **RMG-database** | `rmg_env` | thermo / kinetics / transport data used by both RMG and ARC |

`rmg_env` is pinned in one file, `devtools/rmg_env.yml`: the `rmg` release, its
Python, the `rmgdatabase` package of the same release, and a floor on
`numdifftools`. Both installation options below and T3's CI build `rmg_env`
from that file. RMG-Py and RMG-database are installed as conda packages, so
neither repository is cloned or compiled.

!!! Note "Why two conda environments?"
    T3 itself never imports `rmg-py`; it only shells out to RMG and Arkane.
    Keeping `rmg_env` separate from `t3_env` lets the two have very different
    dependency trees without conflicting, and lets you upgrade T3 without
    rebuilding RMG.  Both environments must exist for T3 to run end-to-end.

Two installation options are available.


### Option 1: Automated installation

Clone T3 and install everything with one command:

```bash
git clone https://github.com/ReactionMechanismGenerator/T3.git
cd T3
make install
```

This creates `rmg_env` from `devtools/rmg_env.yml`, creates `t3_env`, clones
ARC into the parent folder if it is not already there and compiles it, and
installs a custom updated version of PyRDL. Last, it runs
`devtools/check_rmg_env.py`, which stops the install with a list of failures
if `rmg_env` does not match `devtools/rmg_env.yml` (see
[Keep `RMG_PATH` unset](#keep-rmg_path-unset)).


### Option 2: Manual installation

#### Clone repositories

Clone T3 and ARC (e.g., under `~/Code/`). RMG-Py and RMG-database come as
conda packages in `rmg_env`, so do not clone them:

```bash
cd ~/Code
git clone https://github.com/ReactionMechanismGenerator/T3.git
git clone https://github.com/ReactionMechanismGenerator/ARC.git
```

#### Create the T3 conda environment

```bash
cd ~/Code/T3
mamba env create -f environment.yml
conda activate t3_env
```

#### Install PyRDL (required by ARC)

```bash
cd ~/Code/ARC
bash devtools/install_pyrdl.sh t3_env
```

#### Compile ARC extensions

```bash
cd ~/Code/ARC
make compile
```

#### Create the RMG conda environment

T3 spawns RMG and Arkane as subprocesses inside `rmg_env`, so this step is
**required** even though T3 itself never imports `rmg-py`. Create it from T3's
spec file, and add no other package to it:

```bash
cd ~/Code/T3
mamba env create -f devtools/rmg_env.yml
```

!!! Warning "`rmg_env` needs `numdifftools>=0.9.41` from conda-forge"
    Arkane imports `numdifftools`. Version 0.9.39, the one on the `rmg`
    channel, calls `np.finfo().machar`, which the numpy that RMG uses no
    longer has, so with it `import arkane` fails. `devtools/rmg_env.yml`
    already asks for the conda-forge version. If you built `rmg_env` another
    way and see an Arkane import error that mentions `machar`, install the
    newer version into `rmg_env`:

    ```bash
    mamba install -n rmg_env "conda-forge::numdifftools>=0.9.41"
    ```

#### Set environment variables

Add the following to your `~/.bashrc` (or `~/.zshrc` on macOS):

```bash
export PYTHONPATH=$PYTHONPATH:~/Code/T3:~/Code/ARC
```

Then reload your shell:

```bash
source ~/.bashrc
```

Do not set `RMG_DB_PATH`: ARC finds the database that the `rmgdatabase`
package installed in `rmg_env`, which is the one RMG itself reads.

#### Keep `RMG_PATH` unset

ARC runs RMG and Arkane in `rmg_env`. If ARC's `RMG_PATH` setting names an
RMG-Py source tree, ARC puts that tree on the child's `PYTHONPATH`, and it
replaces the pinned `rmgpy` package. ARC sets `RMG_PATH` by itself when it
finds any of these:

- an exported `RMG_PATH` variable;
- a `PYTHONPATH` entry that contains `RMG-Py`;
- a folder `~/Code/RMG-Py` (or `$GITHUB_WORKSPACE/RMG-Py`).

So do not export `RMG_PATH`, and remove RMG-Py from `PYTHONPATH`. If you keep
an RMG-Py clone at `~/Code/RMG-Py` for other work, add this line to
`~/.arc/settings.py` (create the file if it does not exist):

```python
RMG_PATH = None
```

To confirm that `rmg_env` is the pinned environment and that ARC reaches it in
package mode, run the same check that `make install` and CI run:

```bash
conda activate t3_env
cd ~/Code/T3
python devtools/check_rmg_env.py
```

It ends with `rmg_env matches .../devtools/rmg_env.yml (package mode).`, or
lists each failure and exits 1.


## Verify the installation

Activate the environment and run the test suite:

```bash
conda activate t3_env
cd ~/Code/T3
pytest tests/ -x -q
```

All tests should pass. If you encounter import errors, verify that
`PYTHONPATH` is set correctly and run `python devtools/check_rmg_env.py`.


## Optional: convenient aliases

Add these to your `~/.bashrc` (adjust paths accordingly):

```bash
export t3_path=$HOME'/Code/T3'
alias t3e='conda activate t3_env'
alias t3='python $t3_path/T3.py input.yml'
alias t3code='cd $t3_path'
```

After sourcing `~/.bashrc` or restarting the terminal:

- `t3e` activates the environment
- `t3code` changes to the T3 repository directory
- `t3` runs T3 with `input.yml` in the current folder


## Updating T3

Pull the latest changes from the main branch:

```bash
cd ~/Code/T3
git pull origin main
```

You may also want to update ARC, and recompile it:

```bash
cd ~/Code/ARC && git pull origin main
make compile
```

When `devtools/rmg_env.yml` changes, rebuild `rmg_env` from it:

```bash
cd ~/Code/T3
conda env remove -n rmg_env
mamba env create -f devtools/rmg_env.yml
```
