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
| **RMG-database** | (filesystem only) | thermo / kinetics / transport data used by both RMG and ARC |

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

This clones ARC, RMG-Py, and RMG-database into the parent folder if they
are not already there, creates both `t3_env` and `rmg_env`, compiles
RMG-Py and ARC, and installs a custom updated version of PyRDL.
It may take a while (RMG-Py's Cython compile is the slow step).


### Option 2: Manual installation

#### Clone repositories

Clone all four repositories (e.g., under `~/Code/`):

```bash
cd ~/Code
git clone https://github.com/ReactionMechanismGenerator/T3.git
git clone https://github.com/ReactionMechanismGenerator/ARC.git
git clone https://github.com/ReactionMechanismGenerator/RMG-Py.git
git clone https://github.com/ReactionMechanismGenerator/RMG-database.git
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

#### Create the RMG conda environment and compile RMG-Py

T3 spawns RMG and Arkane as subprocesses inside `rmg_env`, so this step is
**required** even though T3 itself never imports `rmg-py`.

```bash
cd ~/Code/RMG-Py
mamba env create -f environment.yml
conda activate rmg_env
make
conda deactivate
```

#### Set environment variables

Add the following to your `~/.bashrc` (or `~/.zshrc` on macOS):

```bash
export PYTHONPATH=$PYTHONPATH:~/Code/T3:~/Code/ARC:~/Code/RMG-Py
export RMG_DB_PATH=~/Code/RMG-database
```

Then reload your shell:

```bash
source ~/.bashrc
```


## Verify the installation

Activate the environment and run the test suite:

```bash
conda activate t3_env
cd ~/Code/T3
pytest tests/ -x -q
```

All tests should pass. If you encounter import errors, verify that
`PYTHONPATH` and `RMG_DB_PATH` are set correctly.


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

You may also want to update ARC, RMG-Py, and RMG-database:

```bash
cd ~/Code/ARC && git pull origin main
cd ~/Code/RMG-Py && git pull origin main
cd ~/Code/RMG-database && git pull origin main
```

Recompile ARC and RMG-Py with the latest changes:

```bash
cd ~/Code/ARC
make compile

cd ~/Code/RMG-Py
conda activate rmg_env
make
conda deactivate
```
