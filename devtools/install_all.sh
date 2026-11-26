#!/bin/bash -l

# Check if Micromamba is installed
if [ -x "$(command -v micromamba)" ]; then
    echo "Micromamba is installed."
    COMMAND_PKG=micromamba
# Check if Mamba is installed
elif [ -x "$(command -v mamba)" ]; then
    echo "Mamba is installed."
    COMMAND_PKG=mamba
# Check if Conda is installed
elif [ -x "$(command -v conda)" ]; then
    echo "Conda is installed."
    COMMAND_PKG=conda
else
    echo "Micromamba, Mamba, and Conda are not installed. Please download and install one of them - we strongly recommend Micromamba or Mamba."
    exit 1
fi

# Set up Conda/Micromamba environment
if [ "$COMMAND_PKG" == "micromamba" ]; then
    eval "$(micromamba shell hook --shell=bash)"
    micromamba activate base
    BASE=$MAMBA_ROOT_PREFIX
    # shellcheck source=/dev/null
    source "$BASE/etc/profile.d/micromamba.sh"
else
    BASE=$(conda info --base)
    # shellcheck source=/dev/null
    source "$BASE/etc/profile.d/conda.sh"
fi

# Create or update rmg_env from T3's own spec (package mode: rmgpy, arkane and the RMG database all
# come from the pinned conda packages, so no RMG-Py or RMG-database clone is needed or compiled).
# devtools/check_rmg_env.py, run at the end, asserts the result matches the spec.
if { $COMMAND_PKG env list | grep 'rmg_env'; } >/dev/null 2>&1; then
    $COMMAND_PKG env update -n rmg_env -f devtools/rmg_env.yml
else
    $COMMAND_PKG env create -f devtools/rmg_env.yml
fi;

# Temporarily change directory to parent to install sibling repos
pushd .

# Create or update t3_env
cd ..
cd T3 || { echo "Cannot find T3 directory"; exit 1; }
if { $COMMAND_PKG env list | grep 't3_env'; } >/dev/null 2>&1; then
    $COMMAND_PKG env update -n t3_env -f environment.yml
else
    $COMMAND_PKG env create -f environment.yml
fi;

# Switch to t3_env and clone/install ARC
if [ "$COMMAND_PKG" == "micromamba" ]; then
    micromamba activate t3_env
else
    conda activate t3_env
fi

cd ..
if [ -d "./ARC" ]; then
    cd ARC || exit
    git pull https://github.com/ReactionMechanismGenerator/ARC
else
    git clone https://github.com/ReactionMechanismGenerator/ARC
    cd ARC || exit
fi;

PYTHONPATH=$PYTHONPATH:$(pwd)
export PYTHONPATH
echo "export PYTHONPATH=$PYTHONPATH:$(pwd)" >> ~/.bashrc
if { $COMMAND_PKG env list | grep 'arc_env'; } >/dev/null 2>&1; then
    $COMMAND_PKG env update -n arc_env -f environment.yml
else
    $COMMAND_PKG env create -f environment.yml
fi;

if [ "$COMMAND_PKG" == "micromamba" ]; then
    micromamba activate arc_env
else
    conda activate arc_env
fi

make install-all

conda deactivate

# Restore original directory
popd || exit

# Install PyRDL into t3_env (delegates to ARC's install_pyrdl.sh)
bash devtools/install_pyrdl.sh

# Assert that rmg_env is the pinned environment, reached the way T3 reaches it (package mode).
if [ "$COMMAND_PKG" == "micromamba" ]; then
    micromamba activate t3_env
else
    conda activate t3_env
fi
python devtools/check_rmg_env.py || exit 1
