#!/bin/bash

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
fi

# Set up Conda/Micromamba environment
if [ "$COMMAND_PKG" == "micromamba" ]; then
    MICROMAMBA_BASE=$(micromamba info --base)
    # shellcheck source=/dev/null
    source "$MICROMAMBA_BASE/etc/profile.d/micromamba.sh"
else
    CONDA_BASE=$(conda info --base)
    # shellcheck source=/dev/null
    source "$CONDA_BASE/etc/profile.d/conda.sh"
fi

# Temporarily change directory to $HOME to install software
pushd .
cd ..

if [ -d "./RMG-database" ]; then
    cd RMG-database || exit
    git pull https://github.com/ReactionMechanismGenerator/RMG-database
    cd ..
else
    git clone https://github.com/ReactionMechanismGenerator/RMG-database
fi;

if [ -d "./RMG-Py" ]; then
    cd RMG-Py || exit
    git pull https://github.com/ReactionMechanismGenerator/RMG-Py
else
    git clone https://github.com/ReactionMechanismGenerator/RMG-Py
    cd RMG-Py || exit
fi;

echo "export PYTHONPATH=$PYTHONPATH:$(pwd)" >> ~/.bashrc
echo "export PATH=$PATH:$(pwd)" >> ~/.bashrc 
PYTHONPATH=$PYTHONPATH:$(pwd)
export PYTHONPATH
PATH=$PATH:$(pwd)
export PATH
# shellcheck source=~/.bashrc
source ~/.bashrc

#Creating rmg_env environment (or updating if it already exists)
if { $COMMAND_PKG env list | grep 'rmg_env'; } >/dev/null 2>&1; then
    $COMMAND_PKG env update -n rmg_env -f environment.yml
else
    $COMMAND_PKG env create -f environment.yml
fi;

#Prior to activation of the environment, the LD_LIBRARY_PATH needs to be set as an environmnet variable when rmg_env is activated.
#This exporting and unsetting will solve the RMS installation during the Julia compile and also when the environment is deactivated
#the original LD_LIBRARY_PATH is set.
echo "export OLD_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >> "$CONDA_BASE/envs/rmg_env/etc/conda/activate.d/env_vars.sh"
echo "export LD_LIBRARY_PATH=$CONDA_BASE/envs/rmg_env/lib:$LD_LIBRARY_PATH" >> "$CONDA_BASE/envs/rmg_env/etc/conda/activate.d/env_vars.sh"
echo "export LD_LIBRARY_PATH=${OLD_LD_LIBRARY_PATH}" >> "$CONDA_BASE/envs/rmg_env/etc/conda/deactivate.d/env_vars.sh"
echo "unset OLD_LD_LIBRARY_PATH" >> "$CONDA_BASE/envs/rmg_env/etc/conda/deactivate.d/env_vars.sh"

#Active rmg_env environment
if [ "$COMMAND_PKG" == "micromamba" ]; then
    micromamba activate rmg_env
else
    conda activate rmg_env
fi

#Compile RMG-Py
make
#Update pyjulia to the latest version
$COMMAND_PKG update pyjulia -y

#Ensure that added paths etc. are set and then reactivate rmg_env
. ~/.bashrc
if [ "$COMMAND_PKG" == "micromamba" ]; then
    micromamba activate rmg_env
else
    conda activate rmg_env
fi

### Install python + julia connection and compile RMS
#
#This code here is for updating julia if ever required:     $(julia -e 'using Pkg; Pkg.add("UpdateJulia"); using UpdateJulia; update_julia()')
python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()"
julia -e 'using Pkg; Pkg.add(PackageSpec(name="ReactionMechanismSimulator",rev="main"));using ReactionMechanismSimulator'

if [ "$COMMAND_PKG" == "micromamba" ]; then
    micromamba activate t3_env
else
    conda activate t3_env
fi

# check that Python and Julia are being accessed from the t3_env
echo checking which python...
which python
echo checking which julia...
which julia
echo linking python-jl to python...
ln -sfn $(which python-jl) $(which python)

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
