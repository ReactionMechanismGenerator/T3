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

#Prior to activation of the environment, the LD_LIBRARY_PATH needs to be set as an environmnet variable when rmg_env is activated.
#This exporting and unsetting will solve the RMS installation during the Julia compile and also when the environment is deactivated
#the original LD_LIBRARY_PATH is set.
echo "export OLD_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >> "$BASE/envs/rmg_env/etc/conda/activate.d/env_vars.sh"
echo "export LD_LIBRARY_PATH=$BASE/envs/rmg_env/lib:$LD_LIBRARY_PATH" >> "$BASE/envs/rmg_env/etc/conda/activate.d/env_vars.sh"
echo "export LD_LIBRARY_PATH=${OLD_LD_LIBRARY_PATH}" >> "$BASE/envs/rmg_env/etc/conda/deactivate.d/env_vars.sh"
echo "unset OLD_LD_LIBRARY_PATH" >> "$BASE/envs/rmg_env/etc/conda/deactivate.d/env_vars.sh"

#Active rmg_env environment
if [ "$COMMAND_PKG" == "micromamba" ]; then
    micromamba activate rmg_env
else
    conda activate rmg_env
fi

#Ensure that added paths etc. are set and then reactivate rmg_env
. ~/.bashrc
if [ "$COMMAND_PKG" == "micromamba" ]; then
    micromamba activate rmg_env
else
    conda activate rmg_env
fi



### Install python + julia connection and compile RMS
echo linking python-jl to python...
ln -sfn $(which python-jl) $(which python)
#
#This code here is for updating julia if ever required:     $(julia -e 'using Pkg; Pkg.add("UpdateJulia"); using UpdateJulia; update_julia()')
python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()"
julia -e 'using Pkg; Pkg.add(PackageSpec(name="ReactionMechanismSimulator",rev="main"));using ReactionMechanismSimulator'
