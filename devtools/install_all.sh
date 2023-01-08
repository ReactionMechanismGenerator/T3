CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
###Check if Mamba is installed####
#Check if mamba/conda is installed
if [ -x "$(command -v mamba)" ]; then
	echo "mamba is installed."
	COMMAND_PKG=mamba
elif [ -x "$(command -v conda)" ]; then
	echo "conda is installed."
	COMMAND_PKG=conda
else
    echo "mamba and conda are not installed. Please download and install mamba or conda - we strongly recommend mamba"
fi

# Temporarily change directory to $HOME to install software
pushd .
cd ..

if [ -d "./RMG-database" ]; then
    cd RMG-database
    git pull https://github.com/ReactionMechanismGenerator/RMG-database
    cd ..
else
    git clone https://github.com/ReactionMechanismGenerator/RMG-database
fi;

if [ -d "./RMG-Py" ]; then
    cd RMG-Py
    git pull https://github.com/ReactionMechanismGenerator/RMG-Py
else
    git clone https://github.com/ReactionMechanismGenerator/RMG-Py
    cd RMG-Py
fi;


echo "export PYTHONPATH=$PYTHONPATH:$(pwd)" >> ~/.bashrc
echo "export PATH=$PATH:$(pwd)" >> ~/.bashrc 
export PYTHONPATH=$PYTHONPATH:$(pwd)
export PATH=$PATH:$(pwd)
source ~/.bashrc

#Creating rmg_env environment (or updating if it already exists)
if { conda env list | grep 'rmg_env'; } >/dev/null 2>&1; then
    $COMMAND_PKG env update -n rmg_env -f environment.yml
else
    $COMMAND_PKG env create -f environment.yml
fi;


#Prior to activation of the environment, the LD_LIBRARY_PATH needs to be set as an environmnet variable when rmg_env is activated.
#This exporting and unsetting will solve the RMS installation during the Julia compile and also when the environment is deactivated
#the original LD_LIBRARY_PATH is set.
echo "export OLD_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >> $CONDA_BASE/envs/rmg_env/etc/conda/activate.d/env_vars.sh
echo "export LD_LIBRARY_PATH=$CONDA_BASE/envs/rmg_env/lib:$LD_LIBRARY_PATH" >> $CONDA_BASE/envs/rmg_env/etc/conda/activate.d/env_vars.sh
echo "export LD_LIBRARY_PATH=${OLD_LD_LIBRARY_PATH}" >> $CONDA_BASE/envs/rmg_env/etc/conda/deactivate.d/env_vars.sh
echo "unset OLD_LD_LIBRARY_PATH" >> $CONDA_BASE/envs/rmg_env/etc/conda/deactivate.d/env_vars.sh


#Active rmg_env environment
conda activate rmg_env

#Compile RMG-Py
make
#Update pyjulia to the latest version
mamba update pyjulia -y

#Ensure that added paths etc. are set and then reactivate rmg_env
. ~/.bashrc
conda activate rmg_env

### Install python + julia connection and compile RMS
#
#This code here is for updating julia if ever required:     $(julia -e 'using Pkg; Pkg.add("UpdateJulia"); using UpdateJulia; update_julia()')
python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()"
julia -e 'using Pkg; Pkg.add(PackageSpec(name="ReactionMechanismSimulator",rev="main"));using ReactionMechanismSimulator'





conda activate t3_env
# check that Python and Julia are being accessed from the t3_env
echo checking which python...
which python
echo checking which julia...
which julia
echo linking python-jl to python...
ln -sfn $(which python-jl) $(which python)



cd ..
if [ -d "./ARC" ]; then
    cd ARC
    git pull https://github.com/ReactionMechanismGenerator/ARC
else
    git clone https://github.com/ReactionMechanismGenerator/ARC
    cd ARC
fi;

export PYTHONPATH=$PYTHONPATH:$(pwd)
echo "export PYTHONPATH=$PYTHONPATH:$(pwd)" >> ~/.bashrc
if { conda env list | grep 'arc_env'; } >/dev/null 2>&1; then
    $COMMAND_PKG env update -n arc_env -f environment.yml
else
    $COMMAND_PKG env create -f environment.yml
fi;
conda activate arc_env
#make install-all
###Temp### To be used until ARC repo install is properly updated
$COMMAND_PKG create -n xtb_env python=3.7 -y
conda activate xtb_env
$COMMAND_PKG install -c conda-forge xtb -y
$COMMAND_PKG install -c anaconda pyyaml -y
##########

conda deactivate

# Restore original directory
popd || exit
