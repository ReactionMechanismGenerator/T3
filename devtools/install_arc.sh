# Properly configure the shell to use 'conda activate'.
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


export PYTHONPATH=$PYTHONPATH:$(pwd)
echo "export PYTHONPATH=$PYTHONPATH:$(pwd)" >> ~/.bashrc


$COMMAND_PKG env update -n rmg_env -f environment.yml
conda activate rmg_env
make
python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()"
julia -e 'using Pkg; Pkg.add(PackageSpec(name="ReactionMechanismSimulator",rev="main")); using ReactionMechanismSimulator;'
conda deactivate

cd ..
git clone https://github.com/ReactionMechanismGenerator/ARC
cd ARC
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
