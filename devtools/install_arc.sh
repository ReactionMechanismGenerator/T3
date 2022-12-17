# Properly configure the shell to use 'conda activate'.
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# Temporarily change directory to $HOME to install software
pushd .

cd ..
git clone https://github.com/ReactionMechanismGenerator/RMG-database
git clone https://github.com/ReactionMechanismGenerator/RMG-Py
cd RMG-Py
export PYTHONPATH=$PYTHONPATH:$(pwd)
echo "export PYTHONPATH=$PYTHONPATH:$(pwd)" >> ~/.bashrc
mamba env update -n rmg_env -f environment.yml
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
    mamba env update -n arc_env -f environment.yml
else
    mamba env create -f environment.yml
fi;
mamba activate arc_env
#make install-all
###Temp###
mamba create -n xtb_env python=3.7 -y
conda activate xtb_env
mamba install -c conda-forge xtb -y
mamba install -c anaconda pyyaml -y
##########

conda deactivate

# Restore original directory
popd || exit
