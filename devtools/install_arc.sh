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

# Restore original directory
popd || exit
