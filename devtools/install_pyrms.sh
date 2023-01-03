# Properly configure the shell to use 'conda activate'.
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

conda activate t3_env
# check that Python and Julia are being accessed from the t3_env
echo checking which python...
which python
echo checking which julia...
which julia

echo linking python-jl to python...
ln -sfn $(which python-jl) $(which python)

#echo installing PyCall...
#julia devtools/install_pycall.jl

echo installing pyrms, RMS, and all required Julia packages...
python -c "import pyrms; pyrms.install(); import julia; julia.install(); import diffeqpy; diffeqpy.install()"
echo installing and building packages for Julia
julia -e 'using Pkg; Pkg.add(PackageSpec(name="ReactionMechanismSimulator",rev="main")); using ReactionMechanism Simulator'
conda deactivate
