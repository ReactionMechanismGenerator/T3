echo "***** in install_pyrms.sh"
pushd .
cd ..
echo "***** clone pyrms"
git clone https://github.com/ReactionMechanismGenerator/pyrms
cd pyrms
echo "checkout pyrms py3 branch"
git checkout py3
echo "julia version"
julia -v
echo "call install_pyrms.py"
python ../T3/devtools/install_pyrms.py
echo "***** source ~/.bashrc"
. ~/.bashrc
echo "***** re-activate t3_env"
source activate t3_env
echo "***** call install_pyrms_1.jl"
julia ../T3/devtools/install_pyrms_1.jl
cd ../T3
echo "***** update new environment"
conda env update -f environment.yml
cd ../pyrms
echo "***** call install_pyrms_1.jl"
julia ../T3/devtools/install_pyrms_2.jl
# install python-jl (simply called julia)
echo "***** install python-jl"
python3 -m pip install julia
echo "***** python -c import julia; julia.install()"
python -c "import julia; julia.install()"
echo "***** pip install diffeqpy"
pip install diffeqpy
echo "***** python -c import diffeqpy; diffeqpy.install()"
python -c "import diffeqpy; diffeqpy.install()"

# uncomment if must
echo "***** lionk python and python-jl"
ln -sfn $(which python-jl) $(which python)

# Restore original directory
popd || exit
