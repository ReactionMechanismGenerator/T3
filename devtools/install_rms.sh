# Store original directory
pushd .

cd ../
echo "***** clone pyrms: cloning into:"
echo "$(pwd)"
git clone https://github.com/ReactionMechanismGenerator/pyrms
cd pyrms
echo "checkout pyrms py3 branch"
git checkout py3

export PYTHONPATH=$PYTHONPATH:$(pwd)
echo "export PYTHONPATH=$PYTHONPATH:$(pwd)" >> ~/.bashrc
export PYTHONPATH=$PYTHONPATH:$(pwd)/pyrms
echo "export PYTHONPATH=$PYTHONPATH:$(pwd)/pyrms" >> ~/.bashrc
echo "***** python path 1 is:"
echo $PYTHONPATH
echo $PATH
cd ../ARC
export PYTHONPATH=$PYTHONPATH:$(pwd)
echo "export PYTHONPATH=$PYTHONPATH:$(pwd)" >> ~/.bashrc
cd ../RMG-Py
export PYTHONPATH=$PYTHONPATH:$(pwd)
echo "export PYTHONPATH=$PYTHONPATH:$(pwd)" >> ~/.bashrc
export PATH="$HOME/julia/bin:$PATH"
echo "export PATH=$HOME/julia/bin:$PATH" >> ~/.bashrc
cd ../pyrms
echo "***** python path 1.5 is:"
echo $PYTHONPATH
echo $PATH
# echo "***** source ~/.bashrc"
# . ~/.bashrc
#echo "***** re-activate t3_env"
#source activate t3_env
echo "***** python path 2 is:"
echo $PYTHONPATH
echo "***** call install_pyrms_1.jl"
julia ../T3/devtools/install_pyrms_1.jl
cd ../T3
echo "***** update new environment"
conda env update -f environment.yml
echo "***** python path 3 is:"
echo $PYTHONPATH
cd ../pyrms
echo "***** call install_pyrms_1.jl"
julia ../T3/devtools/install_pyrms_2.jl
# install python-jl (simply called julia)
echo "***** python path 4 is:"
echo $PYTHONPATH
echo "***** install python-jl"
python3 -m pip install julia
echo "***** python -c import julia; julia.install()"
python -c "import julia; julia.install()"
echo "***** python path 5 is:"
echo $PYTHONPATH
echo "***** pip install diffeqpy"
pip install diffeqpy
echo "***** python -c import diffeqpy; diffeqpy.install()"
python -c "import diffeqpy; diffeqpy.install()"
echo "python -c 'import diffeqpy; diffeqpy.install()'" >> ~/.bashrc
echo "***** python path 6 is:"
echo $PYTHONPATH

echo "***** link python and python-jl"
ln -sfn $(which python-jl) $(which python)
echo "***** python path 7 is:"
echo $PYTHONPATH

# Restore original directory
popd || exit
