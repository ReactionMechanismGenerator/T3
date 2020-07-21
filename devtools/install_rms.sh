# Store original directory
pushd .

# Create an empty environment with the same name (t3_env) before installing RMS
conda deactivate
conda remove --name t3_env --all -y
conda create -n t3_env python=3.7
source activate t3_env

# Install Julia
cd $HOME
mkdir -p ~/temp_Downloads

if [ "$TRAVIS_OS_NAME" == "linux" ]; then curl -L https://julialang-s3.julialang.org/bin/linux/x64/1.4/julia-1.4.2-linux-x86_64.tar.gz -o "$HOME/temp_Downloads/julia.tar.gz"; fi
if [ "$TRAVIS_OS_NAME" == "linux" ]; then tar xzf "$HOME/temp_Downloads/julia.tar.gz" -C "$HOME/temp_Downloads"; fi
if [ "$TRAVIS_OS_NAME" == "linux" ]; then cp -r "$(find "$HOME/temp_Downloads" -maxdepth 2 -name "julia*" -type d | head -n 1)" "$HOME/julia"; fi

if [ "$TRAVIS_OS_NAME" == "osx" ]; then curl -L https://julialang-s3.julialang.org/bin/mac/x64/1.4/julia-1.4.2-mac64.dmg -o "$HOME/temp_Downloads/julia.dmg"; fi
if [ "$TRAVIS_OS_NAME" == "osx" ]; then hdiutil attach ~/temp_Downloads/julia.dmg; fi
if [ "$TRAVIS_OS_NAME" == "osx" ]; then cp -r /Volumes/Julia*/Julia*/Contents/Resources/julia $HOME/julia; fi
if [ "$TRAVIS_OS_NAME" == "osx" ]; then hdiutil detach -force /Volumes/Julia*; fi

rm -rf ~/temp_Downloads/julia*
export PATH="$HOME/julia/bin:$PATH"
echo "export PATH=$HOME/julia/bin:$PATH" >> ~/.bashrc

echo "Julia version:"
julia -v

# Restore original directory
popd || exit

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
echo 'export "PYTHONPATH=$PYTHONPATH:$(pwd)"' >> ~/.bashrc
export PYTHONPATH=$PYTHONPATH:$(pwd)/pyrms
echo 'export "PYTHONPATH=$PYTHONPATH:$(pwd)/pyrms"' >> ~/.bashrc
echo "***** python path is:"
echo $PYTHONPATH
# echo "***** source ~/.bashrc"
# . ~/.bashrc
#echo "***** re-activate t3_env"
#source activate t3_env
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

echo "***** link python and python-jl"
ln -sfn $(which python-jl) $(which python)

# Restore original directory
popd || exit
