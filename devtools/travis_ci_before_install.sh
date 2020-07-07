# Temporarily change directory to $HOME to install software
pushd .
cd $HOME

# make temporary directory
mkdir -p ~/temp_Downloads

# Install Julia
if [ "$TRAVIS_OS_NAME" == "linux" ]; then curl -L https://julialang-s3.julialang.org/bin/linux/x64/1.1/julia-1.1.0-linux-x86_64.tar.gz -o "$HOME/temp_Downloads/julia.tar.gz"; fi
if [ "$TRAVIS_OS_NAME" == "linux" ]; then tar xzf "$HOME/temp_Downloads/julia.tar.gz" -C "$HOME/temp_Downloads"; fi
if [ "$TRAVIS_OS_NAME" == "linux" ]; then cp -r "$(find "$HOME/temp_Downloads" -maxdepth 2 -name "julia*" -type d | head -n 1)" "$HOME/julia"; fi

if [ "$TRAVIS_OS_NAME" == "osx" ]; then curl -L https://julialang-s3.julialang.org/bin/mac/x64/1.1/julia-1.1.0-mac64.dmg -o "$HOME/temp_Downloads/julia.dmg"; fi
if [ "$TRAVIS_OS_NAME" == "osx" ]; then hdiutil attach ~/temp_Downloads/julia.dmg; fi
if [ "$TRAVIS_OS_NAME" == "osx" ]; then cp -r /Volumes/Julia*/Julia*/Contents/Resources/julia $HOME/julia; fi
if [ "$TRAVIS_OS_NAME" == "osx" ]; then hdiutil detach -force /Volumes/Julia*; fi

rm -rf ~/temp_Downloads/julia*
export PATH="$HOME/julia/bin:$PATH"

# Install Miniconda
if [ "$TRAVIS_OS_NAME" == "osx" ]; then
    MINICONDA=Miniconda3-latest-MacOSX-x86_64.sh
else
    MINICONDA=Miniconda3-latest-Linux-x86_64.sh
fi
MINICONDA_HOME=$HOME/miniconda
wget -q https://repo.continuum.io/miniconda/$MINICONDA
bash $MINICONDA -b -p $MINICONDA_HOME

# Configure miniconda
export PIP_ARGS="-U"
export PATH=$MINICONDA_HOME/bin:$PATH
    
conda config --set always_yes yes --set changeps1 no
conda update --q conda
conda info -a

# Restore original directory
popd
