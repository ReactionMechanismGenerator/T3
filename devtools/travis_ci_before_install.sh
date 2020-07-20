# Temporarily change directory to $HOME to install software
pushd .
cd $HOME

# make a temporary directory
mkdir -p ~/temp_Downloads

# Install Miniconda
MINICONDA_HOME=$HOME/miniconda
wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash $MINICONDA -b -p $MINICONDA_HOME

# Configure miniconda
export PIP_ARGS="-U"
export PATH=$MINICONDA_HOME/bin:$PATH
echo 'export "PATH=$MINICONDA_HOME/bin:$PATH"' >> ~/.bashrc
    
conda config --set always_yes yes --set changeps1 no
conda update --q conda
conda info -a
conda init bash

# Restore original directory
popd || exit
