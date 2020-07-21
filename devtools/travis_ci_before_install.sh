# Store original directory
pushd .

# Temporarily change directory to $HOME to install software
cd $HOME

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
echo "export PATH=$MINICONDA_HOME/bin:$PATH" >> ~/.bashrc
#echo "export -f conda" >> ~/.bashrc
#echo "export -f __conda_activate" >> ~/.bashrc
#echo "export -f __conda_reactivate" >> ~/.bashrc
#echo "export -f __conda_hashr" >> ~/.bashrc
#echo "export -f __add_sys_prefix_to_path" >> ~/.bashrc

conda config --set always_yes yes --set changeps1 no
conda update --q conda
conda info -a

# Restore original directory
popd || exit
