# Temporarily change directory to $HOME to install software
pushd .

cd ..
# RMG-database
git clone https://github.com/ReactionMechanismGenerator/RMG-database
# RMG-Py
git clone https://github.com/ReactionMechanismGenerator/RMG-Py
cd RMG-Py
export PYTHONPATH=$PYTHONPATH:$(pwd)
echo "export PYTHONPATH=$PYTHONPATH:$(pwd)" >> ~/.bashrc
# compile RMG
make

cd ..
# ARC
git clone https://github.com/ReactionMechanismGenerator/ARC
cd ARC
export PYTHONPATH=$PYTHONPATH:$(pwd)
echo "export PYTHONPATH=$PYTHONPATH:$(pwd)" >> ~/.bashrc

# Restore original directory
popd || exit
