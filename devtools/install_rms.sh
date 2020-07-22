# Store original directory
pushd .

echo "***** call install_pyrms_1.jl"
julia devtools/install_pyrms_1.jl
echo "***** call install_pyrms_2.jl"
julia devtools/install_pyrms_2.jl
# install python-jl (simply called julia)
echo "***** install python-jl"
python3 -m pip install julia
echo "***** python -c import julia; julia.install()"
python -c "import julia; julia.install()"
#echo "***** pip install diffeqpy"
#pip install diffeqpy
#echo "***** python -c import diffeqpy; diffeqpy.install()"
#python -c "import diffeqpy; diffeqpy.install()"
#python -c "import Pkg; Pkg.add('PyCall')"

# Change python-jl to call python3:
if [ "$TRAVIS_OS_NAME" == "linux" ]; then sed -i 's|bin/python|bin/python3|g' $(which python-jl); fi
if [ "$TRAVIS_OS_NAME" == "osx" ]; then sed -i ".bak" 's|bin/python|bin/python3|g' $(which python-jl); fi

echo "***** link python and python-jl"
ln -sfn $(which python-jl) $(which python)
echo "***** python path is:"
echo $PYTHONPATH

# Restore original directory
popd || exit
