# Store original directory
pushd .

python -c "import pyrms; pyrms.install()"

# Change python-jl to call python3:
if [ "$TRAVIS_OS_NAME" == "linux" ]; then sed -i 's|bin/python|bin/python3|g' $(which python-jl); fi
if [ "$TRAVIS_OS_NAME" == "osx" ]; then sed -i ".bak" 's|bin/python|bin/python3|g' $(which python-jl); fi

ln -sfn $(which python-jl) $(which python)

echo $PYTHONPATH

# Restore original directory
popd || exit
