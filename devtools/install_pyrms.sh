pushd .
cd ..
git clone https://github.com/ReactionMechanismGenerator/pyrms
cd pyrms
git checkout py3
julia -v
#echo "about to install J"
#julia -e 'out = Pipe(); proc = run(pipeline(`which python`,stdout=out)); close(out.in); pypath = chomp(String(read(out))); ENV["CONDA_JL_HOME"] = join(split(pypath,"/")[1:end-2], "/"); ENV["PYTHON"] = pypath; using Pkg; Pkg.add("PyCall"); Pkg.build("PyCall");'
#echo "DONE installing J"
echo "*****% python devtools/install_pyrms.py"
cd ../T3
python devtools/install_pyrms.py
echo "*****% source ~/.bashrc"
source ~/.bashrc
echo "*****% julia devtools/install_pyrms.jl"
julia devtools/install_pyrms.jl
# python3 -m pip install julia
echo "*****% python -c import julia; julia.install()
python -c "import julia; julia.install()"
pip install diffeqpy
python -c "import diffeqpy; diffeqpy.install()"

ln -sfn $(which python-jl) $(which python)

# Restore original directory
popd
