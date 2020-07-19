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
echo $(pwd)
python ../T3/devtools/install_pyrms.py
echo "*****% source ~/.bashrc"
. ~/.bashrc
echo "*****% julia devtools/install_pyrms.jl"
julia ../T3/devtools/install_pyrms_1.jl
cd ../T3
conda env update -f environment.yml
cd ../pyrms
julia ../T3/devtools/install_pyrms_2.jl
# python3 -m pip install julia
echo "*****% python -c import julia; julia.install()"
python -c "import julia; julia.install()"
pip install diffeqpy
python -c "import diffeqpy; diffeqpy.install()"

ln -sfn $(which python-jl) $(which python)

# Restore original directory
popd || exit
