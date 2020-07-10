cd ..
git clone https://github.com/ReactionMechanismGenerator/pyrms
cd pyrms
git checkout py3
conda install -c conda-forge julia
julia -v
julia -e 'out = Pipe(); proc = run(pipeline(`which python`,stdout=out)); close(out.in); pypath = chomp(String(read(out))); ENV["CONDA_JL_HOME"] = join(split(pypath,"/")[1:end-2],"/"); ENV["PYTHON"] = pypath; using Pkg; Pkg.add("PyCall"); Pkg.build("PyCall");'
python devtools/install_pyrms.py
source ~/.bashrc
julia devtools/install_pyrms.jl
ln -sfn $(which python-jl) $(which python)
cd ../T3
