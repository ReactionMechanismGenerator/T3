using Pkg
Pkg.add("PyCall")

out = Pipe()
proc = run(pipeline(`which python`,stdout=out))
close(out.in)
pypath = chomp(String(read(out)))

# set env variables for installing PyCall
ENV["CONDA_JL_HOME"] = join(split(pypath, "/")[1:end-2], "/")
ENV["PYTHON"] = pypath
ENV["PYTHONHOME"] = pypath

Pkg.build("PyCall")
using PyCall
Pkg.add("DifferentialEquations")
using DifferentialEquations
Pkg.add(Pkg.PackageSpec(url="https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl", rev="fixcpkgimp"))
Pkg.build("ReactionMechanismSimulator")
using ReactionMechanismSimulator
Pkg.test("ReactionMechanismSimulator")
