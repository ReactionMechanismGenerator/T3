println("** install_pyrms.jl mssg: in install_pyrms.jl *****")

using Pkg
println("** install_pyrms.jl mssg: in install_pyrms.jl after using Pkg *****")

Pkg.add("PyCall")
println("** link_python 1 *****")
link_python = true

println("PyCall present  *************")
println("Linking PyCall properly *************")
out = Pipe()
println("* 2")
proc = run(pipeline(`which python`,stdout=out))
println("* 3")
close(out.in)
println("* 4")
pypath = chomp(String(read(out)))
println("* 5")

# set env variables for installing PyCall
ENV["CONDA_JL_HOME"] = join(split(pypath, "/")[1:end-2], "/")
ENV["PYTHON"] = pypath
println("* 6")
println(ENV["CONDA_JL_HOME"])
println(ENV["PYTHON"])
println("* 7.1")
Pkg.add("Conda")
println("* 7.2")
Pkg.build("Conda")
println("* 7.3")
using Conda
println("* 7.4")
Pkg.build("PyCall")
println("* 8")

println("adding RMS *************")
Pkg.add("DifferentialEquations")
println("* 9")
using DifferentialEquations
println("* 10")
Pkg.add(PackageSpec(url="https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl", rev="master"))
println("* 11")
Pkg.build("ReactionMechanismSimulator")
println("* 12")
println("importing PyCall and RMS **********")
using PyCall
println("* 13")
using ReactionMechanismSimulator
println("* 14")
