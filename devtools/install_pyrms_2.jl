println("** install_pyrms_2.jl mssg: in install_pyrms_2.jl *****")

println("using Pkg *************")
using Pkg
println("adding RMS *************")
Pkg.add("DifferentialEquations")
println("* 9")
using DifferentialEquations
println("* 10")
Pkg.add(Pkg.PackageSpec(url="https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl", rev="fixcpkgimp"))
println("* 11")
Pkg.build("ReactionMechanismSimulator")
println("* 12")
println("importing PyCall and RMS **********")
using PyCall
println("* 13")
using ReactionMechanismSimulator
println("* 14 run RMS test")
Pkg.test("ReactionMechanismSimulator")