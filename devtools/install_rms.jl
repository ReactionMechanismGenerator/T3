# Install RMS
import Pkg
Pkg.add(Pkg.PackageSpec(url="https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl", rev="fixcpkgimp"))
Pkg.build("ReactionMechanismSimulator")
using ReactionMechanismSimulator
Pkg.test("ReactionMechanismSimulator")
