using Pkg

# get the python path
link_python = false
if !("PyCall" in keys(Pkg.installed()))
    Pkg.add("PyCall")
    println("** link_python 1 *****")
    link_python = true
else
    println("found PyCall importing")
    v = 0
    try
        using PyCall
        sys = pyimport("sys")
        v = sys.version_info[1]
    catch
        v = -1  # PyCall isn't setup right
    end
    if v != 3
        println("removing PyCall and reinstalling")
        Pkg.rm("PyCall")
        if "Conda" in keys(Pkg.installed())
            Pkg.rm("Conda")
        end
        Pkg.add("PyCall")
    println("** link_python 2 *****")
        link_python = true
    end
end

println("PyCall present")
if link_python
    println("Linking PyCall properly")
    out = Pipe()
    proc = run(pipeline(`which python`,stdout=out))
    close(out.in)
    pypath = chomp(String(read(out)))

    # set env variables for installing PyCall
    ENV["CONDA_JL_HOME"] = join(split(pypath, "/")[1:end-2], "/")
    ENV["PYTHON"] = pypath
    println(ENV["CONDA_JL_HOME"])
    println(ENV["PYTHON"])
    Pkg.build("PyCall")
end

println("adding RMS")
Pkg.add("DifferentialEquations")
using DifferentialEquations
Pkg.add(PackageSpec(url="https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl", rev="master"))
Pkg.build("ReactionMechanismSimulator")
println("importing PyCall and RMS")
using PyCall
using ReactionMechanismSimulator
