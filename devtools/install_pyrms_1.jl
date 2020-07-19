println("** install_pyrms_1.jl mssg: in install_pyrms_1.jl *****")

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
ENV["PYTHONHOME"] = pypath
println("* 6")
println(ENV["CONDA_JL_HOME"])
println(ENV["PYTHON"])
println("* 7")
Pkg.build("PyCall")
println("* 8")
# using PyCall
println("* 9")
