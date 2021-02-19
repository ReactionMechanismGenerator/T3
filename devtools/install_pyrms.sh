# check that Python and Julia are being accessed from the t3_env
echo checking which python...
which python
echo checking which julia...
which julia

echo linking python-jl to python...
ln -sfn $(which python-jl) $(which python)

echo installing PyCall...
julia devtools/install_pycall.jl

echo installing pyrms, RMS, and all required Julia packages...
python -c "import pyrms; pyrms.install()"

# the above line installs RMS v0.2.1
# v0.3.2 is needed to use the latest SA features
# once there is a functional package for v0.3.2, 
# the line below can be deleted since pyrms.install() will pull the latest version
echo installing RMS v0.3.2...
julia devtools/install_RMS_v032.jl
