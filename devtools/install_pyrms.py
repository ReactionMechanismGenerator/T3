#!/usr/bin/env python3
# encoding: utf-8

"""
a script for installing pyRMS
"""

import os
import sys
from distutils.spawn import find_executable

print("appending pyrms to path")
original_dir = os.getcwd()
path = os.getcwd()
print(f'install_pyrms.py mssg 1: current path is {path}')
home_path = os.getenv("HOME")

ostype = sys.platform
if "darwin" in ostype:
    bash_path = os.path.join(home_path, ".bash_profile")
    print(f'install_pyrms.py mssg: darwin, OSX')
else:
    bash_path = os.path.join(home_path, ".bashrc")
    print(f'install_pyrms.py mssg: Linux')
st = f'\n#pyrms\nexport PYTHONPATH=$PYTHONPATH:{path}\nexport PYTHONPATH=$PYTHONPATH:{path}/pyrms'
with open(bash_path, 'a') as bash_file:
    bash_file.write(st)

julia_path = find_executable("julia")
print(f'install_pyrms.py mssg: julia_path = {julia_path}')

if not julia_path:
    julia_install_path = os.path.join(os.getenv("HOME"))
    print(f'install_pyrms.py mssg: not julia_path. julia_install_path = {julia_install_path}')

    if sys.platform and "darwin" in sys.platform:
        # Go Mac
        os.system( """curl -L https://julialang-s3.julialang.org/bin/mac/x64/1.3/julia-1.3.0-mac64.dmg -o "$HOME/Downloads/julia.dmg";""")
        os.system("""hdiutil attach ~/Downloads/julia.dmg;""")
        os.system(f"""cp -r /Volumes/Julia*/Julia*/Contents/Resources/julia {julia_install_path};""")
        os.system("""hdiutil detach -force /Volumes/Julia*;""")
        os.system(f"""mv {os.path.join(julia_install_path, "julia*")} {os.path.join(julia_install_path, "julia")}""")
    else:
        # Go Linux
        os.system("""curl -L https://julialang-s3.julialang.org/bin/linux/x64/1.3/julia-1.3.0-linux-x86_64.tar.gz -o "$HOME/Downloads/julia.tar.gz";""")
        os.system("""tar xzf "$HOME/Downloads/julia.tar.gz" -C "$HOME/Downloads";""")
        os.system(f"""cp -r "$(find "$HOME/Downloads" -maxdepth 2 -name "julia*" -type d | head -n 1)" "{julia_install_path}";""")
        os.system(f"""mv {os.path.join(julia_install_path, "julia*")} {os.path.join(julia_install_path, "julia")}""")

    print("appending julia to path julia")
    st = "\n#julia\nexport PATH=\"{0}:$PATH\"".format(os.path.join(home_path, "julia", "bin"))
    with open(bash_path, 'a') as bash_file:
        bash_file.write(st)
    os.system(f'export PATH=\"{os.path.join(home_path, "julia", "bin")}:$PATH\"')

os.chdir(original_dir)
