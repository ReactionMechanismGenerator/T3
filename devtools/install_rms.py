import os
import sys

print("appending pyrms to path")
os.chdir('../../pyrms')
path=os.getcwd()
homepath = os.getenv("HOME")

ostype = sys.platform
if "darwin" in ostype:
    bpath = os.path.join(homepath,".bash_profile")
else:
    bpath = os.path.join(homepath,".bashrc")
st = "\n#pyrms\nexport PYTHONPATH=$PYTHONPATH:{0}".format(path)
with open(bpath, 'a') as bfile:
    bfile.write(st)
    