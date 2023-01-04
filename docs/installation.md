# Installation

!!! Note
    T3 was only tested on Linux (Ubuntu_ 18.04.1 LTS) and Mac machines.
    We don't expect it to work on Windows.
    
    It can be installed on a server,
    as well as on your local desktop / laptop, submitting jobs to the server/s.


## Package Manager and Cloning from Github

### Unix-like platforms
- Install `curl` or `wget`

	```
	sudo apt install curl
	```
	or
	```
	sudo apt install wget
	```

- Install compiler
	Ubuntu or Debian
	```
	sudo apt install git gcc g++ make
	```
	Fedora or Red Hat
	```
	sudo dnf install git gcc gcc-c++ make
	```

- Download the installer using curl or wget or your favorite program download files and run the script.
 
For eg:

    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
    bash Mambaforge-$(uname)-$(uname -m).sh

or

    wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
    bash Mambaforge-$(uname)-$(uname -m).sh
- Get git if you don't have it already by typing
  <code>sudo apt-get install git</code> in a terminal.
- Clone T3's repository by typing the following command in the
  desired folder (e.g., under ~/Code/):

```console
$ git clone https://github.com/ReactionMechanismGenerator/T3.git

---> 100%
```

## Setting up Path

- To set up T3 to your local path in .bashrc, you have two options:

	- Option 1 - Terminal Command (NOTE: Make sure to change "~/Path/to/T3/" accordingly"):

		```echo 'PYTHONPATH=$PYTHONPATH:~/Path/to/T3/' >> ~/.bashrc ```

	- Option 2 - Editing .bashrc directly (NOTE: Make sure to change "~/Path/to/T3/" accordingly"):

		- In terminal, enter the command:
			``` gedit ~/.bashrc ```
		- Then in the opened file, on a new line, enter the following:
			``` export PYTHONPATH=$PYTHONPATH:~/Path/to/T3/ ```


## Install dependencies

- Install the latest versions of RMG-Py and the RMG-database on the same machine where 
  T3 is installed. Follow the instructions on
  <a href="http://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/installation/anacondaDeveloper.html">
  RMG's documentation</a>. Make sure to install RMG's **developer version**
  which includes important recent features. 
  Be sure to add RMG-Py to your PATH and PYTHONPATH as explained in the instructions.
- Install the latest version of ARC on the same machine where T3 is installed.
  Follow the instructions on
  <a href="https://reactionmechanismgenerator.github.io/ARC/installation.html">
  ARC's documentation</a>.
  Make sure to add ARC to your PATH and PYTHONPATH, as well as to define your servers as 
  explained in ARC's documentation.
- Create the Anaconda environment for T3 by executing the following command in the T3 folder:

<div class="termy">

```console
$ mamba env create -f environment.yml

INFO:     Collecting package metadata (repodata.json): done
INFO:     Solving environment: done
INFO:     
INFO:     Downloading and Extracting Packages
INFO:     ...
INFO:     Preparing transaction: done
INFO:     Verifying transaction: done
INFO:     Executing transaction: done
INFO:     #
INFO:     # To activate this environment, use
INFO:     #
INFO:     #     $ conda activate t3_env
INFO:     #
INFO:     # To deactivate an active environment, use
INFO:     #
INFO:     #     $ conda deactivate
```

</div>

- Activate the T3 environment every time before running T3:

    ``conda activate t3_env``


## Add T3 aliases to your .bashrc

Some optional yet convenient aliases are listed below
(make sure to change "/Path/to/T3/" accordingly).
Add these to your ``.bashrc`` file, which can be edited by typing, e.g., ``nano ~./bashrc``:
	
	export t3_path=$HOME`/Path/to/T3'
	alias t3e='source activate t3_env'
	alias t3='python $t3_path/T3.py input.yml'
	alias t3code='cd $t3_path'

Then, after sourcing ``.bashrc`` or restarting the terminal,
simply typing ``t3e`` will activate the environment,
typing ``t3code`` will change the directory into the T3 repository,
and finally typing ``t3`` from any folder with
a valid T3 input file will execute T3 in that folder.


## Updating T3

The T3 repository is being updated frequently.
Make sure to update your instance of T3 to enjoy new features
and get recent bug fixes. To get the most recent developer version,
execute the following commands (make sure to change ```~/Path/to/T3/``` accordingly):

	cd ~/Path/to/T3/
	git fetch origin
	git pull origin main

The above will update your **main** branch of T3.
