# Installation

T3 was only tested on *Linux (Ubuntu_ 22.04.1 LTS)* and *MacOS*.
We don't expect it to work on Windows.

It can be installed on a server, as well as on your local desktop / laptop, in order to submit jobs to the server/s.

## Package Manager and Cloning from Github

### Unix-like platforms

#### 1. Install `curl` or `wget`

You can install either package.

##### *Option 1*

  ``` bash
   sudo apt install curl
  ```

##### *Option 2*

  ``` bash
  sudo apt install wget
  ```

#### 2. Install compiler

- Ubuntu or Debian

    ``` bash
    sudo apt install git gcc g++ make
    ```

- Fedora or Red Hat

    ``` bash
    sudo dnf install git gcc gcc-c++ make
    ```

#### 3. Download Python Package Manager

##### *Option 1*

``` bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
```

##### *Option 2*

``` bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
```

#### 4. Clone T3 Repository

Clone T3's repository by typing the following command in the desired folder (e.g., under ~/Code/):

``` bash
git clone https://github.com/ReactionMechanismGenerator/T3.git
```

## Setting up Path

- To set up T3 to your local path in .bashrc, you have two options:

### Option 1

Terminal Command (NOTE: Make sure to change "~/Path/to/T3/" accordingly"):

  ``` bash
  echo 'PYTHONPATH=$PYTHONPATH:~/Path/to/T3/' >> ~/.bashrc
  ```

### Option 2

Editing .bashrc directly (NOTE: Make sure to change "~/Path/to/T3/" accordingly"):

- In terminal, enter the command:

  ``` bash
  sudo gedit ~/.bashrc
  ```

- Then in the opened file, on a new line, enter the following:

  ``` text
  export PYTHONPATH=$PYTHONPATH:~/Path/to/T3/
  ```

## Install dependencies

T3 requires RMG-Py, RMG-databse and ARC to function correctly. In order to install the necessary dependencies, you can follow either option below.

### Option 1

- Navigate to the T3 folder, depending on where you cloned it to.

- Open a terminal in the T3 folder, and type the following:

  ``` bash
  make install all
  ```

 > **Note**: This can take some time to finish.

- You have now installed all the required dependencies.

### Option 2

- Install the latest versions of RMG-Py and the RMG-database on the same machine where T3 is installed. Follow the instructions on [RMG's Documentation](http://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/installation/anacondaDeveloper.html).
  - Make sure to install RMG's **developer version** which includes important recent features. Note that the installation instructions suggest Anaconda, but Mambaforge can be used in it's place.

  - Be sure to add RMG-Py to your PATH and PYTHONPATH as explained in the instructions.

- Install the latest version of ARC on the same machine where T3 is installed.
  Follow the instructions on [ARC's documentation](https://reactionmechanismgenerator.github.io/ARC/installation.html)

  - Make sure to add ARC to your PATH and PYTHONPATH, as well as to define your servers as explained in ARC's documentation.
  
- Create the Anaconda/Mambaforge environment for T3 by executing the following command in the T3 folder:

```bash
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

- Activate the T3 environment every time before running T3:

  ``` bash
  conda activate t3_env
  ```

## Add T3 aliases to your .bashrc

Some optional yet convenient aliases are listed below
(make sure to change "/Path/to/T3/" accordingly).
Add these to your ``.bashrc`` file, which can be edited by typing, e.g., ``sudo nano ~/.bashrc``:

```bash
export t3_path=$HOME`/Path/to/T3'
alias t3e='conda activate t3_env'
alias t3='python $t3_path/T3.py input.yml'
alias t3code='cd $t3_path'
```

Then, after sourcing ``.bashrc`` or restarting the terminal,
simply typing ``t3e`` will activate the environment,
typing ``t3code`` will change the directory into the T3 repository,
and finally typing ``t3`` from any folder with
a valid T3 input file will execute T3 in that folder.

## Updating T3

The T3 repository is being updated frequently.
Make sure to update your instance of T3 to enjoy new features
and get recent bug fixes. To get the most recent developer version,
execute the following commands (make sure to change ```~/Path/to/T3/``` accordingly:

``` bash
cd ~/Path/to/T3/
git fetch origin
git pull origin main
```

The above will update your **main** branch of T3.
