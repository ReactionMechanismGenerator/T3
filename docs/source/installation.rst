.. _installation:

Installation instructions
=========================

Note: T3 was only tested on Linux (Ubuntu_ 18.04.1 LTS) and Mac machines. We don't expect it to work on Windows.
T3 can be installed on a server, as well as on your local desktop / laptop, submitting jobs to the server/s.
Note that T3 requires installation of ARC_, Automated Rate calculator.


.. _path:

Clone and setup path
^^^^^^^^^^^^^^^^^^^^

- Download and install the `Anaconda Python Platform`__ for Python 3.7 or higher if you haven't already.
- Get git if you don't have it already by typing sudo apt-get install git in a terminal.
- Clone T3's repository by typing the following command in the desired folder (e.g., under `~/home/Code/`)::

    git clone https://github.com/ReactionMechanismGenerator/T3.git
			  
- Add T3 to your local path in .bashrc (make sure to change `~/Path/to/T3/` accordingly)::

    export PYTHONPATH=$PYTHONPATH:~/Path/to/T3/

__ anaconda_


Install dependencies
^^^^^^^^^^^^^^^^^^^^

- Install the latest versions of RMG-Py and the RMG-database on the same machine that 
  T3 is installed on by following the `documented instructions <http://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/installation/anacondaDeveloper.html>`_. 
  Be sure to add RMG-Py to your PATH and PYTHONPATH as explained in the instructions.
- Install the latest version of ARC on the same machine that T3 is installed on 
  by following `ARC's installation instructions 
  <https://reactionmechanismgenerator.github.io/ARC/installation.html>`_. 
  Make sure to add ARC to your PATH and PYTHONPATH variables as well as define servers as 
  explained in ARC's documentation. 
- If you'd like to use `AutoTST <https://github.com/ReactionMechanismGenerator/AutoTST>`_ in ARC (optional),
  clone it in a separate folder and add it to your PYTHONPATH as well.
- Create the Anaconda environment for T3::

    conda env create -f environment.yml

  Activate the T3 environment every time before you run T3::

     source activate t3_env

.. _aliases:

Add T3 aliases to your .bashrc (for convenience)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Some optional, but convenient, aliases are listed below (make sure to change `~/Path/to/T3/` accordingly).
Add these to your ``.bashrc`` file, which can be edited by typing, e.g., ``nano ~./bashrc``::
	
	export t3_path=$HOME`~/Path/to/T3'
	alias t3e='source activate t3_env'
	alias t3='python $t3_path/t3.py -r input.py -a input.yml'
	alias t3code='cd $t3_path'


Updating T3
^^^^^^^^^^^

T3 is being updated frequently. Make sure to update T3 and enjoy new features and bug fixes.
To get the most recent developer version, execute the following commands (make sure to change `~/Path/to/T3/` accordingly)::

	cd ~/Path/to/T3/
	git stash
	git fetch origin
	git pull origin master
	git stash pop

The above will update your `master` branch of T3.

**Note:** This process might cause merge conflicts if the updated version (either the developer version
or a stable version) changes a file you changed locally. Although we try to avoid causing merge conflicts
for T3's users as much as we can, it could still sometimes happen.
You'll identify a merge conflict if git prints a message similar to this::

    $ git merge BRANCH-NAME
    > Auto-merging settings.py
    > CONFLICT (content): Merge conflict in styleguide.md
    > Automatic merge failed; fix conflicts and then commit the result

Detailed steps to resolve a git merge conflict can be found `online`__.

__ mergeConflict_

Principally, you should open the files that have merge conflicts, and look for the following markings::

    <<<<<<< HEAD
    this is some content introduced by updating ARC
    =======
    totally different content the user added, adding different changes
    to the same lines that were also updated remotely
    >>>>>>> new_branch_to_merge_later

Resolving a merge conflict consists of three stages:

- determine which version of the code you'd like to keep
  (usually you should manually append your oun changes to the more
  updated ARC code). Make the changes and get rid of the unneeded ``<<<<<<< HEAD``,
  ``=======``, and ``>>>>>>> new_branch_to_merge_later`` markings. Repeat for all conflicts.
- Stage the changed by typing: ``git add .``
- If you don't plan to commit your changes, unstage them by typing: ``git reset --soft origin/master``


.. include:: links.txt
