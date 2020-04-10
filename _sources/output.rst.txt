.. _output:

Output
======

If T3 is run from the terminal with the RMG-Py and ARC input files,
then the folder in which that file is located becomes the Project's folder.
After running a Project, the local Project folder will contain the following directory tree
(**bold** face represents folders, *italics* face represents files):

- *t3.log*: Details of all project execution procedures.
- **iteration_X**: Contains the output from RMG-Py and from ARC. Please see the `documentation for RMG-Py section 1.7 <https://github.com/ReactionMechanismGenerator/RMG-Py/blob/master/documentation/RMG-Py_and_Arkane_Documentation.pdf>`_ and the `documentation for ARC's output <https://reactionmechanismgenerator.github.io/ARC/output.html>`_ to learn more about their respective outputs and file structure. Here, X ranges from 0 to the number of maximum tandem iterations; the default is a maximum of 10 iterations as described in `tandem.main`__, though the user can specify a different value by adding the following argument to their ``arc.yml`` file::

    max tandem iterations: X 


__ api/main.html

.. include:: links.txt
