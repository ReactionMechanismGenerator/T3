.. _running:

Running T3
==========

Using an input file
^^^^^^^^^^^^^^^^^^^

To run T3, make sure to first activate the T3 environment
(see the :ref:`recommended aliases <aliases>` from the installation instructions).
Then simply type::

    python <path_to_the_T3_folder>/t3.py -r input.py -a input.yml

replacing ``<path_to_the_T3_folder>`` in the above command with your actual local path to T3.
You could of course name the input file whatever you like.
However, if you're using the :ref:`recommended aliases <aliases>`, then simply typing::

    t3

in any folder with a valid RMG-Py ``input.py`` file and ARC ``input.yml`` file will execute T3 using those files.
Please refer to the `documentation for RMG-Py section 1.4 <https://github.com/ReactionMechanismGenerator/RMG-Py/blob/master/documentation/RMG-Py_and_Arkane_Documentation.pdf>`_ and the `documentation for using an input file in ARC <https://reactionmechanismgenerator.github.io/ARC/running.html>`_ to learn more about their respective input file formats and options.


Using the API
^^^^^^^^^^^^^

To run T3, make sure to first activate the T3 environment
(see the :ref:`recommended aliases <aliases>` from the installation instructions).

T3's API can be reached from any python platform, if T3 was added to the PYTHONPATH
(see the :ref:`installation instructions <path>`).

.. include:: links.txt
