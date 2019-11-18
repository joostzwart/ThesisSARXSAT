# HySI-SAT 
Tool for identifying switched ARX and Piecewise ARX systems. Use at your own risk.
This tool was created as part of a thesis project at the TU Delft.
author: Joost Zwart

# Install instructions

## Linux

Run the following command to install:
>python3 setup.py build_ext --inplace

Run the following command to use the GUI:
>python3 gui.py

Run the following command after adding new code:
>python3 python3 setup_cython.py build_ext --inplace

# Dependencies
 * Numpy
 * Matplotlib
 * Pyqt5
 * Z3solve
 * Cplex (A license is required. )
 * Parma polyhedra library
 * pysqlite
 * sklearn
 
 # References
 The method used in this toolbox is based on the work of Ozay, Sznaier and Lagoa:
 N. Ozay, M. Sznaier, C. M. Lagoa and O. I. Camps, "A Sparsification Approach to Set Membership Identification of Switched Affine Systems," in IEEE Transactions on Automatic Control, vol. 57, no. 3, pp. 634-648, March 2012.
doi: 10.1109/TAC.2011.2166295,
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6004819&isnumber=6157057





