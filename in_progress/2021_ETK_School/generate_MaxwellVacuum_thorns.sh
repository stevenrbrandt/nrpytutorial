#!/bin/bash

echo "Generating MaxwellVacuum ETK thorn..."
./run_Jupyter_notebook.sh ETK_Workshop_2021-NRPy_tutorial.ipynb

echo "Generating MaxwellVacuumID ETK thorn..."
python MaxwellVacuumID.py
