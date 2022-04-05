#!/bin/bash

jupyter nbconvert --to python $1 --output=blah
python3 blah.py
rm -f blah.py
