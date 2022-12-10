#!/bin/bash

jupyter nbconvert --to python $1 --output=blah
ipython blah.py && rm -f blah.py
