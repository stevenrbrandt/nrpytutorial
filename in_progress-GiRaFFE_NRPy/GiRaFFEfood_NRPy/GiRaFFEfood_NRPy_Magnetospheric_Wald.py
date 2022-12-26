# Step 0: Import the NRPy+ core modules and set the reference metric to Cartesian
# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)
nrpy_dir_path = os.path.join("../..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)
giraffefood_dir_path = os.path.join("GiRaFFEfood_NRPy")
if giraffefood_dir_path not in sys.path:
    sys.path.append(giraffefood_dir_path)

# Step 0.a: Import the NRPy+ core modules and set the reference metric to Cartesian
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import GiRaFFEfood_NRPy_Common_Functions as gfcf # Some useful functions for GiRaFFE initial data.

import reference_metric as rfm
par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()

# Step 1a: Set commonly used parameters.
thismodule = __name__

C0 = par.Cparameters("REAL",thismodule,"C0",1.0)

# Step 2: Set the vectors A and E in Spherical coordinates
def Ar_MW(r,theta,phi, **params):
    g4DD = params["g4DD"]
    a    = params["a"]
    return sp.Rational(1,2) * C0 * (g4DD[1][3] + 2 * a * g4DD[0][1])

def Ath_MW(r,theta,phi, **params):
    g4DD = params["g4DD"]
    a    = params["a"]
    return sp.Rational(1,2) * C0 * (g4DD[2][3] + 2 * a * g4DD[0][2])

def Aph_MW(r,theta,phi, **params):
    g4DD = params["g4DD"]
    a    = params["a"]
    return sp.Rational(1,2) * C0 * (g4DD[3][3] + 2 * a * g4DD[0][3])

#Step 3: Compute v^i from A_i and E_i
def ValenciavU_func_MW(**params):
    return ixp.zerorank1()