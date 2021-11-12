# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)
nrpy_dir_path = os.path.join("../..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

# Step 0.a: Import the NRPy+ core modules and set the reference metric to Cartesian
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import GiRaFFEfood_NRPy.GiRaFFEfood_NRPy_Common_Functions as gfcf # Some useful functions for GiRaFFE initial data.

import reference_metric as rfm   # NRPy+: Reference metric support
par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()

# Step 1a: Set commonly used parameters.
thismodule = __name__

import Min_Max_and_Piecewise_Expressions as noif

def Ax_TW(x,y,z, **params):
    return sp.sympify(0)

def Ay_TW(x,y,z, **params):
    return sp.Rational(7,2)*x*noif.coord_greater_bound(-x,0) + sp.sympify(3)*x*noif.coord_greater_bound(x,0)

def Az_TW(x,y,z, **params):
    return y-sp.Rational(3,2)*x*noif.coord_greater_bound(-x,0) - sp.sympify(3)*x*noif.coord_greater_bound(x,0)

#Step 3: Compute v^i from B^i and E_i
def ValenciavU_func_TW(**params):
    x = rfm.xx_to_Cart[0]

    B_aU = ixp.zerorank1(DIM=3)
    E_aU = ixp.zerorank1(DIM=3)
    B_pU = ixp.zerorank1(DIM=3)
    E_pU = ixp.zerorank1(DIM=3)
    B_mU = ixp.zerorank1(DIM=3)
    E_mU = ixp.zerorank1(DIM=3)

    B_aU[0] = sp.sympify(1)
    B_aU[1] = noif.coord_leq_bound(x,0) * sp.sympify(1) + noif.coord_greater_bound(x,0) * sp.Rational(3,2)
    B_aU[2] = sp.sympify(2)

    E_aU[0] = noif.coord_leq_bound(x,0) * sp.sympify(-1) + noif.coord_greater_bound(x,0) * sp.Rational(-3,2)
    E_aU[1] = sp.sympify(1)
    E_aU[2] = sp.sympify(0)

    B_pU[0] = sp.sympify(0)
    B_pU[1] = noif.coord_leq_bound(x,0) * sp.sympify(0) + noif.coord_greater_bound(x,0) * sp.Rational(3,2)
    B_pU[2] = noif.coord_leq_bound(x,0) * sp.sympify(0) + noif.coord_greater_bound(x,0) * sp.sympify(1)

    E_pU[0] = sp.sympify(0)
    E_pU[1] = noif.coord_leq_bound(x,0) * sp.sympify(0) + noif.coord_greater_bound(x,0) * sp.sympify(1)
    E_pU[2] = noif.coord_leq_bound(x,0) * sp.sympify(0) + noif.coord_greater_bound(x,0) * sp.Rational(-3,2)

    B_mU[0] = sp.sympify(0)
    B_mU[1] = noif.coord_leq_bound(x,0) * sp.Rational(1,2)  + noif.coord_greater_bound(x,0) * sp.sympify(0)
    B_mU[2] = noif.coord_leq_bound(x,0) * sp.Rational(3,2) + noif.coord_greater_bound(x,0) * sp.sympify(0)

    E_mU[0] = sp.sympify(0)
    E_mU[1] = noif.coord_leq_bound(x,0) * sp.Rational(-3,2) + noif.coord_greater_bound(x,0) * sp.sympify(0)
    E_mU[2] = noif.coord_leq_bound(x,0) * sp.Rational(1,2)  + noif.coord_greater_bound(x,0) * sp.sympify(0)

    BU = ixp.zerorank1(DIM=3)
    EU = ixp.zerorank1(DIM=3)
    for i in range(3):
        BU[i] = B_aU[i] + B_pU[i] + B_mU[i]
        EU[i] = E_aU[i] + E_pU[i] + E_mU[i]

    # In flat space, ED and EU are identical, so we can still use this function.
    return gfcf.compute_ValenciavU_from_ED_and_BU(EU, BU)