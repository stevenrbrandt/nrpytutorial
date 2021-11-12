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
boundL = sp.sympify(0)
boundR = sp.Rational(1,5)

def Ax_FB(x,y,z, **params):
    return sp.sympify(0)

def Ay_FB(x,y,z, **params):
    Ayleft = x - sp.Rational(1,5)
    Aycenter = -sp.sympify(5)*x**2 + x - sp.Rational(1,5)
    Ayright = -x

    out = noif.coord_leq_bound(x,boundL)*Ayleft\
         +noif.coord_greater_bound(x,boundL)*noif.coord_leq_bound(x,boundR)*Aycenter\
         +noif.coord_greater_bound(x,boundR)*Ayright
    return out

def Az_FB(x,y,z, **params):
    return y - Ay_FB(x,y,z, **params)

def zfunc(x):
    # Naming it like this to avoid any possible confusion with the coordinate.
    return -sp.sympify(10)*x + sp.sympify(1)

#Step 3: Compute v^i from B^i and E_i
def ValenciavU_func_FB(**params):
    x = rfm.xx_to_Cart[0]

    BU = ixp.zerorank1(DIM=3)
    BU[0] = sp.sympify(1)
    BU[1] = BU[2] = noif.coord_leq_bound(x,boundL)                                   * sp.sympify(1)\
                   +noif.coord_greater_bound(x,boundL)*noif.coord_leq_bound(x,boundR)* zfunc(x)\
                   +noif.coord_greater_bound(x,boundR)                               * -sp.sympify(1)

    EU = ixp.zerorank1()
    EU[0] = sp.sympify(0)
    EU[1] = sp.Rational(1,2)
    EU[2] = -sp.Rational(1,2)

    # In flat space, ED and EU are identical, so we can still use this function.
    return gfcf.compute_ValenciavU_from_ED_and_BU(EU, BU)