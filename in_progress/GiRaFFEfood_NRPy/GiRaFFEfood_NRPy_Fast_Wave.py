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
bound = sp.Rational(1,10)

def Ax_FW(x,y,z, **params):
    return sp.sympify(0)

def Ay_FW(x,y,z, **params):
    return sp.sympify(0)

def Az_FW(x,y,z, **params):
    # A_z = y+ (-x-0.0075) if x <= -0.1
    #          (0.75x^2 - 0.85x) if -0.1 < x <= 0.1
    #          (-0.7x-0.0075) if x > 0.1
    Azleft = y - x - sp.Rational(75,10000)
    Azcenter = y + sp.Rational(75,100)*x*x - sp.Rational(85,100)*x
    Azright = y - sp.Rational(7,10)*x - sp.Rational(75,10000)

    out = noif.coord_leq_bound(x,-bound)*Azleft\
         +noif.coord_greater_bound(x,-bound)*noif.coord_leq_bound(x,bound)*Azcenter\
         +noif.coord_greater_bound(x,bound)*Azright
    return out

def ValenciavU_func_FW(**params):
    # B^x(0,x) = 1.0
    # B^y(0,x) = 1.0 if x <= -0.1
    #            1.0-1.5(x+0.1) if -0.1 < x <= 0.1
    #            0.7 if x > 0.1
    # B^z(0,x) = 0
    x = rfm.xx_to_Cart[0]
    y = rfm.xx_to_Cart[1]

    Byleft = sp.sympify(1)
    Bycenter = sp.sympify(1) - sp.Rational(15,10)*(x+sp.Rational(1,10))
    Byright = sp.Rational(7,10)

    BU = ixp.zerorank1()
    BU[0] = sp.sympify(1)
    BU[1] = noif.coord_leq_bound(x,-bound)*Byleft\
            +noif.coord_greater_bound(x,-bound)*noif.coord_leq_bound(x,bound)*Bycenter\
            +noif.coord_greater_bound(x,bound)*Byright
    BU[2] = 0

    # E^x(0,x) = 0.0 , E^y(x) = 0.0 , E^z(x) = -B^y(0,x)
    EU = ixp.zerorank1()
    EU[0] = sp.sympify(0)
    EU[1] = sp.sympify(0)
    EU[2] = -BU[1]

    # In flat space, ED and EU are identical, so we can still use this function.
    return gfcf.compute_ValenciavU_from_ED_and_BU(EU, BU)