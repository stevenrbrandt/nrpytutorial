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

mu_DAW = par.Cparameters("REAL",thismodule,["mu_DAW"], -0.5) # The wave speed
M_PI  = par.Cparameters("#define",thismodule,["M_PI"], "")

gammamu = sp.sympify(1)/sp.sqrt(sp.sympify(1)-mu_DAW**2)
bound = sp.Rational(1,10)/gammamu

def h1_DAW(x):
    return sp.cos(sp.Rational(5,2)*M_PI*(gammamu*x+sp.Rational(1,10)))
def h2_DAW(x):
    return sp.sin(sp.Rational(5,2)*M_PI*(gammamu*x+sp.Rational(1,10)))

import Min_Max_and_Piecewise_Expressions as noif

def Ax_DAW(x,y,z, **params):
    return sp.sympify(0)

def Ay_DAW(x,y,z, **params):
    Ayleft = -sp.Rational(4,5)/M_PI
    Aycenter = -sp.Rational(4,5)/M_PI * h1_DAW(x)
    Ayright = sp.sympify(2)*(gammamu*x-sp.Rational(1,10))

    out = noif.coord_leq_bound(x,-bound)*Ayleft\
         +noif.coord_greater_bound(x,-bound)*noif.coord_leq_bound(x,bound)*Aycenter\
         +noif.coord_greater_bound(x,bound)*Ayright
    return out

def Az_DAW(x,y,z, **params):
    Azleft = -sp.sympify(2)*(gammamu*x+sp.Rational(1,10))
    Azcenter = -sp.Rational(4,5)/M_PI * h2_DAW(x)
    Azright = -sp.Rational(4,5)/M_PI

    out = noif.coord_leq_bound(x,-bound)*Azleft\
         +noif.coord_greater_bound(x,-bound)*noif.coord_leq_bound(x,bound)*Azcenter\
         +noif.coord_greater_bound(x,bound)*Azright
    return out

def phi(x):
    xprime = gammamu*x
    bound = sp.Rational(1,10)

    phileft = sp.sympify(0)
    phicenter = sp.Rational(5,2)*M_PI*(xprime+sp.Rational(1,10))
    phiright = sp.Rational(1,2)*M_PI

    out = noif.coord_leq_bound(xprime,-bound)*phileft\
         +noif.coord_greater_bound(xprime,-bound)*noif.coord_leq_bound(x,bound)*phicenter\
         +noif.coord_greater_bound(xprime,bound)*phiright
    return out

#Step 3: Compute v^i from B^i and E_i
def ValenciavU_func_DAW(**params):
    x = rfm.xx_to_Cart[0]

    BpU = ixp.zerorank1()
    BpU[0] = sp.sympify(0)
    BpU[1] = sp.sympify(2)*sp.cos(phi(x))
    BpU[2] = sp.sympify(2)*sp.sin(phi(x))

    EpU = ixp.zerorank1()

    BU = ixp.zerorank1()
    BU[0] = BpU[0]
    BU[1] = gammamu*(BpU[1]-mu_DAW*EpU[2])
    BU[2] = gammamu*(BpU[2]+mu_DAW*EpU[1])

    EU = ixp.zerorank1()
    EU[0] = EpU[0]
    EU[1] = gammamu*(EpU[1]+mu_DAW*BpU[2])
    EU[2] = gammamu*(EpU[2]-mu_DAW*BpU[1])

    # In flat space, ED and EU are identical, so we can still use this function.
    return gfcf.compute_ValenciavU_from_ED_and_BU(EU, BU)