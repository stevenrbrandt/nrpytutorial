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

mu_AW = par.Cparameters("REAL",thismodule,["mu_AW"], -0.5) # The wave speed
M_PI  = par.Cparameters("#define",thismodule,["M_PI"], "")

gammamu = sp.sympify(1)/sp.sqrt(sp.sympify(1)-mu_AW**2)
bound = sp.Rational(1,10)/gammamu
def g_AW(x):
    return sp.cos(sp.sympify(5)*M_PI*gammamu*x)/M_PI

import Min_Max_and_Piecewise_Expressions as noif

def Ax_AW(x,y,z, **params):
    return sp.sympify(0)

def Ay_AW(x,y,z, **params):
    # \gamma_\mu x - 0.015 if x <= -0.1/\gamma_\mu
    # 1.15 \gamma_\mu x - 0.03g(x) if -0.1/\gamma_\mu < x <= 0.1/\gamma_\mu
    # 1.3 \gamma_\mu x - 0.015 if x > 0.1/\gamma_\mu
    Ayleft = gammamu*x - sp.Rational(15,1000)
    Aycenter = sp.Rational(115,100)*gammamu*x - sp.Rational(3,100)*g_AW(x)
    Ayright = sp.Rational(13,10)*gammamu*x - sp.Rational(15,1000)

    out = noif.coord_leq_bound(x,-bound)*Ayleft\
         +noif.coord_greater_bound(x,-bound)*noif.coord_leq_bound(x,bound)*Aycenter\
         +noif.coord_greater_bound(x,bound)*Ayright
    return out

def Az_AW(x,y,z, **params):
    # y - \gamma_\mu (1-\mu)x
    return y-gammamu*(sp.sympify(1)-mu_AW)*x

def f_AW(x):
    xprime = gammamu*x
    return 1 + sp.sin(5*M_PI*xprime)

#Step 3: Compute v^i from B^i and E_i
def ValenciavU_func_AW(**params):
    x = rfm.xx_to_Cart[0]

    Bzleft = sp.sympify(1)
    Bzcenter = sp.sympify(1) + sp.Rational(15,100)*f_AW(x)
    Bzright = sp.Rational(13,10)

    BpU = ixp.zerorank1()
    BpU[0] = sp.sympify(1)
    BpU[1] = sp.sympify(1)
    BpU[2] = noif.coord_leq_bound(x,-bound)*Bzleft\
            +noif.coord_greater_bound(x,-bound)*noif.coord_leq_bound(x,bound)*Bzcenter\
            +noif.coord_greater_bound(x,bound)*Bzright

    EpU = ixp.zerorank1()
    EpU[0] = -BpU[2]
    EpU[1] = sp.sympify(0)
    EpU[2] = sp.sympify(1)

    BU = ixp.zerorank1()
    BU[0] = BpU[0]
    BU[1] = gammamu*(BpU[1]-mu_AW*EpU[2])
    BU[2] = gammamu*(BpU[2]+mu_AW*EpU[1])

    EU = ixp.zerorank1()
    EU[0] = EpU[0]
    EU[1] = gammamu*(EpU[1]+mu_AW*BpU[2])
    EU[2] = gammamu*(EpU[2]-mu_AW*BpU[1])

    # In flat space, ED and EU are identical, so we can still use this function.
    return gfcf.compute_ValenciavU_from_ED_and_BU(EU, BU)