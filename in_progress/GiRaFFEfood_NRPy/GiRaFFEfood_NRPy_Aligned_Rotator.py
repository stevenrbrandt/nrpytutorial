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

# The angular velocity of the "neutron star"
Omega_aligned_rotator = par.Cparameters("REAL",thismodule,"Omega_aligned_rotator",1e3)
B_p_aligned_rotator,R_NS_aligned_rotator = par.Cparameters("REAL",thismodule,
                                                           # B_p_aligned_rotator = the intensity of the magnetic field and
                                                           # R_NS_aligned_rotator= "Neutron star" radius
                                                           ["B_p_aligned_rotator","R_NS_aligned_rotator"],
                                                           [1e-5, 1.0])

# Step 2: Set the vectors A and E in Spherical coordinates
def Ar_AR(r,theta,phi, **params):
    return sp.sympify(0)

def Ath_AR(r,theta,phi, **params):
    return sp.sympify(0)

def Aph_AR(r,theta,phi, **params):
    # \mu \varpi / r^3
    # \varpi = sqrt(x^2+y^2) = r \sin(\theta)
    varpi = r * sp.sin(theta)
    mu = B_p_aligned_rotator * R_NS_aligned_rotator**3 / 2
    return mu * varpi**2 / (r**3)

import Min_Max_and_Piecewise_Expressions as noif
#Step 3: Compute v^i
def ValenciavU_func_AR(**params):
    LeviCivitaSymbolDDD = ixp.LeviCivitaSymbol_dim3_rank3()

    unit_zU = ixp.zerorank1()
    unit_zU[2] = sp.sympify(1)

    r = rfm.xxSph[0]

    ValenciavU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                ValenciavU[i] += noif.coord_leq_bound(r,R_NS_aligned_rotator)*LeviCivitaSymbolDDD[i][j][k] * Omega_aligned_rotator * unit_zU[j] * rfm.xx[k]

    return ValenciavU