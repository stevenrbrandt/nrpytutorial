# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)
giraffefood_dir_path = os.path.join("in_progress","GiRaFFEfood_NRPy")
if giraffefood_dir_path not in sys.path:
    sys.path.append(giraffefood_dir_path)
giraffefood_dir_path = os.path.join("GiRaFFEfood_NRPy")
if giraffefood_dir_path not in sys.path:
    sys.path.append(giraffefood_dir_path)

# Step 0: Import the NRPy+ core modules and set the reference metric to Cartesian
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import reference_metric as rfm   # NRPy+: Reference metric support

# Use the Jacobian matrix to transform the vectors to Cartesian coordinates.
# Construct Jacobian & Inverse Jacobians:
par.set_parval_from_str("reference_metric::CoordSystem","Spherical")
rfm.reference_metric()
Jac_dUCart_dDrfmUD,Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()
# Transform the coordinates of the Jacobian matrix from spherical to Cartesian:
par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()
tmpa,tmpb,tmpc = sp.symbols("tmpa,tmpb,tmpc")
for i in range(3):
    for j in range(3):
        Jac_dUCart_dDrfmUD[i][j] = Jac_dUCart_dDrfmUD[i][j].subs([(rfm.xx[0],tmpa),(rfm.xx[1],tmpb),(rfm.xx[2],tmpc)])
        Jac_dUCart_dDrfmUD[i][j] = Jac_dUCart_dDrfmUD[i][j].subs([(tmpa,rfm.xxSph[0]),(tmpb,rfm.xxSph[1]),(tmpc,rfm.xxSph[2])])
        Jac_dUrfm_dDCartUD[i][j] = Jac_dUrfm_dDCartUD[i][j].subs([(rfm.xx[0],tmpa),(rfm.xx[1],tmpb),(rfm.xx[2],tmpc)])
        Jac_dUrfm_dDCartUD[i][j] = Jac_dUrfm_dDCartUD[i][j].subs([(tmpa,rfm.xxSph[0]),(tmpb,rfm.xxSph[1]),(tmpc,rfm.xxSph[2])])


# Step 1a: Set commonly used parameters.
thismodule = __name__

# Generic function for all 1D tests: Compute Ax,Ay,Az
def Axyz_func_Cartesian(Ax_func,Ay_func,Az_func, stagger_enable, **params):
    x = rfm.xx_to_Cart[0]
    y = rfm.xx_to_Cart[1]
    z = rfm.xx_to_Cart[2]
    AD = ixp.zerorank1()
    # First Ax
    if stagger_enable:
        y += sp.Rational(1,2)*gri.dxx[1]
        z += sp.Rational(1,2)*gri.dxx[2]
    AD[0] = Ax_func(x,y,z, **params)
    # Then Ay
    if stagger_enable:
        x += sp.Rational(1,2)*gri.dxx[0]
        y -= sp.Rational(1,2)*gri.dxx[1]
    AD[1] = Ay_func(x,y,z, **params)
    # Finally Az
    if stagger_enable:
        y += sp.Rational(1,2)*gri.dxx[1]
        z -= sp.Rational(1,2)*gri.dxx[2]
    AD[2] = Az_func(x,y,z, **params)

    return AD

# Generic function for all 1D tests: Compute Ax,Ay,Az
def Axyz_func_spherical(Ar_func,At_func,Ap_func, stagger_enable, **params):
    if "KerrSchild_radial_shift" in params:
        KerrSchild_radial_shift = params["KerrSchild_radial_shift"]
        r = rfm.xxSph[0] + KerrSchild_radial_shift # We are setting the data up in Shifted Kerr-Schild coordinates
    else:
        r = rfm.xxSph[0] # Some other coordinate system
    theta = rfm.xxSph[1]
    phi   = rfm.xxSph[2]
    AsphD = ixp.zerorank1()
    # First Ax
    AsphD[0] = Ar_func(r,theta,phi, **params)
    # Then Ay
    AsphD[1] = At_func(r,theta,phi, **params)
    # Finally Az
    AsphD[2] = Ap_func(r,theta,phi, **params)

    # Use the Jacobian matrix to transform the vectors to Cartesian coordinates.
    AD = rfm.basis_transform_vectorD_from_rfmbasis_to_Cartesian(Jac_dUrfm_dDCartUD, AsphD)
#     M = params["M"]
#     AD[2] = fp_of_r(rfm.xxSph[0] + KerrSchild_radial_shift,M)
    if stagger_enable:
        # First Ax
        AD[0] = AD[0].subs(rfm.xx[1],rfm.xx[1]+sp.Rational(1,2)*gri.dxx[1]).subs(rfm.xx[2],rfm.xx[2]+sp.Rational(1,2)*gri.dxx[2])
        # Then Ay
        AD[1] = AD[1].subs(rfm.xx[0],rfm.xx[0]+sp.Rational(1,2)*gri.dxx[0]).subs(rfm.xx[2],rfm.xx[2]+sp.Rational(1,2)*gri.dxx[2])
        # Finally Az
        AD[2] = AD[2].subs(rfm.xx[0],rfm.xx[0]+sp.Rational(1,2)*gri.dxx[0]).subs(rfm.xx[1],rfm.xx[1]+sp.Rational(1,2)*gri.dxx[1])
    return AD

# Generic function for all 1D tests: Valencia 3-velocity from ED and BU
def compute_ValenciavU_from_ED_and_BU(ED, BU, gammaDD=None):
    # Now, we calculate v^i = ([ijk] E_j B_k) / B^2,
    # where [ijk] is the Levi-Civita symbol and B^2 = \gamma_{ij} B^i B^j$ is a trivial dot product in flat space.

    # In flat spacetime, use the Minkowski metric; otherwise, use the input metric.
    if gammaDD is None:
        gammaDD = ixp.zerorank2()
        for i in range(3):
            gammaDD[i][i] = sp.sympify(1)

    unused_gammaUU,gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    sqrtgammaDET = sp.sqrt(gammaDET)
    LeviCivitaTensorUUU = ixp.LeviCivitaTensorUUU_dim3_rank3(sqrtgammaDET)

    BD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            BD[i] += gammaDD[i][j]*BU[j]
    B2 = sp.sympify(0)
    for i in range(3):
        B2 += BU[i] * BD[i]

    ValenciavU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                ValenciavU[i] += LeviCivitaTensorUUU[i][j][k] * ED[j] * BD[k] / B2

    return ValenciavU
