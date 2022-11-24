# This module performs the conversion from ADM
# spacetime variables in Spherical or Cartesian coordinates,
# to BSSN quantities in any basis supported by NRPy+.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Initialize core Python/NRPy+ modules
from outputC import outputC, lhrh, add_to_Cfunction_dict # NRPy+: Core C code output module
from outputC import outC_NRPy_basic_defines_h_dict
import NRPy_param_funcs as par     # NRPy+: Parameter interface
import finite_difference as fin    # NRPy+: Finite difference C code generation module
import grid as gri                 # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp           # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm     # NRPy+: Reference metric support
import BSSN.BSSN_quantities as Bq  # NRPy+: Computes useful BSSN quantities; e.g., gammabarUU & GammabarUDD needed below
from pickling import pickle_NRPy_env  # NRPy+: Pickle/unpickle NRPy+ environment, for parallel codegen
import sys                         # Standard Python module for multiplatform OS-level functions


def add_to_Cfunction_dict_initial_data_reader__convert_to_BSSN_from_ADM_sph_or_Cart(input_Coord="Spherical"):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    c_type = "void"

    output_Coord = par.parval_from_str("reference_metric::CoordSystem")
    desc = "Read in ADM initial data in the " + input_Coord + " basis, and convert to BSSN data in " + output_Coord + " coordinates"
    name = "initial_data_reader__convert_to_BSSN_from_ADM_" + input_Coord
    params = """griddata_struct *restrict griddata, ID_persist_struct *restrict ID_persist,
                                                             void ID_function(const paramstruct *params, const REAL xCart[3],
                                                                              const ID_persist_struct *restrict ID_persist,
                                                                              ID_output_struct *restrict ID_output)"""

    body = r"""
  const int Nxx_plus_2NGHOSTS0 = griddata.params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata.params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata.params.Nxx_plus_2NGHOSTS2;

  LOOP_OMP("omp parallel for", i0,0,Nxx_plus_2NGHOSTS0, i1,0,Nxx_plus_2NGHOSTS1, i2,0,Nxx_plus_2NGHOSTS2) {
    // xCart is the global Cartesian coordinate, which accounts for any grid offsets from the origin.
    REAL xCart[3];  xx_to_Cart(&griddata.params, griddata.xx, i0,i1,i2, xCart);

    // Read or compute initial data at destination point xCart
    ID_output_struct ID_output;
    ID_function(&griddata.params, xCart, ID_persist, &ID_output);

    // Unpack ID_output for scalar alpha
    const REAL alpha = ID_output.alpha;

    // Unpack ID_output for ADM vectors/tensors
"""
    for i in ["betaSphorCartU", "BSphorCartU"]:
        for j in range(3):
            varname = i + str(j)
            body += "    const REAL " + varname + " = ID_output." + varname + ";\n"
        body += "\n"
    for i in ["gammaSphorCartDD", "KSphorCartDD"]:
        for j in range(3):
            for k in range(j, 3):
                varname = i + str(j) + str(k)
                body += "    const REAL " + varname + " = ID_output." + varname + ";\n"
        body += "\n"

    body += r"""
    // Perform the basis transform on ADM vectors/tensors from """+input_Coord+""" to Cartesian:
    //   (Don't be surprised if it's trivial when ADM quantities are already in the Cartesian basis.)
    REAL betaCartU0,betaCartU1,betaCartU2;
    REAL BCartU0,BCartU1,BCartU2;
    REAL gammaCartDD00,gammaCartDD01,gammaCartDD02,gammaCartDD11,gammaCartDD12,gammaCartDD22;
    REAL KCartDD00,KCartDD01,KCartDD02,KCartDD11,KCartDD12,KCartDD22;
    {
      // Set destination point
      const REAL Cartx = xCart[0];
      const REAL Carty = xCart[1];
      const REAL Cartz = xCart[2];

      // Set destination xx[3]
"""
    # Set reference_metric to the input_Coord
    par.set_parval_from_str("reference_metric::CoordSystem", input_Coord)
    rfm.reference_metric()

    body += outputC(rfm.Cart_to_xx[:3], ["const REAL xx0", "const REAL xx1", "const REAL xx2"],
                    filename="returnstring",
                    params="outCverbose=False,includebraces=False,preindent=3,CSE_varprefix=tmp_xx")

    # Define the input variables:
    gammaSphorCartDD = ixp.declarerank2("gammaSphorCartDD", "sym01")
    KSphorCartDD     = ixp.declarerank2("KSphorCartDD", "sym01")
    betaSphorCartU = ixp.declarerank1("betaSphorCartU")
    BSphorCartU    = ixp.declarerank1("BSphorCartU")

    # Compute Jacobian to convert to Cartesian coordinates
    Jac_dUCart_dDrfmUD, Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

    gammaCartDD = rfm.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(Jac_dUrfm_dDCartUD, gammaSphorCartDD)
    KCartDD = rfm.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(Jac_dUrfm_dDCartUD, KSphorCartDD)
    betaCartU = rfm.basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, betaSphorCartU)
    BCartU = rfm.basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, BSphorCartU)

    list_of_output_exprs    = []
    list_of_output_varnames = []
    for i in range(3):
        list_of_output_exprs += [betaCartU[i]]
        list_of_output_varnames += ["betaCartU" + str(i)]
        list_of_output_exprs += [BCartU[i]]
        list_of_output_varnames += ["BCartU" + str(i)]
        for j in range(i, 3):
            list_of_output_exprs += [gammaCartDD[i][j]]
            list_of_output_varnames += ["gammaCartDD" + str(i) + str(j)]
            list_of_output_exprs += [KCartDD[i][j]]
            list_of_output_varnames += ["KCartDD" + str(i) + str(j)]

    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = \
        (list(t) for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs))))

    body += outputC(list_of_output_exprs, list_of_output_varnames,
                    filename="returnstring", params="outCverbose=False,includebraces=False,preindent=4")
    body += r"""
    }
"""
    # Set reference_metric to Cartesian
    par.set_parval_from_str("reference_metric::CoordSystem", "Cartesian")
    rfm.reference_metric()

    # Reset these variables, as they have been defined above in the C code.
    gammaCartDD = ixp.declarerank2("gammaCartDD", "sym01")
    KCartDD     = ixp.declarerank2("KCartDD", "sym01")

    import BSSN.BSSN_in_terms_of_ADM as BitoA
    BitoA.trK_AbarDD_aDD(gammaCartDD, KCartDD)
    BitoA.gammabarDD_hDD(gammaCartDD)
    BitoA.cf_from_gammaDD(gammaCartDD)

    body += r"""
      // Next convert ADM quantities gammaDD & KDD
      //   into BSSN gammabarDD, AbarDD, cf, and trK, in the Cartesian basis.
      REAL gammabarCartDD00,gammabarCartDD01,gammabarCartDD02,gammabarCartDD11,gammabarCartDD12,gammabarCartDD22;
      REAL AbarCartDD00,AbarCartDD01,AbarCartDD02,AbarCartDD11,AbarCartDD12,AbarCartDD22;
      REAL cf, trK;
      {
"""
    list_of_output_exprs    = [BitoA.cf, BitoA.trK]
    list_of_output_varnames = ["cf", "trK"]
    for i in range(3):
        for j in range(i, 3):
            list_of_output_exprs += [BitoA.gammabarDD[i][j]]
            list_of_output_varnames += ["gammabarCartDD" + str(i) + str(j)]
            list_of_output_exprs += [BitoA.AbarDD[i][j]]
            list_of_output_varnames += ["AbarCartDD" + str(i) + str(j)]
    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = \
        (list(t) for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs))))
    body += outputC(list_of_output_exprs, list_of_output_varnames,
                    filename="returnstring", params="outCverbose=False,includebraces=False,preindent=4")
    body += r"""
      }

      const int idx3 = IDX3S(i0,i1,i2);

      // First set the BSSN scalars, as these don't need a basis transform:
      griddata.gridfuncs.y_n_gfs[IDX4ptS(ALPHAGF, idx3)] = alpha;
      griddata.gridfuncs.y_n_gfs[IDX4ptS(TRKGF, idx3)] = trK;
      griddata.gridfuncs.y_n_gfs[IDX4ptS(CFGF, idx3)] = cf;

      // Then set the BSSN vectors/tensors, which require we perform basis transform & rescaling:
      initial_data_BSSN_basis_transform_Cartesian_to_rfm_and_rescale
        (&griddata.params, griddata.xx[0][i0],griddata.xx[1][i1],griddata.xx[2][i2],
         betaCartU0,betaCartU1,betaCartU2, BCartU0,BCartU1,BCartU2,
         gammabarCartDD00,gammabarCartDD01,gammabarCartDD02,gammabarCartDD11,gammabarCartDD12,gammabarCartDD22,
         AbarCartDD00,AbarCartDD01,AbarCartDD02,AbarCartDD11,AbarCartDD12,AbarCartDD22,
         idx3, griddata.gridfuncs.y_n_gfs);
    } // END LOOP over all gridpoints on given grid

    initial_data_lambdaU_grid_interior(&griddata.params, griddata.xx,
                                               griddata.gridfuncs.y_n_gfs);

  }  // END LOOP over all grids
"""

    # Restore reference_metric to output_Coord
    par.set_parval_from_str("reference_metric::CoordSystem", output_Coord)
    rfm.reference_metric()

    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=False)
    return pickle_NRPy_env()


# By the time this function is called, all BSSN tensors and vectors are in the Cartesian
# coordinate basis $x^i_{\rm Cart} = (x,y,z)$, but we need them in the curvilinear
# coordinate basis $x^i_{\rm rfm}=$`(xx0,xx1,xx2)` set by the
#  `"reference_metric::CoordSystem"` variable.
def add_to_Cfunction_dict_initial_data_BSSN_basis_transform_Cartesian_to_rfm_and_rescale():
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]  # Need NRPy_function_prototypes.h for ChebSpherical
    c_type = "void"

    output_Coord = par.parval_from_str("reference_metric::CoordSystem")
    desc = "Basis transform AbarDD and gammabarDD from Cartesian to " + output_Coord + " coordinates"
    name = "initial_data_BSSN_basis_transform_Cartesian_to_rfm_and_rescale"
    params = """const paramstruct *restrict params, const REAL xx0,const REAL xx1,const REAL xx2,
         const REAL betaCartU0,const REAL betaCartU1,const REAL betaCartU2,
         const REAL BCartU0,const REAL BCartU1,const REAL BCartU2,
         const REAL gammabarCartDD00,const REAL gammabarCartDD01,const REAL gammabarCartDD02,
         const REAL gammabarCartDD11,const REAL gammabarCartDD12,const REAL gammabarCartDD22,
         const REAL AbarCartDD00,const REAL AbarCartDD01,const REAL AbarCartDD02,
         const REAL AbarCartDD11,const REAL AbarCartDD12,const REAL AbarCartDD22,
         const int idx3, REAL *restrict y_n_gfs"""

    # Define the input variables:
    gammabarCartDD = ixp.declarerank2("gammabarCartDD", "sym01")
    AbarCartDD     = ixp.declarerank2("AbarCartDD", "sym01")
    betaCartU = ixp.declarerank1("betaCartU")
    BCartU    = ixp.declarerank1("BCartU")

    # Set reference_metric to the output_Coord
    par.set_parval_from_str("reference_metric::CoordSystem", output_Coord)
    rfm.reference_metric()
    # Compute Jacobian to convert to Cartesian coordinates
    Jac_dUCart_dDrfmUD, Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

    gammabarDD = rfm.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(Jac_dUCart_dDrfmUD, gammabarCartDD)
    AbarDD = rfm.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(Jac_dUCart_dDrfmUD, AbarCartDD)
    betaU = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, betaCartU)
    BU = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, BCartU)

    # Next rescale:
    vetU = ixp.zerorank1()
    betU = ixp.zerorank1()
    hDD  = ixp.zerorank2()
    aDD  = ixp.zerorank2()
    for i in range(3):
        vetU[i] = betaU[i] / rfm.ReU[i]
        betU[i] =    BU[i] / rfm.ReU[i]
        for j in range(3):
            hDD[i][j] = (gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[i][j]
            aDD[i][j] = AbarDD[i][j] / rfm.ReDD[i][j]

    def gfaccess(gfname):
        return "y_n_gfs[IDX4ptS("+gfname.upper()+"GF, idx3)]"
    list_of_output_exprs    = []
    list_of_output_varnames = []
    for i in range(3):
        list_of_output_exprs += [vetU[i]]
        list_of_output_varnames += [gfaccess("vetU" + str(i))]
        list_of_output_exprs += [betU[i]]
        list_of_output_varnames += [gfaccess("betU" + str(i))]
        for j in range(i, 3):
            list_of_output_exprs += [hDD[i][j]]
            list_of_output_varnames += [gfaccess("hDD" + str(i) + str(j))]
            list_of_output_exprs += [aDD[i][j]]
            list_of_output_varnames += [gfaccess("aDD" + str(i) + str(j))]
    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = \
        (list(t) for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs))))

    body = outputC(list_of_output_exprs, list_of_output_varnames,
                    filename="returnstring", params="outCverbose=False,includebraces=False,preindent=1")

    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=True)
    return pickle_NRPy_env()


# initial_data_lambdaU_grid_interior() computes lambdaU from
#  finite-difference derivatives of rescaled metric quantities
def add_to_Cfunction_dict_initial_data_lambdaU_grid_interior():
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]  # Need NRPy_function_prototypes.h for ChebSpherical
    c_type = "void"

    output_Coord = par.parval_from_str("reference_metric::CoordSystem")
    desc = "Compute lambdaU in " + output_Coord + " coordinates"
    name = "initial_data_lambdaU_grid_interior"
    params = """const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict in_gfs"""
    # Step 7: Compute $\bar{\Lambda}^i$ from finite-difference derivatives of rescaled metric quantities

    # We will need all BSSN gridfunctions to be defined, as well as
    #     expressions for gammabarDD_dD in terms of exact derivatives of
    #     the rescaling matrix and finite-difference derivatives of
    #     hDD's. This functionality is provided by BSSN.BSSN_unrescaled_and_barred_vars,
    #     which we call here to overwrite above definitions of gammabarDD,gammabarUU, etc.
    Bq.gammabar__inverse_and_derivs()  # Provides gammabarUU and GammabarUDD
    gammabarUU    = Bq.gammabarUU
    GammabarUDD   = Bq.GammabarUDD

    # Next evaluate \bar{\Lambda}^i, based on GammabarUDD above and GammahatUDD
    #       (from the reference metric):
    LambdabarU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                LambdabarU[i] += gammabarUU[j][k] * (GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k])

    # Finally apply rescaling:
    # lambda^i = Lambdabar^i/\text{ReU[i]}
    lambdaU = ixp.zerorank1()
    for i in range(3):
        lambdaU[i] = LambdabarU[i] / rfm.ReU[i]

    lambdaU_expressions = [lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU0"), rhs=lambdaU[0]),
                           lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU1"), rhs=lambdaU[1]),
                           lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU2"), rhs=lambdaU[2])]
    body = fin.FD_outputC("returnstring", lambdaU_expressions,
                           params="outCverbose=False,includebraces=False,preindent=0")
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        loopopts="InteriorPoints,Read_xxs",
        enableCparameters=True)
    return pickle_NRPy_env()


def add_to_Cfunction_dict_exact_ADM_ID_function(IDtype, IDCoordSystem, alpha, betaU, BU, gammaDD, KDD):
    includes = ["NRPy_basic_defines.h"]
    desc = IDtype + " initial data"
    c_type = "void"
    name = IDtype
    params = "const paramstruct *params, const REAL xCart[3], const ID_persist_struct *restrict ID_persist, ID_output_struct *restrict ID_output"
    orig_Coord = par.parval_from_str("reference_metric::CoordSystem")
    par.set_parval_from_str("reference_metric::CoordSystem", IDCoordSystem)
    rfm.reference_metric()
    body = ""
    if IDCoordSystem == "Spherical":
        body += r"""  const REAL Cartx=xCart[0], Carty=xCart[1], Cartz=xCart[2];
  REAL xx0,xx1,xx2 __attribute__((unused));  // xx2 might be unused in the case of axisymmetric initial data.
  {
""" + outputC(rfm.Cart_to_xx[:3], ["xx0", "xx1", "xx2"], filename="returnstring",
              params="outCverbose=False,includebraces=False,preindent=2") + """
  }
"""
    elif IDCoordSystem == "Cartesian":
        body += r"""  const REAL xx0=xCart[0], xx1=xCart[1], xx2=xCart[2];
"""
    else:
        print("add_to_Cfunction_dict_exact_ADM_ID_function() Error: IDCoordSystem == " + IDCoordSystem + " unsupported")
        sys.exit(1)
    list_of_output_exprs = [alpha]
    list_of_output_varnames = ["ID_output->alpha"]
    for i in range(3):
        list_of_output_exprs += [betaU[i]]
        list_of_output_varnames += ["ID_output->betaSphorCartU" + str(i)]
        list_of_output_exprs += [BU[i]]
        list_of_output_varnames += ["ID_output->BSphorCartU" + str(i)]
        for j in range(i, 3):
            list_of_output_exprs += [gammaDD[i][j]]
            list_of_output_varnames += ["ID_output->gammaSphorCartDD" + str(i) + str(j)]
            list_of_output_exprs += [KDD[i][j]]
            list_of_output_varnames += ["ID_output->KSphorCartDD" + str(i) + str(j)]
    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = \
        (list(t) for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs))))

    body += outputC(list_of_output_exprs, list_of_output_varnames,
                    filename="returnstring", params="outCverbose=False,includebraces=False,preindent=1")

    # Restore CoordSystem:
    par.set_parval_from_str("reference_metric::CoordSystem", orig_Coord)
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc, c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=True)
    return pickle_NRPy_env()


# Other than its core use as a means to store ADM input quantities,
# `ID_output_struct` is designed to be extensible. For example, it may be
# used to store e.g., pseudospectral coefficients for TwoPunctures,
# initial data gridfunctions from NRPyElliptic, pointers to TOV 1D data
# from the TOV solver, etc.
def register_NRPy_basic_defines(ID_persist_struct_contents_str=""):
    Nbd = r"""typedef struct __ID_output_struct__ {
  REAL alpha;

  REAL betaSphorCartU0, betaSphorCartU1, betaSphorCartU2;
  REAL BSphorCartU0, BSphorCartU1, BSphorCartU2;

  REAL gammaSphorCartDD00, gammaSphorCartDD01, gammaSphorCartDD02;
  REAL gammaSphorCartDD11, gammaSphorCartDD12, gammaSphorCartDD22;

  REAL KSphorCartDD00, KSphorCartDD01, KSphorCartDD02;
  REAL KSphorCartDD11, KSphorCartDD12, KSphorCartDD22;
} ID_output_struct;
"""
    if ID_persist_struct_contents_str == "":
        Nbd += "typedef struct __ID_persist_struct__ {\n"
        Nbd += "} ID_persist_struct;\n"
    else:
        Nbd += ID_persist_struct_contents_str + "\n"
    outC_NRPy_basic_defines_h_dict["BSSN_initial_data"] = Nbd
