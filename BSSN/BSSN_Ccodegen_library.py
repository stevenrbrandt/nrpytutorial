# Generating full C library for solving
#  the BSSN equations in
#  ***curvilinear*** coordinates, using a
#  reference-metric formalism
#
# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
#
# License: BSD 2-Clause

# COMPLETE DOCUMENTATION (JUPYTER NOTEBOOKS):
# START PAGE (start here!):  ../NRPy+_Tutorial.ipynb
# THIS MODULE: ../Tutorial-BSSN_time_evolution-C_codegen_library.ipynb

# RULES FOR ADDING FUNCTIONS TO THIS ROUTINE:
# 1. The function must be runnable from a multiprocessing environment,
#    which means that the function
# 1.a: cannot depend on previous function calls.
# 1.b: cannot create directories (this is not multiproc friendly)


# Step P1: Import needed NRPy+ core modules:
from outputC import lhrh, add_to_Cfunction_dict  # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
from pickling import pickle_NRPy_env   # NRPy+: Pickle/unpickle NRPy+ environment, for parallel codegen
import os, time             # Standard Python modules for multiplatform OS-level functions, benchmarking
import BSSN.BSSN_RHSs as rhs
import BSSN.BSSN_gauge_RHSs as gaugerhs
import BSSN.Enforce_Detgammahat_Constraint as EGC
import loop as lp

###############################################
# print_msg_with_timing() gives the user an idea of what's going on/taking so long. Also outputs timing info.
def print_msg_with_timing(desc, msg="Symbolic", startstop="start", starttime=0.0):
    CoordSystem = par.parval_from_str("reference_metric::CoordSystem")
    elapsed = time.time()-starttime
    if msg == "Symbolic":
        if startstop == "start":
            print("Generating symbolic expressions for " + desc + " (%s coords)..." % CoordSystem)
            return time.time()
        else:
            print("Finished generating symbolic expressions for "+desc+
                  " (%s coords) in %.1f seconds. Next up: C codegen..." % (CoordSystem, elapsed))
    elif msg == "Ccodegen":
        if startstop == "start":
            print("Generating C code for "+desc+" (%s coords)..." % CoordSystem)
            return time.time()
        else:
            print("Finished generating C code for "+desc+" (%s coords) in %.1f seconds." % (CoordSystem, elapsed))


# get_loopopts() sets up options for NRPy+'s loop module
def get_loopopts(points_to_update, enable_SIMD, enable_rfm_precompute, gridsuffix):
    loopopts = points_to_update
    if enable_SIMD:
        loopopts += ",EnableSIMD"
    if enable_rfm_precompute:
        loopopts += ",Enable_rfm_precompute"
    else:
        loopopts += ",Read_xxs"
    if gridsuffix != "":
        loopopts += ","+gridsuffix.replace("_", "")
    return loopopts


# register_stress_energy_source_terms_return_T4UU() registers gridfunctions
#        for T4UU if needed and not yet registered.
def register_stress_energy_source_terms_return_T4UU(enable_stress_energy_source_terms):
    if enable_stress_energy_source_terms:
        registered_already = False
        for i in range(len(gri.glb_gridfcs_list)):
            if gri.glb_gridfcs_list[i].name == "T4UU00":
                registered_already = True
        if not registered_already:
            return ixp.register_gridfunctions_for_single_rank2("AUXEVOL", "T4UU", "sym01", DIM=4)
        else:
            return ixp.declarerank2("T4UU", "sym01", DIM=4)
    return None

#################################
#################################

def BSSN_RHSs__generate_symbolic_expressions(LapseCondition="OnePlusLog",
                                             ShiftCondition="GammaDriving2ndOrder_Covariant",
                                             enable_KreissOliger_dissipation=True,
                                             enable_stress_energy_source_terms=False,
                                             leave_Ricci_symbolic=True):
    ######################################
    # START: GENERATE SYMBOLIC EXPRESSIONS
    starttime = print_msg_with_timing("BSSN_RHSs", msg="Symbolic", startstop="start")

    # Returns None if enable_stress_energy_source_terms==False; otherwise returns symb expressions for T4UU
    T4UU = register_stress_energy_source_terms_return_T4UU(enable_stress_energy_source_terms)

    # Evaluate BSSN RHSs:
    import BSSN.BSSN_quantities as Bq
    par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic", str(leave_Ricci_symbolic))
    rhs.BSSN_RHSs()

    if enable_stress_energy_source_terms:
        import BSSN.BSSN_stress_energy_source_terms as Bsest
        Bsest.BSSN_source_terms_for_BSSN_RHSs(T4UU)
        rhs.trK_rhs += Bsest.sourceterm_trK_rhs
        for i in range(3):
            # Needed for Gamma-driving shift RHSs:
            rhs.Lambdabar_rhsU[i] += Bsest.sourceterm_Lambdabar_rhsU[i]
            # Needed for BSSN RHSs:
            rhs.lambda_rhsU[i] += Bsest.sourceterm_lambda_rhsU[i]
            for j in range(3):
                rhs.a_rhsDD[i][j] += Bsest.sourceterm_a_rhsDD[i][j]

    par.set_parval_from_str("BSSN.BSSN_gauge_RHSs::LapseEvolutionOption", LapseCondition)
    par.set_parval_from_str("BSSN.BSSN_gauge_RHSs::ShiftEvolutionOption", ShiftCondition)
    gaugerhs.BSSN_gauge_RHSs()  # Can depend on above RHSs
    # Restore BSSN.BSSN_quantities::LeaveRicciSymbolic to False
    par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic", "False")

    # Add Kreiss-Oliger dissipation to the BSSN RHSs:
    if enable_KreissOliger_dissipation:
        thismodule = "KO_Dissipation"
        # diss_strength = Bq.cf  # <- another attempt. Maybe should multiply by diss_strength?
        diss_strength = par.Cparameters("REAL", thismodule, "diss_strength", 0.2)*Bq.cf # *Bq.cf*Bq.cf*Bq.cf # cf**1 is found better than cf**4 over the long term.

        alpha_dKOD = ixp.declarerank1("alpha_dKOD")
        cf_dKOD = ixp.declarerank1("cf_dKOD")
        trK_dKOD = ixp.declarerank1("trK_dKOD")
        betU_dKOD = ixp.declarerank2("betU_dKOD", "nosym")
        vetU_dKOD = ixp.declarerank2("vetU_dKOD", "nosym")
        lambdaU_dKOD = ixp.declarerank2("lambdaU_dKOD", "nosym")
        aDD_dKOD = ixp.declarerank3("aDD_dKOD", "sym01")
        hDD_dKOD = ixp.declarerank3("hDD_dKOD", "sym01")
        for k in range(3):
            gaugerhs.alpha_rhs += diss_strength * alpha_dKOD[k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            rhs.cf_rhs += diss_strength * cf_dKOD[k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            rhs.trK_rhs += diss_strength * trK_dKOD[k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            for i in range(3):
                if "2ndOrder" in ShiftCondition:
                    gaugerhs.bet_rhsU[i] += diss_strength * betU_dKOD[i][k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                gaugerhs.vet_rhsU[i] += diss_strength * vetU_dKOD[i][k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                rhs.lambda_rhsU[i] += diss_strength * lambdaU_dKOD[i][k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                for j in range(3):
                    rhs.a_rhsDD[i][j] += diss_strength * aDD_dKOD[i][j][k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                    rhs.h_rhsDD[i][j] += diss_strength * hDD_dKOD[i][j][k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]

    # We use betaU as our upwinding control vector:
    Bq.BSSN_basic_tensors()
    betaU = Bq.betaU

    # END: GENERATE SYMBOLIC EXPRESSIONS
    ######################################

    lhs_names = ["alpha", "cf", "trK"]
    rhs_exprs = [gaugerhs.alpha_rhs, rhs.cf_rhs, rhs.trK_rhs]
    for i in range(3):
        lhs_names.append("betU" + str(i))
        rhs_exprs.append(gaugerhs.bet_rhsU[i])
        lhs_names.append("lambdaU" + str(i))
        rhs_exprs.append(rhs.lambda_rhsU[i])
        lhs_names.append("vetU" + str(i))
        rhs_exprs.append(gaugerhs.vet_rhsU[i])
        for j in range(i, 3):
            lhs_names.append("aDD" + str(i) + str(j))
            rhs_exprs.append(rhs.a_rhsDD[i][j])
            lhs_names.append("hDD" + str(i) + str(j))
            rhs_exprs.append(rhs.h_rhsDD[i][j])

    # Sort the lhss list alphabetically, and rhss to match.
    #   This ensures the RHSs are evaluated in the same order
    #   they're allocated in memory:
    lhs_names, rhs_exprs = [list(x) for x in zip(*sorted(zip(lhs_names, rhs_exprs), key=lambda pair: pair[0]))]

    # Declare the list of lhrh's
    BSSN_RHSs_SymbExpressions = []
    for var in range(len(lhs_names)):
        BSSN_RHSs_SymbExpressions.append(lhrh(lhs=gri.gfaccess("rhs_gfs", lhs_names[var]), rhs=rhs_exprs[var]))

    print_msg_with_timing("BSSN_RHSs", msg="Symbolic", startstop="stop", starttime=starttime)
    return [betaU, BSSN_RHSs_SymbExpressions]


def add_rhs_eval_to_Cfunction_dict(includes=None, rel_path_to_Cparams=os.path.join("."),
                                   enable_rfm_precompute=True, enable_golden_kernels=False,
                                   enable_SIMD=True, enable_split_for_optimizations_doesnt_help=False,
                                   LapseCondition="OnePlusLog", ShiftCondition="GammaDriving2ndOrder_Covariant",
                                   enable_KreissOliger_dissipation=False, enable_stress_energy_source_terms=False,
                                   leave_Ricci_symbolic=True):
    gridsuffix = par.parval_from_str("grid::current_gridsuffix")

    if includes is None:
        includes = []
    if enable_SIMD:
        includes += [os.path.join("SIMD", "SIMD_intrinsics.h")]
    FD_functions_enable = bool(par.parval_from_str("finite_difference::FD_functions_enable"))
    if FD_functions_enable:
        includes += ["finite_difference_functions.h"]

    # Set up the C function for the BSSN RHSs
    desc = "Evaluate the BSSN RHSs"
    name = "rhs_eval" + gridsuffix
    params = "const paramstruct *restrict params, "
    if enable_rfm_precompute:
        params += "const rfm_struct" + gridsuffix + " *restrict rfmstruct, "
    params += """
              const REAL *restrict auxevol_gfs,const REAL *restrict in_gfs,REAL *restrict rhs_gfs"""

    # Construct body:
    betaU, BSSN_RHSs_SymbExpressions = \
        BSSN_RHSs__generate_symbolic_expressions(LapseCondition=LapseCondition, ShiftCondition=ShiftCondition,
                                                 enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
                                                 enable_stress_energy_source_terms=enable_stress_energy_source_terms,
                                                 leave_Ricci_symbolic=leave_Ricci_symbolic)

    starttime = print_msg_with_timing("BSSN_RHSs", msg="Ccodegen", startstop="start")

    FD_outCparams = "outCverbose=False,SIMD_enable=" + str(enable_SIMD)
    FD_outCparams += ",GoldenKernelsEnable=" + str(enable_golden_kernels)
    if gridsuffix != "":
        FD_outCparams += ",gridsuffix=" + gridsuffix

    loopopts = get_loopopts("InteriorPoints", enable_SIMD, enable_rfm_precompute, gridsuffix)
    FDorder = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
    if enable_split_for_optimizations_doesnt_help and FDorder == 6:
        loopopts += ",DisableOpenMP"
        BSSN_RHSs_SymbExpressions_pt1 = []
        BSSN_RHSs_SymbExpressions_pt2 = []
        for lhsrhs in BSSN_RHSs_SymbExpressions:
            if "BETU" in lhsrhs.lhs or "LAMBDAU" in lhsrhs.lhs:
                BSSN_RHSs_SymbExpressions_pt1.append(lhrh(lhs=lhsrhs.lhs, rhs=lhsrhs.rhs))
            else:
                BSSN_RHSs_SymbExpressions_pt2.append(lhrh(lhs=lhsrhs.lhs, rhs=lhsrhs.rhs))
        preloop = """#pragma omp parallel
    {
"""
        preloopbody = fin.FD_outputC("returnstring", BSSN_RHSs_SymbExpressions_pt1,
                                     params=FD_outCparams,
                                     upwindcontrolvec=betaU)
        preloop += "\n#pragma omp for\n" + lp.simple_loop(loopopts, preloopbody)
        preloop += "\n#pragma omp for\n"
        body = fin.FD_outputC("returnstring", BSSN_RHSs_SymbExpressions_pt2,
                              params=FD_outCparams,
                              upwindcontrolvec=betaU)
        postloop = "\n    } // END #pragma omp parallel\n"
    else:
        preloop = ""
        body = fin.FD_outputC("returnstring", BSSN_RHSs_SymbExpressions,
                              params=FD_outCparams,
                              upwindcontrolvec=betaU)
        postloop = ""
    print_msg_with_timing("BSSN_RHSs", msg="Ccodegen", startstop="stop", starttime=starttime)

    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        preloop=preloop, body=body, loopopts=loopopts, postloop=postloop,
        rel_path_to_Cparams=rel_path_to_Cparams)
    return pickle_NRPy_env()


def Ricci__generate_symbolic_expressions():
    ######################################
    # START: GENERATE SYMBOLIC EXPRESSIONS
    starttime = print_msg_with_timing("3-Ricci tensor", msg="Symbolic", startstop="start")

    # Evaluate 3-Ricci tensor:
    import BSSN.BSSN_quantities as Bq
    par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic", "False")

    # We use betaU as our upwinding control vector:
    Bq.BSSN_basic_tensors()

    # Next compute Ricci tensor
    Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
    # END: GENERATE SYMBOLIC EXPRESSIONS
    ######################################
    # Must register RbarDD as gridfunctions, as we're outputting them to gridfunctions here:
    foundit = False
    for i in range(len(gri.glb_gridfcs_list)):
        if "RbarDD00" in gri.glb_gridfcs_list[i].name:
            foundit = True
    if not foundit:
        ixp.register_gridfunctions_for_single_rank2("AUXEVOL", "RbarDD", "sym01")

    Ricci_SymbExpressions = [lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD00"), rhs=Bq.RbarDD[0][0]),
                             lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD01"), rhs=Bq.RbarDD[0][1]),
                             lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD02"), rhs=Bq.RbarDD[0][2]),
                             lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD11"), rhs=Bq.RbarDD[1][1]),
                             lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD12"), rhs=Bq.RbarDD[1][2]),
                             lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD22"), rhs=Bq.RbarDD[2][2])]
    print_msg_with_timing("3-Ricci tensor", msg="Symbolic", startstop="stop", starttime=starttime)

    return Ricci_SymbExpressions


def add_Ricci_eval_to_Cfunction_dict(includes=None, rel_path_to_Cparams=os.path.join("."),
                                     enable_rfm_precompute=True, enable_golden_kernels=False, enable_SIMD=True,
                                     enable_split_for_optimizations_doesnt_help=False):
    gridsuffix = par.parval_from_str("grid::current_gridsuffix")

    if includes is None:
        includes = []
    if enable_SIMD:
        includes += [os.path.join("SIMD", "SIMD_intrinsics.h")]
    FD_functions_enable = bool(par.parval_from_str("finite_difference::FD_functions_enable"))
    if FD_functions_enable:
        includes += ["finite_difference_functions.h"]

    # Set up the C function for the 3-Ricci tensor
    desc = "Evaluate the 3-Ricci tensor"
    name = "Ricci_eval" + gridsuffix
    params = "const paramstruct *restrict params, "
    if enable_rfm_precompute:
        params += "const rfm_struct" + gridsuffix + " *restrict rfmstruct, "
    params += "const REAL *restrict in_gfs, REAL *restrict auxevol_gfs"

    # Construct body:
    Ricci_SymbExpressions = Ricci__generate_symbolic_expressions()
    FD_outCparams = "outCverbose=False,SIMD_enable=" + str(enable_SIMD)
    FD_outCparams += ",GoldenKernelsEnable=" + str(enable_golden_kernels)
    if gridsuffix != "":
        FD_outCparams += ",gridsuffix="+gridsuffix
    starttime = print_msg_with_timing("3-Ricci tensor", msg="Ccodegen", startstop="start")
    loopopts = get_loopopts("InteriorPoints", enable_SIMD, enable_rfm_precompute, gridsuffix)
    preloop = ""
    FDorder = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
    if enable_split_for_optimizations_doesnt_help and FDorder >= 8:
        loopopts += ",DisableOpenMP"
        Ricci_SymbExpressions_pt1 = []
        Ricci_SymbExpressions_pt2 = []
        for lhsrhs in Ricci_SymbExpressions:
            if "RBARDD00" in lhsrhs.lhs or "RBARDD11" in lhsrhs.lhs or "RBARDD22" in lhsrhs.lhs:
                Ricci_SymbExpressions_pt1.append(lhrh(lhs=lhsrhs.lhs, rhs=lhsrhs.rhs))
            else:
                Ricci_SymbExpressions_pt2.append(lhrh(lhs=lhsrhs.lhs, rhs=lhsrhs.rhs))
        preloop = """#pragma omp parallel
    {
#pragma omp for
"""
        preloopbody = fin.FD_outputC("returnstring", Ricci_SymbExpressions_pt1,
                                     params=FD_outCparams)
        preloop += lp.simple_loop(loopopts, preloopbody)
        preloop += "#pragma omp for\n"
        body = fin.FD_outputC("returnstring", Ricci_SymbExpressions_pt2,
                              params=FD_outCparams)
        postloop = "\n    } // END #pragma omp parallel\n"
    else:
        body = fin.FD_outputC("returnstring", Ricci_SymbExpressions,
                              params=FD_outCparams)
        postloop = ""
    print_msg_with_timing("3-Ricci tensor", msg="Ccodegen", startstop="stop", starttime=starttime)

    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        preloop=preloop, body=body, loopopts=loopopts, postloop=postloop,
        rel_path_to_Cparams=rel_path_to_Cparams)
    return pickle_NRPy_env()


def BSSN_constraints__generate_symbolic_expressions(enable_stress_energy_source_terms=False, output_H_only=False):
    ######################################
    # START: GENERATE SYMBOLIC EXPRESSIONS
    starttime = print_msg_with_timing("BSSN constraints", msg="Symbolic", startstop="start")

    # Define the Hamiltonian constraint and output the optimized C code.
    par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic", "True")
    import BSSN.BSSN_constraints as bssncon

    # Returns None if enable_stress_energy_source_terms==False; otherwise returns symb expressions for T4UU
    T4UU = register_stress_energy_source_terms_return_T4UU(enable_stress_energy_source_terms)

    bssncon.BSSN_constraints(add_T4UUmunu_source_terms=False, output_H_only=output_H_only)  # We'll add them below if desired.
    if enable_stress_energy_source_terms:
        import BSSN.BSSN_stress_energy_source_terms as Bsest
        Bsest.BSSN_source_terms_for_BSSN_constraints(T4UU)
        bssncon.H += Bsest.sourceterm_H
        for i in range(3):
            bssncon.MU[i] += Bsest.sourceterm_MU[i]

    BSSN_constraints_SymbExpressions = [lhrh(lhs=gri.gfaccess("aux_gfs", "H"), rhs=bssncon.H)]
    if not output_H_only:
        BSSN_constraints_SymbExpressions += [lhrh(lhs=gri.gfaccess("aux_gfs", "MU0"), rhs=bssncon.MU[0]),
                                             lhrh(lhs=gri.gfaccess("aux_gfs", "MU1"), rhs=bssncon.MU[1]),
                                             lhrh(lhs=gri.gfaccess("aux_gfs", "MU2"), rhs=bssncon.MU[2])]
    par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic", "False")
    print_msg_with_timing("BSSN constraints", msg="Symbolic", startstop="stop", starttime=starttime)
    # END: GENERATE SYMBOLIC EXPRESSIONS
    ######################################
    return BSSN_constraints_SymbExpressions


def add_BSSN_constraints_to_Cfunction_dict(includes=None, rel_path_to_Cparams=os.path.join("."),
                                           enable_rfm_precompute=True, enable_golden_kernels=False, enable_SIMD=True,
                                           enable_stress_energy_source_terms=False,
                                           output_H_only=False):

    gridsuffix = par.parval_from_str("grid::current_gridsuffix")

    if includes is None:
        includes = []
    if enable_SIMD:
        includes += [os.path.join("SIMD", "SIMD_intrinsics.h")]
    FD_functions_enable = bool(par.parval_from_str("finite_difference::FD_functions_enable"))
    if FD_functions_enable:
        includes += ["finite_difference_functions.h"]

    # Set up the C function for the BSSN constraints
    desc = "Evaluate the BSSN constraints"
    name = "BSSN_constraints" + gridsuffix
    params = "const paramstruct *restrict params, "
    if enable_rfm_precompute:
        params += "const rfm_struct" + gridsuffix + " *restrict rfmstruct, "
    params += """
                 const REAL *restrict in_gfs, const REAL *restrict auxevol_gfs, REAL *restrict aux_gfs"""

    # Construct body:
    BSSN_constraints_SymbExpressions = BSSN_constraints__generate_symbolic_expressions(enable_stress_energy_source_terms,
                                                                                       output_H_only=output_H_only)

    FD_outCparams = "outCverbose=False,SIMD_enable=" + str(enable_SIMD)
    FD_outCparams += ",GoldenKernelsEnable=" + str(enable_golden_kernels)
    if gridsuffix != "":
        FD_outCparams += ",gridsuffix="+gridsuffix
    starttime = print_msg_with_timing("BSSN constraints", msg="Ccodegen", startstop="start")
    body = fin.FD_outputC("returnstring", BSSN_constraints_SymbExpressions,
                          params=FD_outCparams)
    print_msg_with_timing("BSSN constraints", msg="Ccodegen", startstop="stop", starttime=starttime)

    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        body=body,
        loopopts=get_loopopts("InteriorPoints", enable_SIMD, enable_rfm_precompute, gridsuffix),
        rel_path_to_Cparams=rel_path_to_Cparams)
    return pickle_NRPy_env()


def add_enforce_detgammahat_constraint_to_Cfunction_dict(includes=None, rel_path_to_Cparams=os.path.join("."),
                                                         enable_rfm_precompute=True, enable_golden_kernels=False):
    # This function disables SIMD, as it includes cbrt() and abs() functions.
    gridsuffix = par.parval_from_str("grid::current_gridsuffix")

    if includes is None:
        includes = []
    FD_functions_enable = bool(par.parval_from_str("finite_difference::FD_functions_enable"))
    if FD_functions_enable:
        includes += ["finite_difference_functions.h"]

    # Set up the C function for enforcing the det(gammabar) = det(gammahat) BSSN algebraic constraint
    desc = "Enforce the det(gammabar) = det(gammahat) (algebraic) constraint"
    name = "enforce_detgammahat_constraint" + gridsuffix
    params = "const paramstruct *restrict params, "
    if enable_rfm_precompute:
        params += "const rfm_struct" + gridsuffix + " *restrict rfmstruct, "
    params += "REAL *restrict in_gfs"

    # Construct body:
    enforce_detg_constraint_symb_expressions = EGC.Enforce_Detgammahat_Constraint_symb_expressions()

    FD_outCparams = "outCverbose=False,SIMD_enable=False"
    FD_outCparams += ",GoldenKernelsEnable=" + str(enable_golden_kernels)
    if gridsuffix != "":
        FD_outCparams += ",gridsuffix="+gridsuffix
    starttime = print_msg_with_timing("Enforcing det(gammabar)=det(gammahat) constraint", msg="Ccodegen", startstop="start")
    body = fin.FD_outputC("returnstring", enforce_detg_constraint_symb_expressions,
                          params=FD_outCparams)
    print_msg_with_timing("Enforcing det(gammabar)=det(gammahat) constraint", msg="Ccodegen", startstop="stop", starttime=starttime)

    enable_SIMD = False
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        body=body,
        loopopts=get_loopopts("AllPoints", enable_SIMD, enable_rfm_precompute, gridsuffix),
        rel_path_to_Cparams=rel_path_to_Cparams)
    return pickle_NRPy_env()
