# hermite_interpolator.py:
#  As documented in the NRPy+ tutorial notebook:
#    Tutorial-Hermite_Interpolator.ipynb,
#  This module generates C kernels for numerically
#   interpolating for general uniform grids
#
# Depends primarily on: outputC.py and grid.py.

# Author: Maria C. Babiuc Hamilton template courtesy Zachariah B. Etienne
#         babiuc **at** marshall **dot* edu

from outputC import parse_outCparams_string, outC_function_dict, outC_function_prototype_dict, outC_NRPy_basic_defines_h_dict, outC_function_master_list  # NRPy+: Core C code output module
import NRPy_param_funcs as par   # NRPy+: parameter interface
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import grid as gri               # NRPy+: Functions having to do with numerical grids
import os, sys                   # Standard Python module for multiplatform OS-level functions
from hermite_interpolator_helpers import extract_from_list_of_interp_vars__base_gfs_and_interp_ops_lists
from hermite_interpolator_helpers import generate_list_of_interp_vars_from_lhrh_sympyexpr_list
from hermite_interpolator_helpers import read_gfs_from_memory, HIparams, construct_Ccode

# Step 1: Initialize free parameters for this module:
modulename = __name__
# Hermite interpolator dimension order
par.initialize_param(par.glb_param("int",  modulename, "HI_DIMENSIONS_ORDER",          3))
par.initialize_param(par.glb_param("bool", modulename, "enable_HI_functions",      False))

def HI_outputC(filename, sympyexpr_list, params="", upwindcontrolvec=""):
    outCparams = parse_outCparams_string(params)

    # Step 0.a:
    # In case sympyexpr_list is a single sympy expression,
    #     convert it to a list with just one element.
    #     This enables the rest of the routine to assume
    #     sympyexpr_list is indeed a list.
    if not isinstance(sympyexpr_list, list):
        sympyexpr_list = [sympyexpr_list]

    # Step 0.b:
    # hermite_interpolator.py takes control over outCparams.includebraces here,
    #     which is necessary because outputC() is called twice:
    #     first for the reads from main memory and interpolator dimension order,
    #     and second for the SymPy expressions, and writes to main memory.
    # If outCparams.includebraces==True, then it will close off the braces
    #     after the finite difference stencil expressions and start new ones
    #     for the SymPy expressions and writes to main memory, resulting
    #     in a non-functioning C code.
    # To get around this issue, we create braces around the entire
    #     string of C output from this function, only if
    #     outCparams.includebraces==True.
    # See Step 5 for open and close braces
    if outCparams.includebraces == "True":
        indent = "  "
    else:
        indent = ""

    # Step 0.c: HIparams named tuple stores parameters used in the codegen of the Hermite interpolator
    HIparams.enable_SIMD         = outCparams.enable_SIMD
    HIparams.PRECISION           = par.parval_from_str("PRECISION")
    HIparams.HI_CD_order         = par.parval_from_str("HI_DIMENSIONS_ORDER")
    HIparams.enable_HI_functions = par.parval_from_str("enable_HI_functions")
    HIparams.DIM                 = par.parval_from_str("DIM")
    HIparams.MemAllocStyle       = par.parval_from_str("MemAllocStyle")
    HIparams.upwindcontrolvec    = upwindcontrolvec
    HIparams.fullindent          = indent + outCparams.preindent
    HIparams.outCparams          = params

    # Step 1: Generate from list of SymPy expressions in the form
    #     [lhrh(lhs=var, rhs=expr),lhrh(...),...]
    #     all interpolation expressions, which we will process next.
    list_of_interp_vars = generate_list_of_interp_vars_from_lhrh_sympyexpr_list(sympyexpr_list, HIparams)

    # Step 2a: Extract from list_of_interp_vars a list of base gridfunctions
    #         and a list of interpolation  operators. Usually takes list of SymPy
    #         symbols as input, but could just take strings, as this function
    #         does only string manipulations.
    # Example:
    # >>> extract_from_list_of_interp_vars__base_gfs_and_interp_ops_lists(["aDD_dD012","hDD_dDD0112"])
    # (['aDD01', 'aDD01', 'vetU2', 'hDD01'], ['dD2', 'dDD12'])
    list_of_base_gridfunction_names_in_interp, list_of_interp_operators = \
        extract_from_list_of_interp_vars__base_gfs_and_interp_ops_lists(list_of_interp_vars)

    # Step 2b:
    # Next, check each base gridfunction to determine whether
    #     it is indeed registered as a gridfunction.
    #     If not, exit with error.
    for basegf in list_of_base_gridfunction_names_in_interp:
        is_gf = False
        for gf in gri.glb_gridfcs_list:
            if basegf == str(gf.name):
                is_gf = True
        if not is_gf:
            print("Error: Attempting to apply the interpolation of "+basegf+", which is not a registered gridfunction.")
            print("       Make sure your gridfunction name does not have any underscores in it!")
            sys.exit(1)

    # Step 2c:
    # Check each interpolation operator to make sure it is
    #     supported. If not, error out.
    for i in range(len(list_of_interp_operators)):
        found_interpID = False
        for interpID in ["dD", "dupD", "ddnD"]:
            if interpID in list_of_interp_operators[i]:
                found_interpID = True
        if not found_interpID:
            print("Error: Valid interpolator operator in "+list_of_interp_operators[i]+" not found.")
            sys.exit(1)

    # Step 3:
    # Evaluate the interpolator stencil for each
    #     interpolation operator, being careful not to
    #     needlessly recompute.
    # Note: Each interpolator stencil consists
    #     of two parts:
    #     1) The coefficient, and
    #     2) The index corresponding to the coefficient.
    #     The former is stored as a rational number, and
    #     the latter as a simple string, such that e.g.,
    #     in 3D, the empty string corresponds to (i,j,k),
    #     the string "ip1" corresponds to (i+1,j,k),
    #     the string "ip1kp1" corresponds to (i+1,j,k+1),
    #     etc.
    hicoeffs = [[] for i in range(len(list_of_interp_operators))]
    histencl = [[[] for i in range(4)] for j in range(len(list_of_interp_operators))]
    for i in range(len(list_of_interp_operators)):
        hicoeffs[i], histencl[i] = compute_hicoeffs_histencl(list_of_interp_operators[i])

    # Step 4: Create C code to read gridfunctions from memory
    read_from_memory_Ccode = read_gfs_from_memory(list_of_base_gridfunction_names_in_interp, histencl, sympyexpr_list,
                                                  HIparams)

    # Step 5: construct C code.
    Coutput = ""
    if outCparams.includebraces == "True":
        Coutput = outCparams.preindent + "{\n"
    Coutput = construct_Ccode(sympyexpr_list, list_of_interp_vars,
                           list_of_base_gridfunction_names_in_interp, list_of_interp_operators,
                           hicoeffs, histencl, read_from_memory_Ccode, HIparams, Coutput)
    if outCparams.includebraces == "True":
        Coutput += outCparams.preindent+"}"

    # Step 6: Output the C code in desired format: stdout, string, or file.
    if filename == "stdout":
        print(Coutput)
    elif filename == "returnstring":
        return Coutput
    else:
        # Output to the file specified by outCfilename
        with open(filename, outCparams.outCfileaccess) as file:
            file.write(Coutput)
        successstr = ""
        if outCparams.outCfileaccess == "a":
            successstr = "Appended "
        elif outCparams.outCfileaccess == "w":
            successstr = "Wrote "
        print(successstr + "to file \"" + filename + "\"")

################
# TO BE DEPRECATED:
def output_hermite_interpolator_functions_h(path=os.path.join(".")):
    with open(os.path.join(path, "hermite_interpolator_functions.h"), "w") as file:
        file.write("""
#ifndef __HI_FUNCTIONS_H__
#define __HI_FUNCTIONS_H__
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
""")
        UNUSED   = "__attribute__((unused))"
        NOINLINE = "__attribute__((noinline))"
        if par.parval_from_str("grid::GridFuncMemAccess") == "ETK":
            UNUSED   = "CCTK_ATTRIBUTE_UNUSED"
            NOINLINE = "CCTK_ATTRIBUTE_NOINLINE"
        file.write("#define _UNUSED   " + UNUSED   + "\n")
        file.write("#define _NOINLINE " + NOINLINE + "\n")

        for key, item in outC_function_dict.items():
            if "__HI_OPERATOR_FUNC__" in item:
                file.write(item.replace("const REAL_SIMD_ARRAY _NegativeOne_ =",
                                        "const REAL_SIMD_ARRAY "+UNUSED+" _NegativeOne_ =")) # Many of the NegativeOne's get optimized away in the SIMD postprocessing step. No need for all the warnings

        # Clear all HI functions from outC_function_dict after outputting to hermite_interpolator_functions.h.
        #   Otherwise outputC will be outputting these as separate individual C codes & attempting to build them in Makefile.
        key_list_del = []
        element_del = []
        for i, func in enumerate(outC_function_master_list):
            if "__HI_OPERATOR_FUNC__" in func.desc:
                if func.name not in key_list_del:
                    key_list_del += [func.name]
                if func not in element_del:
                    element_del += [func]
        for func in element_del:
            outC_function_master_list.remove(func)
        for key in key_list_del:
            outC_function_dict.pop(key)
            if key in outC_function_prototype_dict:
                outC_function_prototype_dict.pop(key)
        file.write("#endif // #ifndef __HI_FUNCTIONS_H__\n")
################


def register_C_functions_and_NRPy_basic_defines(NGHOSTS_account_for_onezone_upwind=False, enable_SIMD=True):
    # First register C functions needed by finite_difference

    # Then set up the dictionary entry for hermite_interpolator in NRPy_basic_defines
    NGHOSTS = int(par.parval_from_str("hermite_interpolator::HI_DIMENSIONS_ORDER")/2)
    if NGHOSTS_account_for_onezone_upwind:
        NGHOSTS += 1
    Nbd_str = """
// Set the number of ghost zones
// Note that upwinding in e.g., BSSN requires that NGHOSTS = HI_DIMENSIONS_ORDER/2 + 1 <- Notice the +1.
"""
    Nbd_str += "#define NGHOSTS " + str(NGHOSTS)+"\n"
    if not enable_SIMD:
        Nbd_str += """
// When enable_SIMD = False, this is the UPWIND_ALG() macro:
#define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0\n"""
    outC_NRPy_basic_defines_h_dict["finite_difference"] = Nbd_str


#######################################################
#  HERMITE INTERPOLATOR COEFFICIENT ALGORITHM

#  Define the to-be-inverted matrix, A.
#  We define A row-by-row, via the following pattern
#  that applies for arbitrary order.
#
#  As an example, consider a 2-D hermite interpolator
#  where we wish to interpolate x between xk and xk+h, at point s h
#
#  Then A is given by:
#
#  The general pattern is as follows:
#
#  1) 
#  2) 
#  3) 
def setup_HI_matrix__return_inverse(STENCILWIDTH, UPDOWNWIND_stencil_shift):
    # Set up matrix based on the stencil size (HIORDER+1).
    #         See documentation above for details on how this
    #         matrix is set up.
    M = sp.zeros(STENCILWIDTH, STENCILWIDTH)
    for i in range(STENCILWIDTH):
        for j in range(STENCILWIDTH):
            if i == 0:
                M[(i, j)] = 1  # Setting n^0 = 1 for all n, including n=0, because this matches the pattern
            else:
                dist_from_xeq0_col = j - sp.Rational((STENCILWIDTH - 1), 2) + UPDOWNWIND_stencil_shift
                if dist_from_xeq0_col == 0:
                    M[(i, j)] = 0
                else:
                    M[(i, j)] = dist_from_xeq0_col ** i
    return M**(-1)


def compute_hicoeffs_histencl(interpstring, HIORDER=-1):
    # Step 0: Set hermite interpolator order, stencil size, and up/downwinding
    if HIORDER == -1:
        HIORDER = par.parval_from_str("HI_DIMENSIONS_ORDER")

    STENCILWIDTH = HIORDER+1
    UPDOWNWIND_stencil_shift = 0
    # dup/dnD = single-point-offset upwind/downwinding.
    if "dupD" in interpstring:
        UPDOWNWIND_stencil_shift =  1
    elif "ddnD" in interpstring:
        UPDOWNWIND_stencil_shift = -1
    # dfullup/dnD = full upwind/downwinding.
    elif "dfullupD" in interpstring:
        UPDOWNWIND_stencil_shift =  int(HIORDER/2)
    elif "dfulldnD" in interpstring:
        UPDOWNWIND_stencil_shift = -int(HIORDER/2)

    # Step 1: Set up HI matrix and return the inverse, as documented above.
    Minv = setup_HI_matrix__return_inverse(STENCILWIDTH, UPDOWNWIND_stencil_shift)

    # Step 2:
    #     Based on the input interpolator string,
    #     pick out the relevant row of the matrix
    #     inverse, as outlined in the detailed code
    #     comments prior to this function definition.
    interptype = "FirstInterp"
    matrixrow = 1
    if "DDD" in interpstring:
        print("Error: Only 2D interpolation currently supported.")
        print("       Feel free to contribute to NRPy+ to extend its functionality!")
        sys.exit(1)
    elif "DD" in interpstring:

        if interpstring[len(interpstring)-1] == interpstring[len(interpstring)-2]:
            # Assuming i==j, we call "SecondInterp":
            interptype = "SecondInterp"
            matrixrow = 2
        else:
            # Assuming i!=j, we call a MIXED second interpolator,
            #     which is computed using a composite of first interpolator operations.
            interptype = "MixedSecondInterp"
    else:
        pass

    # Step 3:
    #     Set interpolation  coefficients
    #     and stencil points corresponding to
    #     each interpolation coefficient.
    hicoeffs = []
    histencl = []
    if interptype != "MixedSecondInterp":
        for i in range(STENCILWIDTH):
            idx4 = [0, 0, 0, 0]
            # First compute interpolation coefficient.
            hicoeff = sp.factorial(matrixrow)*Minv[(i, matrixrow)]
            # Do not store hicoeff or histencil if
            # interpolation coefficient is zero.
            if hicoeff != 0:
                hicoeffs.append(hicoeff)

                # Next store interpolation stencil point
                # corresponding to coefficient.
                gridpt_posn = i - int((STENCILWIDTH-1)/2) + UPDOWNWIND_stencil_shift
                if gridpt_posn != 0:
                    dirn = int(interpstring[len(interpstring)-1])
                    idx4[dirn] = gridpt_posn
                histencl.append(idx4)
    else:
        # Mixed second derivative interpolation coeffs
        #     consist of products of first deriv coeffs,
        #     defined in first Minv matrix row.
        for i in range(STENCILWIDTH):
            for j in range(STENCILWIDTH):
                idx4 = [0, 0, 0, 0]

                # First compute interpolation coefficient.
                hicoeff = (sp.factorial(matrixrow)*Minv[(i, matrixrow)]) * \
                          (sp.factorial(matrixrow)*Minv[(j, matrixrow)])

                # Do not store hicoeff or histencil if
                # interpolator coefficient is zero.
                if hicoeff != 0:
                    hicoeffs.append(hicoeff)

                    # Next store interpolator stencil point
                    # corresponding to coefficient.
                    gridpt_posn1 = i - int((STENCILWIDTH - 1) / 2)
                    gridpt_posn2 = j - int((STENCILWIDTH - 1) / 2)
                    dirn1 = int(interpstring[len(interpstring) - 1])
                    dirn2 = int(interpstring[len(interpstring) - 2])
                    idx4[dirn1] = gridpt_posn1
                    idx4[dirn2] = gridpt_posn2
                    histencl.append(idx4)
    return hicoeffs, histencl
