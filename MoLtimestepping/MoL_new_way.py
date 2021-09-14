# As documented in the NRPy+ tutorial module
#   Tutorial-RK_Butcher_Table_Generating_C_Code_new_way.ipynb,
#   this module will produce the required C codes for
#   allocating required memory Method of Lines (MoL) timestepping,
#   implementing MoL timestepping, and deallocating memory

# Authors: Brandon Clark
#          Zachariah B. Etienne (maintainer)
#          zachetie **at** gmail **dot* com

import sympy as sp  # Import SymPy, a computer algebra system written entirely in Python
import os, sys      # Standard Python modules for multiplatform OS-level functions
from MoLtimestepping.RK_Butcher_Table_Dictionary import Butcher_dict
from outputC import add_to_Cfunction_dict, indent_Ccode  # NRPy+: Basic C code output functionality


# Check if Butcher Table is diagonal
def diagonal(key):
    Butcher = Butcher_dict[key][0]
    L = len(Butcher)-1  # Establish the number of rows to check for diagonal trait, all bust last row
    row_idx = 0  # Initialize the Butcher table row index
    for i in range(L):  # Check all the desired rows
        for j in range(1,row_idx):  # Check each element before the diagonal element in a row
            if Butcher[i][j] != sp.sympify(0):  # If any non-diagonal coeffcient is non-zero,
                                                # then the table is not diagonal
                return False
        row_idx += 1  # Update to check the next row
    return True


# Each MoL method has its own set of names for groups of gridfunctions,
#   aiming to be sufficiently descriptive. So for example a set of
#   gridfunctions that store "k_1" in an RK-like method could be called
#   "k1_gfs".
def generate_gridfunction_names(MoL_method = "RK4"):
    # Step 3.a: MoL gridfunctions fall into 3 overlapping categories:
    #           1) y_n=y_i(t_n) gridfunctions y_n_gfs, which stores data for the vector of gridfunctions y_i at t_n,
    #              the start of each MoL timestep.
    #           2) non-y_n gridfunctions, needed to compute the data at t_{n+1}. Often labeled with k_i in the name,
    #              these gridfunctions are *not* needed at the start of each timestep, so are available for temporary
    #              storage when gridfunctions needed for diagnostics are computed at the start of each timestep.
    #              These gridfunctions can also be freed during a regrid, to enable storage for the post-regrid
    #              destination y_n_gfs.
    #           3) Diagnostic output gridfunctions diagnostic_output_gfs, which simply uses the memory from auxiliary
    #              gridfunctions at one auxiliary time to compute diagnostics at t_n.

    # Here we specify which gridfunctions fall into each category, starting with the obvious: y_n_gridfunctions
    y_n_gridfunctions = "y_n_gfs"

    # Next the less-obvious, which depend on non-y_n_gfs
    non_y_n_gridfunctions_list = []

    # No matter the method we define gridfunctions "y_n_gfs" to store the initial data
    if diagonal(MoL_method) and "RK3" in MoL_method:
        non_y_n_gridfunctions_list.append("k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs")
        non_y_n_gridfunctions_list.append("k2_or_y_nplus_a32_k2_gfs")
        diagnostic_gridfunctions_point_to = "k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs"
    else:
        if not diagonal(MoL_method):  # Allocate memory for non-diagonal Butcher tables
            # Determine the number of k_i steps based on length of Butcher Table
            num_k = len(Butcher_dict[MoL_method][0])-1
            # For non-diagonal tables an intermediate gridfunction "next_y_input" is used for rhs evaluations
            non_y_n_gridfunctions_list.append("next_y_input_gfs")
            for i in range(num_k): # Need to allocate all k_i steps for a given method
                non_y_n_gridfunctions_list.append("k" + str(i + 1) + "_gfs")
            diagnostic_gridfunctions_point_to = "k1_gfs"
        else:  # Allocate memory for diagonal Butcher tables, which use a "y_nplus1_running_total gridfunction"
            non_y_n_gridfunctions_list.append("y_nplus1_running_total_gfs")
            if MoL_method != 'Euler':  # Allocate memory for diagonal Butcher tables that aren't Euler
                # Need k_odd for k_1,3,5... and k_even for k_2,4,6...
                non_y_n_gridfunctions_list.append("k_odd_gfs")
                non_y_n_gridfunctions_list.append("k_even_gfs")
            diagnostic_gridfunctions_point_to = "y_nplus1_running_total_gfs"
    non_y_n_gridfunctions_list.append("auxevol_gfs")

    return y_n_gridfunctions, non_y_n_gridfunctions_list, diagnostic_gridfunctions_point_to


# add_to_Cfunction_dict_MoL_malloc() registers
#           MoL_malloc_y_n_gfs() and
#           MoL_malloc_non_y_n_gfs(), which allocate memory for
#           the indicated sets of gridfunctions
def add_to_Cfunction_dict_MoL_malloc(MoL_method, which_gfs):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc  = "Method of Lines (MoL) for \"" + MoL_method + "\" method: Allocate memory for \""+which_gfs+"\" gridfunctions\n"
    desc += "   * y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep\n"
    desc += "   * non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method\n"
    c_type = "void"

    y_n_gridfunctions, non_y_n_gridfunctions_list, diagnostic_gridfunctions_point_to = \
        generate_gridfunction_names(MoL_method = MoL_method)

    gridfunctions_list = []
    if which_gfs == "y_n_gfs":
        gridfunctions_list = [y_n_gridfunctions]
    elif which_gfs == "non_y_n_gfs":
        gridfunctions_list = non_y_n_gridfunctions_list
    else:
        print("ERROR: which_gfs = \"" + which_gfs + "\" unrecognized.")
        sys.exit(1)
    name = "MoL_malloc_" + which_gfs
    params = "const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs"
    body = "const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;\n"
    for gridfunctions in gridfunctions_list:
        num_gfs = "NUM_EVOL_GFS"
        if gridfunctions == "auxevol_gfs":
            num_gfs = "NUM_AUXEVOL_GFS"
        body += "gridfuncs." + gridfunctions + " = (REAL *restrict)malloc(sizeof(REAL) * " + num_gfs + " * Nxx_plus_2NGHOSTS_tot);\n"
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=indent_Ccode(body, "  "),
        rel_path_to_Cparams=os.path.join("."))

# add_to_Cfunction_dict_MoL_free_memory() registers
#           MoL_free_memory_y_n_gfs() and
#           MoL_free_memory_non_y_n_gfs(), which free memory for
#           the indicated sets of gridfunctions
def add_to_Cfunction_dict_MoL_free_memory(MoL_method, which_gfs):
        includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
        desc = "Method of Lines (MoL) for \"" + MoL_method + "\" method: Free memory for \"" + which_gfs + "\" gridfunctions\n"
        desc += "   - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep\n"
        desc += "   - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method\n"
        c_type = "void"

        y_n_gridfunctions, non_y_n_gridfunctions_list, diagnostic_gridfunctions_point_to = \
            generate_gridfunction_names(MoL_method=MoL_method)

        gridfunctions_list = []
        if which_gfs == "y_n_gfs":
            gridfunctions_list = [y_n_gridfunctions]
        elif which_gfs == "non_y_n_gfs":
            gridfunctions_list = non_y_n_gridfunctions_list
        else:
            print("ERROR: which_gfs = \"" + which_gfs + "\" unrecognized.")
            sys.exit(1)
        name = "MoL_free_memory_" + which_gfs
        params = "const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs"
        body = ""
        for gridfunctions in gridfunctions_list:
            body += "    free(gridfuncs." + gridfunctions + ");\n"
        add_to_Cfunction_dict(
            includes=includes,
            desc=desc,
            c_type=c_type, name=name, params=params,
            body=indent_Ccode(body, "  "),
            rel_path_to_Cparams=os.path.join("."))


# single_RK_substep() performs necessary replacements to
#   define C code for a single RK substep
#   (e.g., computing k_1 and then updating the outer boundaries)
def single_RK_substep(commentblock, RHS_str, RHS_input_str, RHS_output_str, RK_lhss_list, RK_rhss_list,
                      post_RHS_list, post_RHS_output_list, indent="  "):
    addl_indent = ""
    return_str  = commentblock + "\n"
    if not isinstance(RK_lhss_list, list):
        RK_lhss_list = [RK_lhss_list]
    if not isinstance(RK_rhss_list, list):
        RK_rhss_list = [RK_rhss_list]

    if not isinstance(post_RHS_list, list):
        post_RHS_list = [post_RHS_list]
    if not isinstance(post_RHS_output_list, list):
        post_RHS_output_list = [post_RHS_output_list]

    # Part 1: RHS evaluation:
    return_str += indent_Ccode(RHS_str.replace("RK_INPUT_GFS",  RHS_input_str).
                                       replace("RK_OUTPUT_GFS", RHS_output_str)+"\n", indent=addl_indent)

    # Part 2: RK update
    return_str += addl_indent + "LOOP_ALL_GFS_GPS(i) {\n"
    for lhs, rhs in zip(RK_lhss_list, RK_rhss_list):
        return_str += addl_indent + indent + lhs + "[i] = " + rhs + ";\n"
    return_str += addl_indent + "}\n"

    # Part 3: Call post-RHS functions
    for post_RHS, post_RHS_output in zip(post_RHS_list, post_RHS_output_list):
        return_str += indent_Ccode(post_RHS.replace("RK_OUTPUT_GFS", post_RHS_output), indent=addl_indent)

    return return_str


########################################################################################################################
# EXAMPLE
# ODE: y' = f(t,y), y(t_0) = y_0
# Starting at time t_n with solution having value y_n and trying to update to y_nplus1 with timestep dt

# Example of scheme for RK4 with k_1, k_2, k_3, k_4 (Using non-diagonal algorithm) Notice this requires storage of
# y_n, y_nplus1, k_1 through k_4

# k_1      = dt*f(t_n, y_n)
# k_2      = dt*f(t_n + 1/2*dt, y_n + 1/2*k_1)
# k_3      = dt*f(t_n + 1/2*dt, y_n + 1/2*k_2)
# k_4      = dt*f(t_n + dt, y_n + k_3)
# y_nplus1 = y_n + 1/3k_1 + 1/6k_2 + 1/6k_3 + 1/3k_4

# Example of scheme RK4 using only k_odd and k_even (Diagonal algroithm) Notice that this only requires storage

# k_odd     = dt*f(t_n, y_n)
# y_nplus1  = 1/3*k_odd
# k_even    = dt*f(t_n + 1/2*dt, y_n + 1/2*k_odd)
# y_nplus1 += 1/6*k_even
# k_odd     = dt*f(t_n + 1/2*dt, y_n + 1/2*k_even)
# y_nplus1 += 1/6*k_odd
# k_even    = dt*f(t_n + dt, y_n + k_odd)
# y_nplus1 += 1/3*k_even
########################################################################################################################
def add_to_Cfunction_dict_MoL_step_forward_one_timestep(MoL_method, RHS_string = "", post_RHS_string = ""):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc  = "Method of Lines (MoL) for \"" + MoL_method + "\" method: Step forward one full timestep.\n"
    c_type = "void"
    name = "MoL_step_forward_one_timestep"
    params = "const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs"

    indent = ""  # We don't bother with an indent here.

    body = indent + "// C code implementation of -={ " + MoL_method + " }=- Method of Lines timestepping.\n\n"

    y_n_gridfunctions, non_y_n_gridfunctions_list, _throwaway = generate_gridfunction_names(MoL_method)
    body += "// First set gridfunction aliases from gridfuncs struct\n\n"
    body += "// y_n gridfunctions:\n"
    body += "REAL restrict *" + y_n_gridfunctions + " = gridfuncs." + y_n_gridfunctions + ";\n"
    body += "\n"
    body += "// Temporary timelevel & AUXEVOL gridfunctions:\n"
    for gf in non_y_n_gridfunctions_list:
        body += "REAL restrict *" + gf + " = gridfuncs." + gf + ";\n"
    body += "\n"
    body += "// Next perform a full step forward in time\n"

    # Implement Method of Lines (MoL) Timestepping
    Butcher = Butcher_dict[MoL_method][0] # Get the desired Butcher table from the dictionary
    num_steps = len(Butcher)-1 # Specify the number of required steps to update solution

    # Diagonal RK3 only!!!
    if diagonal(MoL_method) and "RK3" in MoL_method:
        #  In a diagonal RK3 method, only 3 gridfunctions need be defined. Below implements this approach.

        # k_1
        body += """
// In a diagonal RK3 method like this one, only 3 gridfunctions need be defined. Below implements this approach.
// Using y_n_gfs as input, k1 and apply boundary conditions\n"""

        body += single_RK_substep(
            commentblock = """// -={ START k1 substep }=-
// RHS evaluation:
//  1. We will store k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs now as
//     ...  the update for the next rhs evaluation y_n + a21*k1*dt
// Post-RHS evaluation:
//  1. Apply post-RHS to y_n + a21*k1*dt""",
            RHS_str = RHS_string,
            RHS_input_str = "y_n_gfs", RHS_output_str = "k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs",
            RK_lhss_list = ["k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs"],
            RK_rhss_list = ["("+sp.ccode(Butcher[1][1]).replace("L","")+")*k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i]*dt + y_n_gfs[i]"],
            post_RHS_list = [post_RHS_string], post_RHS_output_list = ["k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs"]) + "// -={ END k1 substep }=-\n\n"

        # k_2
        body += single_RK_substep(
            commentblock="""// -={ START k2 substep }=-
// RHS evaluation:
//    1. Reassign k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs to be the running total y_{n+1}; a32*k2*dt to the running total
//    2. Store k2_or_y_nplus_a32_k2_gfs now as y_n + a32*k2*dt
// Post-RHS evaluation:
//    1. Apply post-RHS to both y_n + a32*k2 (stored in k2_or_y_nplus_a32_k2_gfs)
//       ... and the y_{n+1} running total, as they have not been applied yet to k2-related gridfunctions""",
            RHS_str=RHS_string,
            RHS_input_str="k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs", RHS_output_str="k2_or_y_nplus_a32_k2_gfs",
            RK_lhss_list=["k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs","k2_or_y_nplus_a32_k2_gfs"],
            RK_rhss_list=["("+sp.ccode(Butcher[3][1]).replace("L","")+")*(k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i] - y_n_gfs[i])/("+sp.ccode(Butcher[1][1]).replace("L","")+") + y_n_gfs[i] + ("+sp.ccode(Butcher[3][2]).replace("L","")+")*k2_or_y_nplus_a32_k2_gfs[i]*dt",
                          "("+sp.ccode(Butcher[2][2]).replace("L","")+")*k2_or_y_nplus_a32_k2_gfs[i]*dt + y_n_gfs[i]"],
            post_RHS_list=[post_RHS_string,post_RHS_string],
            post_RHS_output_list=["k2_or_y_nplus_a32_k2_gfs","k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs"]) + "// -={ END k2 substep }=-\n\n"

        # k_3
        body += single_RK_substep(
            commentblock="""// -={ START k3 substep }=-
// RHS evaluation:
//    1. Add k3 to the running total and save to y_n
// Post-RHS evaluation:
//    1. Apply post-RHS to y_n""",
            RHS_str=RHS_string,
            RHS_input_str="k2_or_y_nplus_a32_k2_gfs", RHS_output_str="y_n_gfs",
            RK_lhss_list=["y_n_gfs","k2_or_y_nplus_a32_k2_gfs"],
            RK_rhss_list=["k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i] + ("+sp.ccode(Butcher[3][3]).replace("L","")+")*y_n_gfs[i]*dt"],
            post_RHS_list=[post_RHS_string],
            post_RHS_output_list=["y_n_gfs"]) + "// -={ END k3 substep }=-\n\n"
    else:
        y_n = "y_n_gfs"
        if not diagonal(MoL_method):
            for s in range(num_steps):
                next_y_input = "next_y_input_gfs"

                # If we're on the first step (s=0), we use y_n gridfunction as input.
                #      Otherwise next_y_input is input. Output is just the reverse.
                if s == 0:  # If on first step:
                    RHS_input = y_n
                else:    # If on second step or later:
                    RHS_input = next_y_input
                RHS_output = "k" + str(s + 1) + "_gfs"
                if s == num_steps-1: # If on final step:
                    RK_lhs = y_n
                    RK_rhs = y_n + "[i] + dt*("
                else:                # If on anything but the final step:
                    RK_lhs = next_y_input
                    RK_rhs = y_n + "[i] + dt*("
                for m in range(s+1):
                    if Butcher[s+1][m+1] != 0:
                        if Butcher[s+1][m+1] != 1:
                            RK_rhs += " + k"+str(m+1)+"_gfs[i]*("+sp.ccode(Butcher[s+1][m+1]).replace("L","")+")"
                        else:
                            RK_rhs += " + k"+str(m+1)+"_gfs[i]"
                RK_rhs += " )"

                post_RHS = post_RHS_string
                if s == num_steps-1: # If on final step:
                    post_RHS_output = y_n
                else:                # If on anything but the final step:
                    post_RHS_output = next_y_input

                body += single_RK_substep(
                    commentblock="// -={ START k" + str(s + 1) + " substep }=-",
                    RHS_str=RHS_string,
                    RHS_input_str=RHS_input, RHS_output_str=RHS_output,
                    RK_lhss_list=[RK_lhs],   RK_rhss_list=[RK_rhs],
                    post_RHS_list=[post_RHS],
                    post_RHS_output_list=[post_RHS_output]) + "// -={ END k" + str(s + 1) + " substep }=-\n\n"
        else:  # diagonal case:
            y_nplus1_running_total = "y_nplus1_running_total_gfs"
            if MoL_method == 'Euler': # Euler's method doesn't require any k_i, and gets its own unique algorithm
                body += single_RK_substep(
                    commentblock=indent + "// ***Euler timestepping only requires one RHS evaluation***",
                    RHS_str=RHS_string,
                    RHS_input_str=y_n, RHS_output_str=y_nplus1_running_total,
                    RK_lhss_list=[y_n],   RK_rhss_list=[y_n+"[i] + "+y_nplus1_running_total+"[i]*dt"],
                    post_RHS_list=[post_RHS_string],
                    post_RHS_output_list=[y_n])
            else:
                for s in range(num_steps):
                    # If we're on the first step (s=0), we use y_n gridfunction as input.
                    # and k_odd as output.
                    if s == 0:
                        RHS_input  = "y_n_gfs"
                        RHS_output = "k_odd_gfs"
                    # For the remaining steps the inputs and ouputs alternate between k_odd and k_even
                    elif s % 2 == 0:
                        RHS_input = "k_even_gfs"
                        RHS_output = "k_odd_gfs"
                    else:
                        RHS_input = "k_odd_gfs"
                        RHS_output = "k_even_gfs"

                    RK_lhs_list = []
                    RK_rhs_list = []
                    if s != num_steps-1:  # For anything besides the final step
                        if s == 0:  # The first RK step
                            RK_lhs_list.append(y_nplus1_running_total)
                            RK_rhs_list.append(RHS_output+"[i]*dt*("+sp.ccode(Butcher[num_steps][s+1]).replace("L","")+")")

                            RK_lhs_list.append(RHS_output)
                            RK_rhs_list.append(y_n+"[i] + "+RHS_output+"[i]*dt*("+sp.ccode(Butcher[s+1][s+1]).replace("L","")+")")
                        else:
                            if Butcher[num_steps][s+1] != 0:
                                RK_lhs_list.append(y_nplus1_running_total)
                                if Butcher[num_steps][s+1] != 1:
                                    RK_rhs_list.append(y_nplus1_running_total+"[i] + "+RHS_output+"[i]*dt*("+sp.ccode(Butcher[num_steps][s+1]).replace("L","")+")")
                                else:
                                    RK_rhs_list.append(y_nplus1_running_total+"[i] + "+RHS_output+"[i]*dt")
                            if Butcher[s+1][s+1] != 0:
                                RK_lhs_list.append(RHS_output)
                                if Butcher[s+1][s+1] != 1:
                                    RK_rhs_list.append(y_n+"[i] + "+RHS_output+"[i]*dt*("+sp.ccode(Butcher[s+1][s+1]).replace("L","")+")")
                                else:
                                    RK_rhs_list.append(y_n+"[i] + "+RHS_output+"[i]*dt")
                        post_RHS_output = RHS_output
                    if s == num_steps-1:  # If on the final step
                        if Butcher[num_steps][s+1] != 0:
                            RK_lhs_list.append(y_n)
                            if Butcher[num_steps][s+1] != 1:
                                RK_rhs_list.append(y_n+"[i] + "+y_nplus1_running_total+"[i] + "+RHS_output+"[i]*dt*("+sp.ccode(Butcher[num_steps][s+1]).replace("L","")+")")
                            else:
                                RK_rhs_list.append(y_n+"[i] + "+y_nplus1_running_total+"[i] + "+RHS_output+"[i]*dt)")
                        post_RHS_output = y_n
                    body += single_RK_substep(
                        commentblock=indent + "// -={ START k" + str(s + 1) + " substep }=-",
                        RHS_str=RHS_string,
                        RHS_input_str=RHS_input, RHS_output_str=RHS_output,
                        RK_lhss_list=RK_lhs_list, RK_rhss_list=RK_rhs_list,
                        post_RHS_list=[post_RHS_string],
                        post_RHS_output_list=[post_RHS_output]) + "// -={ END k" + str(s + 1) + " substep }=-\n\n"

    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=indent_Ccode(body, "  "),
        rel_path_to_Cparams=os.path.join("."))


# Set NRPy_basic_defines
def NRPy_basic_defines_MoL_timestepping_struct(MoL_method="RK4"):
    y_n_gridfunctions, non_y_n_gridfunctions_list, diagnostic_gridfunctions_point_to = \
        generate_gridfunction_names(MoL_method=MoL_method)
    # Step 3.b: Create MoL_timestepping struct:
    indent = "  "
    Nbd = "typedef struct __MoL_gridfunctions_struct__ {\n"
    Nbd += indent + "REAL *restrict " + y_n_gridfunctions + ";\n"
    for gfs in non_y_n_gridfunctions_list:
        Nbd += indent + "REAL *restrict " + gfs + ";\n"
    Nbd += indent + "REAL *restrict diagnostic_output_gfs;\n"
    Nbd += "} MoL_gridfunctions_struct;\n"

    return Nbd


# Finally declare the master registration function
def MoL_register_C_functions_and_NRPy_basic_defines(MoL_method = "RK4",
            RHS_string =  "rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, RK_INPUT_GFS, RK_OUTPUT_GFS);",
       post_RHS_string = "apply_bcs(Nxx,Nxx_plus_2NGHOSTS, RK_OUTPUT_GFS);"):
    NRPy_basic_defines_MoL_timestepping_struct(MoL_method = MoL_method)
    for which_gfs in ["y_n_gfs", "non_y_n_gfs"]:
        add_to_Cfunction_dict_MoL_malloc(MoL_method, which_gfs)
        add_to_Cfunction_dict_MoL_free_memory(MoL_method, which_gfs)
    add_to_Cfunction_dict_MoL_step_forward_one_timestep(MoL_method, RHS_string, post_RHS_string)
