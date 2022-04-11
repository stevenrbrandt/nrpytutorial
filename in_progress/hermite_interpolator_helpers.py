# hermite_interpolator_helpers.py:
#  This module provides supporting functions for
#    hermite_interpolator.py, which is documented in
#    the NRPy+ tutorial notebook:
#   Tutorial-Hermite_Interpolator.ipynb ,
#
# Depends primarily on: outputC.py and grid.py.

# Author: Maria C. Babiuc Hamilton template courtesy Zachariah B. Etienne
#         babiuc **at** marshall **dot* edu
from outputC import superfast_uniq, outputC, outC_function_dict, add_to_Cfunction_dict  # NRPy+: Core C code output module
import NRPy_param_funcs as par      # NRPy+: parameter interface
import sympy as sp                  # SymPy: The Python computer algebra package upon which NRPy+ depends
import grid as gri                  # NRPy+: Functions having to do with numerical grids
import sys                          # Standard Python module for multiplatform OS-level functions
from collections import namedtuple  # Standard Python: Enable namedtuple data type

HIparams = namedtuple('HIparams', 'PRECISION HI_DM_order enable_HI_functions enable_SIMD DIM MemAllocStyle upwindcontrolvec fullindent outCparams')

################################################
# STEP 1: COMPUTE FROM LIST OF SYMPY EXPRESSIONS

def generate_list_of_interp_vars_from_lhrh_sympyexpr_list(sympyexpr_list,HIparams):
    """
    Generate from list of SymPy expressions in the form
    [lhrh(lhs=var, rhs=expr),lhrh(...),...]
    all interpolator expressions.
    :param sympyexpr_list <- list of SymPy expressions in the form [lhrh(lhs=var, rhs=expr),lhrh(...),...]:
    :return list of interpolator variables; creating _ddnD in case upwinding is enabled with control vector:
    >>> from outputC import lhrh
    >>> import indexedexp as ixp
    >>> import grid as gri
    >>> import NRPy_param_funcs as par
    >>> from hermite_interpolator_helpers import generate_list_of_interp_vars_from_lhrh_sympyexpr_list,HIparams
    >>> aDD     = ixp.register_gridfunctions_for_single_rank2("EVOL","aDD","sym01")
    >>> aDD_dDD = ixp.declarerank4("aDD_dDD","sym01_sym23")
    >>> aDD_dupD = ixp.declarerank3("aDD_dupD","sym01")
    >>> betaU   = ixp.register_gridfunctions_for_single_rank1("EVOL","betaU")
    >>> a0,a1,b,c = par.Cparameters("REAL",__name__,["a0","a1","b","c"],1)
    >>> HIparams.upwindcontrolvec=betaU
    >>> exprlist = [lhrh(lhs=a0,rhs=b*aDD[1][0] + b*aDD_dDD[2][1][2][1] + c*aDD_dDD[0][1][1][0]), \
                    lhrh(lhs=a1,rhs=aDD_dDD[1][0][0][1] + c*aDD_dupD[0][2][1]*betaU[1])]
    >>> generate_list_of_interp_vars_from_lhrh_sympyexpr_list(exprlist,HIparams)
    [aDD_dDD0101, aDD_dDD1212, aDD_ddnD021, aDD_dupD021]
    """
    # Step 1a:
    # Create a list of free symbols in the sympy expr list
    #     that are registered neither as gridfunctions nor
    #     as C parameters. These *must* be interpolators,
    #     so we call the list "list_of_interp_vars"
    list_of_interp_vars_with_duplicates = []
    for expr in sympyexpr_list:
        for var in expr.rhs.free_symbols:
            vartype = gri.variable_type(var)
            if vartype == "other":
                # vartype=="other" should ONLY refer to interpolators, so
                #    if "_dD" or variants do not appear in a variable classified
                #    neither as a gridfunction nor a Cparameter, then error out.
                if ("_dD"   in str(var)) or \
                   ("_dKOD" in str(var)) or \
                   ("_dupD" in str(var)) or \
                   ("_ddnD" in str(var)):
                    list_of_interp_vars_with_duplicates.append(var)
                else:
                    print("Error: Unregistered variable \""+str(var)+"\" in SymPy expression for "+expr.lhs)
                    print("All variables in SymPy expressions passed to HI_outputC() must be registered")
                    print("in NRPy+ as either a gridfunction or Cparameter, by calling")
                    print(str(var)+" = register_gridfunctions...() (in ixp/grid) if \""+str(var)+"\" is a gridfunction, or")
                    print(str(var)+" = Cparameters() (in par) otherwise (e.g., if it is a free parameter set at C runtime).")
                    sys.exit(1)
    list_of_interp_vars = superfast_uniq(list_of_interp_vars_with_duplicates)

    # Step 1b: For each variable with suffix _dupD, append to
    #          the list_of_interp_vars the corresponding _ddnD.
    if HIparams.upwindcontrolvec != "":
        for var in list_of_interp_vars:
            if "_dupD" in str(var):
                list_of_interp_vars.append(sp.sympify(str(var).replace("_dupD", "_ddnD")))

    # Finally, sort the list_of_interp_vars. This ensures
    #     consistency in the C code output, and might even be
    #     tuned to reduce cache misses.
    #     Thanks to Aaron Meurer for this nice one-liner!
    return sorted(list_of_interp_vars, key=sp.default_sort_key)
############################################################

############################################################
# STEP 2: EXTRACT INFORMATION FROM LIST OF SYMPY EXPRESSIONS

def extract_from_list_of_interp_vars__base_gfs_and_interp_ops_lists(list_of_interp_vars):
    """ Extract from list_of_interp_vars a list of base gridfunctions
        and a list of interpolator operators.
    :param list_of_interp_vars:
    :return list_of_base_gridfunctions, list_of_interp_operators:
    >>> from hermite_interpolator_helpers import extract_from_list_of_interp_vars__base_gfs_and_interp_ops_lists
    >>> extract_from_list_of_interp_vars__base_gfs_and_interp_ops_lists(["aDD_dD012","hDD_dDD0112"])
    (['aDD01', 'aDD01', 'vetU2', 'hDD01'], ['dD2', 'dDD12'])
    """
    list_of_base_gridfunction_names_in_interps = []
    list_of_interp_operators = []
    # Step 2a:
    # For each var in "list_of_interp_vars", determine the
    #     base gridfunction name and interpolator operator.
    for var in list_of_interp_vars:
        # Step 2a.1: Check that the number of integers appearing
        #            in the suffix of a variable name matches the
        #            number of U's + D's in the variable name:
        varstr = str(var)
        num_UDs = 0
        for i in range(len(varstr)):
            if varstr[i] == 'D' or varstr[i] == 'U':
                num_UDs += 1
        num_digits = 0
        i = len(varstr) - 1
        while varstr[i].isdigit():
            num_digits += 1
            i -= 1
        if num_UDs != num_digits:
            print("Error: " + varstr + " has " + str(num_UDs) + " U's and D's, but ")
            print(str(num_digits) + " integers at the end. These must be equal.")
            print("Please rename your gridfunction.")
            sys.exit(1)
        # Step 2a.2: Based on the variable name, find the rank of
        #            the underlying gridfunction of which we're
        #            trying to apply the interpolator.
        rank = 0  # rank = "number of juxtaposed U's and D's before the underscore in a interpolator expression"
        underscore_position = -1
        for i in range(len(varstr) - 1, -1, -1):
            if underscore_position > 0 and (varstr[i] == "U" or varstr[i] == "D"):
                rank += 1
            if varstr[i] == "_":
                underscore_position = i

        # Step 2a.3: Based on the variable name, find the order
        #            of the interpolator we're trying to apply.
        interp_order = 0  # interp_order = "number of Dimensions in a interpolator expression"
        for i in range(underscore_position + 1, len(varstr)):
            if (varstr[i] == "D"):
                interp_order += 1

        # Step 2a.4: Based on interpolator order and rank,
        #            store the base gridfunction name in
        #            list_of_base_gridfunction_names_in_interps[]
        list_of_base_gridfunction_names_in_interps.append(varstr[0:underscore_position] +
                                             varstr[len(varstr) - interp_order - rank:len(varstr) - interp_order])
        list_of_interp_operators.append(varstr[underscore_position + 1:len(varstr) - interp_order - rank] +
                               varstr[len(varstr) - interp_order:len(varstr)])
    return list_of_base_gridfunction_names_in_interps, list_of_interp_operators
###################################################################

###################################################################
# STEP 4: GENERATE C CODE FOR READING NEEDED GRIDPOINTS FROM MEMORY

from operator import itemgetter

def type__var(in_var,HIparams, AddPrefix_for_UpDownWindVars=True):
    """ Outputs [type] [variable name]; e.g.,
    "const double variable"

    :param in_var: Variable name
    :param HIparams: Parameters used in the hermite_interpolator codegen
    :param AddPrefix_for_UpDownWindVars: Boolean -- add a prefix to up/downwind variables?
    :return: Return [type] [variable name]
    >>> from hermite_interpolator_helpers import type__var, HIparams
    >>> HIparams.enable_SIMD = "True"
    >>> type__var("aDD00",HIparams)
    \'const REAL_SIMD_ARRAY aDD00\'

    >>> from hermite_interpolator_helpers import type__var, HIparams
    >>> HIparams.enable_SIMD = "False"
    >>> HIparams.PRECISION = "double"
    >>> type__var("variable",HIparams)
    \'const double variable\'
    """
    varname = str(in_var)
    # Disable prefixing upwinded and downwinded variables
    # if the upwind control vector algorithm is disabled.
    if HIparams.upwindcontrolvec == "":
        AddPrefix_for_UpDownWindVars = False
    if AddPrefix_for_UpDownWindVars:
        if "_dupD" in varname:  # Variables suffixed with "_dupD" are set
            #                    to be the "pure" upwinded interpolator,
            #                    before the upwinding algorithm has been
            #                    applied. However, when they are used
            #                    in the RHS expressions, it is assumed
            #                    that the up. algorithm has been applied.
            #                    To ensure consistency we rename all
            #                    _dupD suffixed variables as
            #                    _dupDPUREUPWIND, and use them as input
            #                    into the upwinding algorithm. The output
            #                    will be the original _dupD variable.
            varname = "UpwindAlgInput"+varname
        if "_ddnD" in varname:  # For consistency with _dupD
            varname = "UpwindAlgInput"+varname
    if HIparams.enable_SIMD == "True":
        return "const REAL_SIMD_ARRAY " + varname
    return "const " + HIparams.PRECISION + " " + varname

def read_from_memory_Ccode_onept(gfname,idx, HIparams):
    """

    :param gfname: gridfunction name; a string
    :param idx: Grid index relative to (i0,i1,i2,i3); e.g., "0,1,2,3" -> (i0,i1+1,i2+2,i3+3); later indices ignored for DIM<4
    :param HIparams: Parameters used in the hermite_interpolator codegen
    :return: C code string for reading in this gridfunction at point idx from memory
    >>> import indexedexp as ixp
    >>> from hermite_interpolator_helpers import HIparams, read_from_memory_Ccode_onept
    >>> HIparams.DIM = 3
    >>> HIparams.upwindcontrolvec = ""
    >>> HIparams.enable_SIMD = "True"
    >>> vetU = ixp.register_gridfunctions_for_single_rank1("EVOL","vetU",HIparams.DIM)
    >>> read_from_memory_Ccode_onept("vetU0","0,1,-2,300",HIparams)
    \'const REAL_SIMD_ARRAY vetU0_i0_i1p1_i2m2 = ReadSIMD(&in_gfs[IDX4S(VETU0GF, i0,i1+1,i2-2)]);\\n\'
    """
    idxsplit = idx.split(',')
    idx4 = [int(idxsplit[0]),int(idxsplit[1]),int(idxsplit[2]),int(idxsplit[3])]
    gf_array_name = "in_gfs" # Default array name.
    gfaccess_str = gri.gfaccess(gf_array_name,gfname,ijkl_string(idx4, HIparams))
    if HIparams.enable_SIMD == "True":
        retstring = type__var(gfname, HIparams) + varsuffix(idx4, HIparams) + " = ReadSIMD(&" + gfaccess_str + ");"
    else:
        retstring = type__var(gfname, HIparams) + varsuffix(idx4, HIparams) + " = " + gfaccess_str + ";"
    return retstring+"\n"

def ijkl_string(idx4, HIparams):
    """Generate string for reading gridfunction from specific location in memory
    if DIM==4:
        input: [i,j,k,l]
        output: "i0+i,i1+j,i2+k,i3+l"
    if DIM==3:
        input: [i,j,k,l]
        output: "i0+i,i1+j,i2+k"
    etc.
    :param idx4: An array of 4 integers, indicating a Hermite interpolator grid index relative to where the HI is being computed
    :param HIparams: Parameters used in the hermite_interpolator codegen
    :return: DIM==3 input [i,j,k,l] -> output "i0+i,i1+j,i2+k"
    >>> from hermite_interpolator_helpers import ijkl_string, HIparams
    >>> HIparams.DIM = 4
    >>> ijkl_string([-2,1,0,-1], HIparams)
    \'i0-2,i1+1,i2,i3-1\'

    >>> from hermite_interpolator_helpers import ijkl_string, HIparams
    >>> HIparams.DIM = 3
    >>> ijkl_string([-2,-1,-1,-300], HIparams)
    \'i0-2,i1-1,i2-1\'
    """
    retstring = ""
    for i in range(HIparams.DIM):
        if i > 0:
            # Add a comma
            retstring += ","
        retstring += "i" + str(i) + "+" + str(idx4[i])
    return retstring.replace("+-", "-").replace("+0", "")

def varsuffix(idx4, HIparams):
    """Generate string for suffixing single point read in from memory
    Example: If a gridfunction is named hDD00, and we want to read from memory data at i0+1,i1,i2-1,
    we store the value of this gridfunction as hDD00_i0p1_i1_i2m1; this function provides the suffix.
    if DIM==3:
    input: [0,2,1,-100]
    output: "_i0_i1p2_i2p1"

    :param idx4: An array of 4 integers, indicating a Hermite interpolator grid index relative to where the HI is being computed
    :param HIparams: Parameters used in the hermite_interpolator codegen
    :return: returns suffix to uniquely name a point of data for a gridfunction
    >>> from hermite_interpolator_helpers import varsuffix, HIparams
    >>> HIparams.DIM=3
    >>> varsuffix([-2,0,-1,-300], HIparams)
    \'_i0m2_i1_i2m1\'
    """
    if idx4 == [0, 0, 0, 0]:
        return ""
    return "_" + ijkl_string(idx4, HIparams).replace(",", "_").replace("+", "p").replace("-", "m")

def read_gfs_from_memory(list_of_base_gridfunction_names_in_interps, histencl, sympyexpr_list, HIparams):
    # with open(list_of_base_gridfunction_names_in_interps[0]+".txt","w") as file:
    #     file.write(str(list_of_base_gridfunction_names_in_interps))
    #     file.write(str(histencl))
    #     file.write(str(sympyexpr_list))
    #     file.write(str(HIparams))
    """

    :param list_of_base_gridfunction_names_in_interps:
    :param histencl:
    :param sympyexpr_list:
    :param HIparams:
    :return:
    >>> from outputC import lhrh
    >>> import indexedexp as ixp
    >>> import NRPy_param_funcs as par
    >>> from hermite_interpolator_helpers import generate_list_of_interp_vars_from_lhrh_sympyexpr_list,HIparams
    >>> from hermite_interpolator_helpers import extract_from_list_of_interp_vars__base_gfs_and_interp_ops_lists
    >>> from hermite_interpolator_helpers import read_gfs_from_memory
    >>> from hermite_interpolator import compute_hicoeffs_histencl
    >>> import grid as gri
    >>> gri.glb_gridfcs_list = []
    >>> hDD      = ixp.register_gridfunctions_for_single_rank2("EVOL","hDD","sym01")
    >>> hDD_dD   = ixp.declarerank3("hDD_dD","sym01")
    >>> hDD_dupD = ixp.declarerank3("hDD_dupD","sym01")
    >>> vU       = ixp.register_gridfunctions_for_single_rank1("EVOL","vU")
    >>> a0,a1,b,c = par.Cparameters("REAL",__name__,["a0","a1","b","c"],1)
    >>> par.set_parval_from_str("hermite_interpolator::HI_DIMENSIONS_ORDER",2)
    >>> HIparams.DIM=3
    >>> HIparams.enable_SIMD="False"
    >>> HIparams.PRECISION="double"
    >>> HIparams.MemAllocStyle="012"
    >>> HIparams.upwindcontrolvec=vU
    >>> exprlist = [lhrh(lhs=a0,rhs=b*hDD[1][0] + c*hDD_dD[0][1][1]), \
                    lhrh(lhs=a1,rhs=c*hDD_dupD[0][2][1]*vU[1])]
    >>> list_of_interp_vars = generate_list_of_interp_vars_from_lhrh_sympyexpr_list(exprlist,HIparams)
    >>> list_of_base_gridfunction_names_in_interps, list_of_interp_operators = extract_from_list_of_interp_vars__base_gfs_and_interp_ops_lists(list_of_interp_vars)
    >>> hicoeffs = [[] for i in range(len(list_of_interp_operators))]
    >>> histencl = [[[] for i in range(4)] for j in range(len(list_of_interp_operators))]
    >>> for i in range(len(list_of_interp_operators)): hicoeffs[i], histencl[i] = compute_hicoeffs_histencl(list_of_interp_operators[i])
    >>> print(read_gfs_from_memory(list_of_base_gridfunction_names_in_interps, histencl, exprlist, HIparams))
    const double hDD01_i0_i1m1_i2 = in_gfs[IDX4S(HDD01GF, i0,i1-1,i2)];
    const double hDD01 = in_gfs[IDX4S(HDD01GF, i0,i1,i2)];
    const double hDD01_i0_i1p1_i2 = in_gfs[IDX4S(HDD01GF, i0,i1+1,i2)];
    const double hDD02_i0_i1m2_i2 = in_gfs[IDX4S(HDD02GF, i0,i1-2,i2)];
    const double hDD02_i0_i1m1_i2 = in_gfs[IDX4S(HDD02GF, i0,i1-1,i2)];
    const double hDD02 = in_gfs[IDX4S(HDD02GF, i0,i1,i2)];
    const double hDD02_i0_i1p1_i2 = in_gfs[IDX4S(HDD02GF, i0,i1+1,i2)];
    const double hDD02_i0_i1p2_i2 = in_gfs[IDX4S(HDD02GF, i0,i1+2,i2)];
    const double vU1 = in_gfs[IDX4S(VU1GF, i0,i1,i2)];
    <BLANKLINE>
    """

    # Step 4a: Compile list of points to read from memory
    #          for each gridfunction i, based on list
    #          provided in histencil[i][].
    list_of_points_read_from_memory_with_duplicates = [[] for i in range(len(gri.glb_gridfcs_list))]
    for j in range(len(list_of_base_gridfunction_names_in_interps)):
        interpgfname = list_of_base_gridfunction_names_in_interps[j]
        # Next find the corresponding gridfunction index:
        for i in range(len(gri.glb_gridfcs_list)):
            gfname = gri.glb_gridfcs_list[i].name
            # If the gridfunction for the interpolator matches, then
            #    add to the list of points read from memory:
            if interpgfname == gfname:
                for k in range(len(histencl[j])):
                    list_of_points_read_from_memory_with_duplicates[i].append(str(histencl[j][k][0]) + "," +
                                                                              str(histencl[j][k][1]) + "," +
                                                                              str(histencl[j][k][2]) + "," +
                                                                              str(histencl[j][k][3]))

    # Step 4b: "Zeroth interpolator" case:
    #     If gridfunction appears in expression not
    #     as interpolator (i.e., by itself), it must
    #     be read from memory as well.
    for expr in range(len(sympyexpr_list)):
        for var in sympyexpr_list[expr].rhs.free_symbols:
            vartype = gri.variable_type(var)
            if vartype == "gridfunction":
                for i in range(len(gri.glb_gridfcs_list)):
                    gfname = gri.glb_gridfcs_list[i].name
                    if gfname == str(var):
                        list_of_points_read_from_memory_with_duplicates[i].append("0,0,0,0")


    # Step 4c: Remove duplicates when reading from memory;
    #     do not needlessly read the same variable
    #     from memory twice.
    list_of_points_read_from_memory = [[] for i in range(len(gri.glb_gridfcs_list))]
    for i in range(len(gri.glb_gridfcs_list)):
        list_of_points_read_from_memory[i] = superfast_uniq(list_of_points_read_from_memory_with_duplicates[i])

    # Step 4d: Minimize cache misses:
    #      Sort the list of points read from
    #      main memory by how they are stored
    #      in memory.

    # Step 4d.i: Define a function that maps a gridpoint
    #     index (i,j,k,l) to a unique memory "address",
    #     which will correspond to the correct ordering
    #     of actual memory addresses.
    #
    #     Input: a list of 4 indices, e.g., (i,j,k,l)
    #            corresponding to a gridpoint's *spatial*
    #            index in memory (thus we support up to
    #            4D in space). If spatial dimension is
    #            less than 4D, then just set latter
    #            index/indices to zero. E.g., for 2D
    #            spatial indexing, set (i,j,0,0).
    #     Output: a single number, which when sorted
    #            will yield a unique "address" in memory
    #            such that consecutive addresses are
    #            consecutive in memory.
    def unique_idx(idx4,HIparams):
        # os and sz are set *just for the purposes of ensuring indices are ordered in memory*
        #    Do not modify the values of os and sz.
        os = 50  # offset
        sz = 100 # assumed size in each direction
        if HIparams.MemAllocStyle == "210":
            return str(int(idx4[0])+os + sz*( (int(idx4[1])+os) + sz*( (int(idx4[2])+os) + sz*( int(idx4[3])+os ) ) ))
        if HIparams.MemAllocStyle == "012":
            return str(int(idx4[3])+os + sz*( (int(idx4[2])+os) + sz*( (int(idx4[1])+os) + sz*( int(idx4[0])+os ) ) ))
        print("Error: MemAllocStyle = "+HIparams.MemAllocStyle+" unsupported.")
        sys.exit(1)

    # Step 4d.ii: For each gridfunction and
    #      point read from memory, call unique_idx,
    #      then sort according to memory "address"
    # Input: list_of_points_read_from_memory[gridfunction][point],
    #        gri.glb_gridfcs_list[gridfunction]
    # Output: 1) A list of points to be read from
    #            memory, sorted according to memory
    #            "address":
    #            sorted_list_of_points_read_from_memory[gridfunction][point]
    #        2) A list containing the gridfunction
    #           read at each point, with the number
    #           of elements corresponding exactly
    #           to the total number of points read
    #           from memory for all gridfunctions:
    #           read_from_memory_gf[]
    read_from_memory_gf    = []
    sorted_list_of_points_read_from_memory = [[] for i in range(len(gri.glb_gridfcs_list))]
    for gfidx in range(len(gri.glb_gridfcs_list)):
        # Continue only if reading at least one point of gfidx from memory.
        #     The sorting algorithm at the end of this code block is not
        #     well-defined (will throw an error) if no points of gfidx are
        #     read from memory.
        if len(list_of_points_read_from_memory[gfidx]) > 0:
            read_from_memory_index = []
            for idx in list_of_points_read_from_memory[gfidx]:
                read_from_memory_gf.append(gri.glb_gridfcs_list[gfidx])
                idxsplit = idx.split(',')
                idx4 = [int(idxsplit[0]),int(idxsplit[1]),int(idxsplit[2]),int(idxsplit[3])]
                read_from_memory_index.append(unique_idx(idx4, HIparams))
                # https://stackoverflow.com/questions/13668393/python-sorting-two-lists
                _unused_list, sorted_list_of_points_read_from_memory[gfidx] = \
                    [list(x) for x in zip(*sorted(zip(read_from_memory_index, list_of_points_read_from_memory[gfidx]),
                                                  key=itemgetter(0)))]
    # Step 4e: Create the full C code string
    #      for reading from memory:

    read_from_memory_Ccode = ""
    count = 0
    for gfidx in range(len(gri.glb_gridfcs_list)):
        for pt in range(len(sorted_list_of_points_read_from_memory[gfidx])):
            read_from_memory_Ccode += read_from_memory_Ccode_onept(read_from_memory_gf[count].name,
                                                                   sorted_list_of_points_read_from_memory[gfidx][pt],
                                                                   HIparams)
            count += 1
    return read_from_memory_Ccode
#################################

#################################
# STEP 5: C CODE OUTPUT ROUTINES
def construct_HI_exprs_as_SymPy_exprs(list_of_interp_vars,
                                      list_of_base_gridfunction_names_in_interps, list_of_interp_operators,
                                      hicoeffs, histencl):
    HIexprs = []
    HIlhsvarnames = []
    # Step 5.a.ii.A: Output Hermite interpolator expressions to
    #                Coutput string
    for i in range(len(list_of_interp_vars)):
        HIexprs.append(sp.sympify(0))  # Append a new element to the list of interpolator expressions.
        HIlhsvarnames.append(type__var(list_of_interp_vars[i], HIparams))
        var = list_of_base_gridfunction_names_in_interps[i]
        for j in range(len(hicoeffs[i])):
            varname = str(var) + varsuffix(histencl[i][j], HIparams)
            HIexprs[i] += hicoeffs[i][j] * sp.sympify(varname)

        # Multiply each expression by the appropriate power
        #   of 1/dx[i]
        invdx = []
        for d in range(HIparams.DIM):
            invdx.append(sp.sympify("invdx" + str(d)))
        # First-order or Kreiss-Oliger interpolators:
        if (len(list_of_interp_operators[i]) == 5 and "dKOD" in list_of_interp_operators[i]) or \
                (len(list_of_interp_operators[i]) == 3 and "dD" in list_of_interp_operators[i]) or \
                (len(list_of_interp_operators[i]) == 5 and (
                        "dupD" in list_of_interp_operators[i] or "ddnD" in list_of_interp_operators[i])):
            dirn = int(list_of_interp_operators[i][len(list_of_interp_operators[i]) - 1])
            HIexprs[i] *= invdx[dirn]
        # Second-order interps:
        elif len(list_of_interp_operators[i]) == 5 and "dDD" in list_of_interp_operators[i]:
            dirn1 = int(list_of_interp_operators[i][len(list_of_interp_operators[i]) - 2])
            dirn2 = int(list_of_interp_operators[i][len(list_of_interp_operators[i]) - 1])
            HIexprs[i] *= invdx[dirn1] * invdx[dirn2]
        else:
            print("Error: was unable to parse interpolator operator: ", list_of_interp_operators[i])
            sys.exit(1)
    return HIexprs, HIlhsvarnames

def find_which_op_idx(op, list_of_interp_operators):
    for j in range(len(list_of_interp_operators)):
        if op == list_of_interp_operators[j]:
            return j
    print("Error: could not find operator "+str(op)+" in ",list_of_interp_operators)
    sys.exit(1)


def add_HI_func_to_outC_function_dict(list_of_interp_vars,
                                      list_of_base_gridfunction_names_in_interps, list_of_interp_operators,
                                      hicoeffs, histencl):
    # Step 5.a.ii.A: First construct a list of all the unique Hermite interpolator functions
    list_of_uniq_interp_operators = superfast_uniq(list_of_interp_operators)
    c_type = "REAL"
    if par.parval_from_str("grid::GridFuncMemAccess") == "ETK":
        c_type = "CCTK_REAL"
    func_prefix = "order_"+str(HIparams.HI_DM_order)+"_"
    if HIparams.enable_SIMD == "True":
        c_type = "REAL_SIMD_ARRAY"
        func_prefix = "SIMD_"+func_prefix

    # Stores the needed calls to the functions we're adding to outC_function_dict:
    HIfunccall_list = []
    for op in list_of_uniq_interp_operators:
        which_op_idx = find_which_op_idx(op, list_of_interp_operators)

        rhs_expr = sp.sympify(0)
        for j in range(len(hicoeffs[which_op_idx])):
            var = sp.sympify("f" + varsuffix(histencl[which_op_idx][j], HIparams))
            rhs_expr += hicoeffs[which_op_idx][j] * var

        # Multiply each expression by the appropriate power
        #   of 1/dx[i]
        invdx = []
        used_invdx = [False, False, False, False]
        for d in range(HIparams.DIM):
            invdx.append(sp.sympify("invdx" + str(d)))
        # First-order or Kreiss-Oliger interpolators:
        if ( (len(op) == 5 and "dKOD" in op) or
             (len(op) == 3 and   "dD" in op) or
             (len(op) == 5 and ("dupD" in op or "ddnD" in op)) ):
            dirn = int(op[len(op) - 1])
            rhs_expr *= invdx[dirn]
            used_invdx[dirn] = True
        # Second-order interps:
        elif len(op) == 5 and "dDD" in op:
            dirn1 = int(op[len(op) - 2])
            dirn2 = int(op[len(op) - 1])
            used_invdx[dirn1] = used_invdx[dirn2] = True
            rhs_expr *= invdx[dirn1]*invdx[dirn2]
        else:
            print("Error: was unable to parse interpolator operator: ", op)
            sys.exit(1)

        outfunc_params = ""
        for d in range(HIparams.DIM):
            if used_invdx[d]:
                outfunc_params += "const " + c_type + " invdx" + str(d) + ","

        for j in range(len(hicoeffs[which_op_idx])):
            var = sp.sympify("f" + varsuffix(histencl[which_op_idx][j], HIparams))
            outfunc_params += "const " + c_type + " " + str(var)
            if j != len(hicoeffs[which_op_idx])-1:
                outfunc_params += ","

        for i in range(len(list_of_interp_operators)):
            # print("comparing ",list_of_interp_operators[i],op)
            if list_of_interp_operators[i] == op:
                funccall = type__var(list_of_interp_vars[i], HIparams) + " = " + func_prefix + "f_" + str(op) + "("
                for d in range(HIparams.DIM):
                    if used_invdx[d]:
                        funccall += "invdx" + str(d) + ","
                gfname = list_of_base_gridfunction_names_in_interps[i]
                for j in range(len(hicoeffs[which_op_idx])):
                    funccall += gfname + varsuffix(histencl[which_op_idx][j], HIparams)
                    if j != len(hicoeffs[which_op_idx])-1:
                        funccall += ","
                funccall += ");"
                HIfunccall_list.append(funccall)

        # If the function already exists in the outC_function_dict, then do not add it; move to the next op.
        if func_prefix + "f_" + str(op) not in outC_function_dict:
            p = "preindent=1,enable_SIMD="+HIparams.enable_SIMD+",outCverbose=False,CSE_preprocess=True,includebraces=False"
            outHIstr = outputC(rhs_expr, "retval", "returnstring", params=p)
            outHIstr = outHIstr.replace("retval = ", "return ")
            add_to_Cfunction_dict(desc=" * (__HI_OPERATOR_FUNC__) Hermite interpolator operator for "+str(op).replace("dDD", "second interpolator: ").
                                  replace("dD", "first interpolator: ").replace("dKOD", "Kreiss-Oliger interpolator: ").
                                  replace("dupD", "upwinded interpolator: ").replace("ddnD", "downwinded interpolator: ") + " direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.",
                                  c_type="static " + c_type + " _NOINLINE _UNUSED",
                                  name=func_prefix+"f_" + str(op), enableCparameters=False,
                                  params=outfunc_params, preloop="", body=outHIstr)
    return HIfunccall_list

def construct_Ccode(sympyexpr_list, list_of_interp_vars,
                    list_of_base_gridfunction_names_in_interps,list_of_interp_operators,
                    hicoeffs, histencl, read_from_memory_Ccode, HIparams, Coutput):
    """
    C code is constructed in *up to* 3 parts:
     5.a) Read gridfunctions from memory at needed pts
          for Hermite interpolator; compute hermite_interpolator
          stencils.
     5.b) Implement upwinding algorithm (if relevant)
     5.c) Evaluate SymPy expressions and write to main
          memory
    """
    # Failed Doctest. However, mathematically equivalent with Sympy 1.3
    # :param sympyexpr_list:
    # :param list_of_interp_vars:
    # :param list_of_base_gridfunction_names_in_interps:
    # :param list_of_interp_operators:
    # :param hicoeffs:
    # :param histencl:
    # :param read_from_memory_Ccode:
    # :param HIparams:
    # :param Coutput: The start of the Coutput string; this function's output will be pasted to a copy of Coutput
    # :return: Returns a C code string
    # >>> from outputC import lhrh
    # >>> import indexedexp as ixp
    # >>> import NRPy_param_funcs as par
    # >>> from hermite_interpolator_helpers import generate_list_of_interp_vars_from_lhrh_sympyexpr_list,HIparams
    # >>> from hermite_interpolator_helpers import extract_from_list_of_interp_vars__base_gfs_and_interp_ops_lists
    # >>> from hermite_interpolator_helpers import read_gfs_from_memory, construct_Ccode
    # >>> from hermite_interpolator import compute_hicoeffs_histencl
    # >>> import grid as gri
    # >>> gri.glb_gridfcs_list = []
    # >>> hDD      = ixp.register_gridfunctions_for_single_rank2("EVOL","hDD","sym01")
    # >>> hDD_dD   = ixp.declarerank3("hDD_dD","sym01")
    # >>> hDD_dupD = ixp.declarerank3("hDD_dupD","sym01")
    # >>> vU       = ixp.register_gridfunctions_for_single_rank1("EVOL","vU")
    # >>> a0,a1,b,c = par.Cparameters("REAL",__name__,["a0","a1","b","c"],1)
    # >>> par.set_parval_from_str("hermite_interpolator::HI_DIMENSIONS_ORDER",2)
    # >>> HIparams.DIM=3
    # >>> HIparams.enable_SIMD="False"
    # >>> HIparams.enable_HI_functions=False
    # >>> HIparams.PRECISION="double"
    # >>> HIparams.MemAllocStyle="012"
    # >>> HIparams.upwindcontrolvec=vU
    # >>> HIparams.fullindent=""
    # >>> HIparams.outCparams="outCverbose=False"
    # >>> exprlist = [lhrh(lhs=a0,rhs=b*hDD[1][0] + c*hDD_dD[0][1][1]), \
    #                 lhrh(lhs=a1,rhs=c*hDD_dupD[0][2][1]*vU[1])]
    # >>> list_of_interp_vars = generate_list_of_interp_vars_from_lhrh_sympyexpr_list(exprlist,HIparams)
    # >>> list_of_base_gridfunction_names_in_interps, list_of_interp_operators = extract_from_list_of_interp_vars__base_gfs_and_interp_ops_lists(list_of_interp_vars)
    # >>> hicoeffs = [[] for i in range(len(list_of_interp_operators))]
    # >>> histencl = [[[] for i in range(4)] for j in range(len(list_of_interp_operators))]
    # >>> for i in range(len(list_of_interp_operators)): hicoeffs[i], histencl[i] = compute_hicoeffs_histencl(list_of_interp_operators[i])
    # >>> memread_Ccode = read_gfs_from_memory(list_of_base_gridfunction_names_in_interps, histencl, exprlist, HIparams)
    # >>> print(construct_Ccode(exprlist, list_of_interp_vars, \
    #           list_of_base_gridfunction_names_in_interps, list_of_interp_operators, \
    #           hicoeffs, histencl, memread_Ccode, HIparams, ""))
    # /*
    #  * NRPy+ Hermite Difference Code Generation, Step 1 of 3: Read from main memory and compute Hermite interpolator stencils:
    #  */
    # const double hDD01_i0_i1m1_i2 = in_gfs[IDX4S(HDD01GF, i0,i1-1,i2)];
    # const double hDD01 = in_gfs[IDX4S(HDD01GF, i0,i1,i2)];
    # const double hDD01_i0_i1p1_i2 = in_gfs[IDX4S(HDD01GF, i0,i1+1,i2)];
    # const double hDD02_i0_i1m2_i2 = in_gfs[IDX4S(HDD02GF, i0,i1-2,i2)];
    # const double hDD02_i0_i1m1_i2 = in_gfs[IDX4S(HDD02GF, i0,i1-1,i2)];
    # const double hDD02 = in_gfs[IDX4S(HDD02GF, i0,i1,i2)];
    # const double hDD02_i0_i1p1_i2 = in_gfs[IDX4S(HDD02GF, i0,i1+1,i2)];
    # const double hDD02_i0_i1p2_i2 = in_gfs[IDX4S(HDD02GF, i0,i1+2,i2)];
    # const double vU1 = in_gfs[IDX4S(VU1GF, i0,i1,i2)];
    # const double HIPart1_Rational_1_2 = 1.0/2.0;
    # const double HIPart1_Integer_2 = 2.0;
    # const double HIPart1_Rational_3_2 = 3.0/2.0;
    # const double hDD_dD011 = HIPart1_Rational_1_2*invdx1*(-hDD01_i0_i1m1_i2 + hDD01_i0_i1p1_i2);
    # const double UpwindAlgInputhDD_ddnD021 = invdx1*(-HIPart1_Integer_2*hDD02_i0_i1m1_i2 + HIPart1_Rational_1_2*hDD02_i0_i1m2_i2 + HIPart1_Rational_3_2*hDD02);
    # const double UpwindAlgInputhDD_dupD021 = invdx1*(HIPart1_Integer_2*hDD02_i0_i1p1_i2 - HIPart1_Rational_1_2*hDD02_i0_i1p2_i2 - HIPart1_Rational_3_2*hDD02);
    # const double UpwindControlVectorU1 = vU1;
    # /*
    #  * NRPy+ Hermite Difference Code Generation, Step 2 of 3: Implement upwinding algorithm:
    #  */
    # const double UpWind1 = UPWIND_ALG(UpwindControlVectorU1);
    # const double hDD_dupD021 = UpWind1*(-UpwindAlgInputhDD_ddnD021 + UpwindAlgInputhDD_dupD021) + UpwindAlgInputhDD_ddnD021;
    # /*
    #  * NRPy+ Hermite Difference Code Generation, Step 3 of 3: Evaluate SymPy expressions and write to main memory:
    #  */
    # a0 = b*hDD01 + c*hDD_dD011;
    # a1 = c*hDD_dupD021*vU1;
    # <BLANKLINE>

    def indent_Ccode(Ccode):
        Ccodesplit = Ccode.splitlines()
        outstring = ""
        for i in range(len(Ccodesplit)):
         if Ccodesplit[i] != "":
            if Ccodesplit[i].lstrip().startswith("#"):
                # Remove all indentation from preprocessor statements (lines that start with "#")
                outstring += Ccodesplit[i].lstrip() + '\n'
            else:
                outstring += HIparams.fullindent + Ccodesplit[i] + '\n'
        return outstring.rstrip(" ")  # make sure to remove trailing whitespace!

    # Step 5.a.i: Read gridfunctions from memory at needed pts.
    # *** No need to do anything here; already set in
    #     string "read_from_memory_Ccode". ***

    # Step 5.a.ii: Perform arithmetic needed for Hermite interpolators
    #              associated with input expressions provided in
    #              sympyexpr_list[].rhs.
    #           Note that HIexprs and HIlhsvarnames contain
    #          A) Hermite interpolator expressions (constructed
    #             in steps above) and associated variable names,
    #             and
    #          B) Input expressions sympyexpr_list[], which
    #             in general depend on Hermite interpolator  
    #             variables.
    HIexprs = []
    HIlhsvarnames = []
    if not HIparams.enable_HI_functions:
        HIexprs, HIlhsvarnames = \
            construct_HI_exprs_as_SymPy_exprs(list_of_interp_vars,
                                              list_of_base_gridfunction_names_in_interps, list_of_interp_operators,
                                              hicoeffs, histencl)

    # Compute Hermite interpolators using function calls (instead of inlined calculations)?
    if HIparams.enable_HI_functions:
        # If so, add HI functions to outputC's outC_function_dict (C function dictionary),
        #   AND return the full set of needed calls to these functions (to funccall_list)
        funccall_list = \
            add_HI_func_to_outC_function_dict(list_of_interp_vars,
                                              list_of_base_gridfunction_names_in_interps, list_of_interp_operators,
                                              hicoeffs, histencl)

    # Step 5.b.i: (Upwinded interpolators algorithm, part 1):
    # If an upwinding control vector is specified, determine
    #    which of the elements of the vector will be required.
    #    This ensures that those elements are read from memory.
    # For example, if a symmetry axis is specified,
    #     upwind interpolators with respect to only
    #     two of the three dimensions are used. Here
    #     we find all directions used for upwinding.
    upwind_directions = []
    if HIparams.upwindcontrolvec != "":
        upwind_directions_unsorted_withdups = []
        for interp_op in list_of_interp_operators:
            if "dupD" in interp_op:
                if interp_op[len(interp_op)-1].isdigit():
                    dirn = int(interp_op[len(interp_op)-1])
                    upwind_directions_unsorted_withdups.append(dirn)
                else:
                    print("Error: Derivative operator "+interp_op+" does not contain a direction")
                    sys.exit(1)
        if len(upwind_directions_unsorted_withdups) > 0:
            upwind_directions = superfast_uniq(upwind_directions_unsorted_withdups)
            upwind_directions = sorted(upwind_directions,key=sp.default_sort_key)
        #   If upwind control vector is specified,
        #        add upwind control vectors to the
        #        interpolator expression list, so its
        #        needed elements are read from memory.
        for dirn in upwind_directions:
            HIexprs.append(HIparams.upwindcontrolvec[dirn])
            HIlhsvarnames.append(type__var("UpwindControlVectorU" + str(dirn), HIparams))

    # Step 5.x: Output useful code comment regarding
    #           which step we are on. *At most* this
    #           is a 3-step process:
    #        1. Read from memory & compute HI stencils,
    #        2. Perform upwinding, and
    #        3. Evaluate remaining expressions+write
    #           results to main memory.
    NRPy_HI_StepNumber = 1
    NRPy_HI__Number_of_Steps = 1
    if len(read_from_memory_Ccode) > 0:
        NRPy_HI__Number_of_Steps += 1
    if HIparams.upwindcontrolvec != "" and len(upwind_directions) > 0:
        NRPy_HI__Number_of_Steps += 1

    if len(read_from_memory_Ccode) > 0:
        Coutput += indent_Ccode("/*\n * NRPy+ Hermite Difference Code Generation, Step "
                                + str(NRPy_HI_StepNumber) + " of " + str(NRPy_HI__Number_of_Steps) +
                                ": Read from main memory and compute Hermite interpolator stencils:\n */\n")
        NRPy_HI_StepNumber = NRPy_HI_StepNumber + 1
        if HIparams.enable_HI_functions:
            # Compute Hermite interpolators using function calls (instead of inlined calculations)
            Coutput += indent_Ccode(read_from_memory_Ccode)
            for funccall in funccall_list:
                Coutput += indent_Ccode(funccall)
            if HIparams.upwindcontrolvec != "":
                # Compute Hermite interpolators using inlined calculations
                params = HIparams.outCparams
                # We choose the CSE temporary variable prefix "HIpart1" for the Hermite interpolator coefficients:
                params += ",CSE_varprefix=HIPart1,includebraces=False,CSE_preprocess=True,SIMD_find_more_subs=True"
                Coutput += indent_Ccode(outputC(HIexprs, HIlhsvarnames, "returnstring", params=params))

        else:
            # Compute Hermite interpolators using inlined calculations
            params = HIparams.outCparams.replace("preindent=1", "preindent=0")  # Remove an unnecessary indentation
            # We choose the CSE temporary variable prefix "HIpart1" for the Hermite interpolator coefficients:
            params += ",CSE_varprefix=HIPart1,includebraces=False,CSE_preprocess=True,SIMD_find_more_subs=True"
            Coutput += indent_Ccode(outputC(HIexprs, HIlhsvarnames, "returnstring",params=params,
                                            prestring=read_from_memory_Ccode))

    # Step 5.b.ii: Implement control-vector upwinding algorithm.
    if HIparams.upwindcontrolvec != "":
        if len(upwind_directions) > 0:
            Coutput += indent_Ccode("/*\n * NRPy+ Hermite Difference Code Generation, Step "
                                    + str(NRPy_HI_StepNumber) + " of " + str(NRPy_HI__Number_of_Steps) +
                                    ": Implement upwinding algorithm:\n */\n")
            NRPy_HI_StepNumber = NRPy_HI_StepNumber + 1
            if HIparams.enable_SIMD == "True":
                for n in ["0", "1"]:
                    Coutput += indent_Ccode("const double tmp_upwind_Integer_"+n+" = "+n+".000000000000000000000000000000000;\n")
                    Coutput += indent_Ccode("const REAL_SIMD_ARRAY upwind_Integer_"+n+" = ConstSIMD(tmp_upwind_Integer_"+n+");\n")
            for dirn in upwind_directions:
                Coutput += indent_Ccode(type__var("UpWind" + str(dirn), HIparams) +
                                        " = UPWIND_ALG(UpwindControlVectorU" + str(dirn) + ");\n")
        upwindU = [sp.sympify(0) for i in range(HIparams.DIM)]
        for dirn in upwind_directions:
            upwindU[dirn] = sp.sympify("UpWind" + str(dirn))
        upwind_expr_list, var_list = [], []
        for i in range(len(list_of_interp_vars)):
            if len(list_of_interp_operators[i]) == 5 and ("dupD" in list_of_interp_operators[i]):
                var_dupD = sp.sympify("UpwindAlgInput" + str(list_of_interp_vars[i]))
                var_ddnD = sp.sympify("UpwindAlgInput" + str(list_of_interp_vars[i]).replace("_dupD", "_ddnD"))
                upwind_dirn = int(list_of_interp_operators[i][len(list_of_interp_operators[i]) - 1])
                upwind_expr = upwindU[upwind_dirn] * (var_dupD - var_ddnD) + var_ddnD
                upwind_expr_list.append(upwind_expr)
                var_list.append(type__var(str(list_of_interp_vars[i]), HIparams, AddPrefix_for_UpDownWindVars=False))
        # For convenience, we require type__var() above to
        # prefix up/downwinded variables with "UpwindAlgInput".
        # Here we do not wish to have this prefix.
        Coutput += indent_Ccode(outputC(upwind_expr_list, var_list,
                                        "returnstring", params=HIparams.outCparams + ",CSE_varprefix=HIPart2,includebraces=False"))

    # Step 5.c.i: Add input RHS & LHS expressions from
    #             sympyexpr_list[]
    Coutput += indent_Ccode("/*\n * NRPy+ Hermite Difference Code Generation, Step "
                            + str(NRPy_HI_StepNumber) + " of " + str(NRPy_HI__Number_of_Steps) +
                            ": Evaluate SymPy expressions and write to main memory:\n */\n")
    exprs = []
    lhsvarnames = []
    for i in range(len(sympyexpr_list)):
        exprs.append(sympyexpr_list[i].rhs)
        if HIparams.enable_SIMD == "True":
            lhsvarnames.append("const REAL_SIMD_ARRAY __RHS_exp_" + str(i))
        else:
            lhsvarnames.append(sympyexpr_list[i].lhs)

    # Step 5.c.ii: Write output to gridfunctions specified in
    #              sympyexpr_list[].lhs.
    write_to_mem_string = ""
    if HIparams.enable_SIMD == "True":
        for i in range(len(sympyexpr_list)):
            write_to_mem_string += "WriteSIMD(&" + sympyexpr_list[i].lhs + ", __RHS_exp_" + str(i) + ");\n"

    # outputC requires as its second argument a list of strings.
    #   Sometimes when the lhs's are simple constants, but the inputs
    #   contain gridfunctions, it is necessary to convert the lhs's
    #   to strings:
    lhsvarnamestrings = []
    for lhs in lhsvarnames:
        lhsvarnamestrings.append(str(lhs))

    Coutput += indent_Ccode(outputC(exprs, lhsvarnamestrings, "returnstring",
                                    params=HIparams.outCparams + ",CSE_varprefix=HIPart3,includebraces=False,preindent=0",
                                    prestring="", poststring=write_to_mem_string))

    return Coutput
#################################

if __name__ == "__main__":
    import doctest
    sys.exit(doctest.testmod()[0])
