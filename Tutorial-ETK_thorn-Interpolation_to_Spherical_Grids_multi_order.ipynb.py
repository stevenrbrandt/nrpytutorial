#!/usr/bin/env python






import cmdline_helper as cmd  # NRPy+: Multi-platform Python command-line interface
import shutil, os             # Standard Python modules for multiplatform OS-level functions, benchmarking

Ccodesdir = "interp_sphgrid_MO_ETK"
shutil.rmtree(Ccodesdir, ignore_errors=True)
cmd.mkdir(Ccodesdir)
cmd.mkdir(os.path.join(Ccodesdir,"src/"))





get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/Interpolate_to_sph_grid.h', '\nvoid Interpolate_to_sph_grid(cGH *cctkGH,CCTK_INT interp_num_points, CCTK_INT interp_order,\n                             CCTK_REAL *point_x_temp,CCTK_REAL *point_y_temp,CCTK_REAL *point_z_temp,\n                             const CCTK_STRING input_array_names[1], CCTK_REAL *output_f[1]) {\n  DECLARE_CCTK_PARAMETERS;\n  CCTK_INT ierr;\n\n  const CCTK_INT NUM_INPUT_ARRAYS=1;\n  const CCTK_INT NUM_OUTPUT_ARRAYS=1;\n\n  CCTK_STRING coord_system = "cart3d";\n\n  // Set up handles\n  const CCTK_INT coord_system_handle = CCTK_CoordSystemHandle(coord_system);\n  if (coord_system_handle < 0) {\n    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,\n        "can\'t get coordinate system handle for coordinate system \\"%s\\"!",\n               coord_system);\n  }\n\n  const CCTK_INT operator_handle = CCTK_InterpHandle(interpolator_name);\n  if (operator_handle < 0)\n    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,\n               "couldn\'t find interpolator \\"%s\\"!",\n               interpolator_name);\n\n  char interp_order_string[10];\n  snprintf(interp_order_string, 10, "order=%d", interp_order);\n  CCTK_STRING interpolator_pars = interp_order_string;\n  CCTK_INT param_table_handle = Util_TableCreateFromString(interpolator_pars);\n  if (param_table_handle < 0) {\n    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,\n               "bad interpolator parameter(s) \\"%s\\"!",\n               interpolator_pars);\n  }\n\n  CCTK_INT operand_indices[NUM_INPUT_ARRAYS]; //NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS];\n  for(int i = 0 ; i < NUM_INPUT_ARRAYS  ; i++) {\n    operand_indices[i] = i;\n  }\n  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,\n                        operand_indices, "operand_indices");\n\n\n  CCTK_INT operation_codes[NUM_INPUT_ARRAYS];\n  for(int i = 0 ; i < NUM_INPUT_ARRAYS  ; i++) {\n    operation_codes[i] = 0;\n  }\n  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,\n                        operation_codes, "operation_codes");\n\n  const void* interp_coords[3]\n    = { (const void *) point_x_temp,\n        (const void *) point_y_temp,\n        (const void *) point_z_temp };\n\n  CCTK_INT input_array_indices[NUM_INPUT_ARRAYS];\n  for(int i = 0 ; i < NUM_INPUT_ARRAYS ; i++) {\n    input_array_indices[i] = CCTK_VarIndex(input_array_names[i]);\n    if(input_array_indices[i] < 0) {\n      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,\n        "COULD NOT FIND VARIABLE \'%s\'.",\n        input_array_names[i]);\n      exit(1);\n    }\n  }\n\n  CCTK_INT output_array_types[NUM_OUTPUT_ARRAYS];\n  for(int i = 0 ; i < NUM_OUTPUT_ARRAYS ; i++) {\n    output_array_types[i] = CCTK_VARIABLE_REAL;\n  }\n\n  void * output_arrays[NUM_OUTPUT_ARRAYS]\n    = { (void *) output_f[0] };\n\n  // actual interpolation call\n  ierr = CCTK_InterpGridArrays(cctkGH,\n                               3, // number of dimensions\n                               operator_handle,\n                               param_table_handle,\n                               coord_system_handle,\n                               interp_num_points,\n                               CCTK_VARIABLE_REAL,\n                               interp_coords,\n                               NUM_INPUT_ARRAYS, // Number of input arrays\n                               input_array_indices,\n                               NUM_OUTPUT_ARRAYS, // Number of output arrays\n                               output_array_types,\n                               output_arrays);\n  if (ierr<0) {\n    CCTK_WARN(1,"interpolation screwed up");\n    Util_TableDestroy(param_table_handle);\n    exit(1);\n  }\n\n  ierr = Util_TableDestroy(param_table_handle);\n  if (ierr != 0) {\n    CCTK_WARN(1,"Could not destroy table");\n    exit(1);\n  }\n}\n')





get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/Set_up_interp_points_on_sph_grid.h', '\nvoid sph_grid_Interpolate_many_pts__set_interp_pts(CCTK_ARGUMENTS) {\n  DECLARE_CCTK_ARGUMENTS;\n  DECLARE_CCTK_PARAMETERS;\n\n  CCTK_REAL dx0 = log( (Rout - R0) / (Rin - R0) ) / ((CCTK_REAL)N0);\n  CCTK_REAL dx1 =      1.0 / ((CCTK_REAL)N1);\n  CCTK_REAL dx2 = 2.0*M_PI / ((CCTK_REAL)N2);\n  CCTK_REAL x0_beg = log( Rin - R0 );\n  CCTK_INT which_pt = 0;\n  for(CCTK_INT k=0;k<N2;k++) for(CCTK_INT j=0;j<N1;j++) for(CCTK_INT i=0;i<N0;i++) {\n    CCTK_REAL x0_i = x0_beg + ((CCTK_REAL)i + 0.5)*dx0;\n    CCTK_REAL rr = R0 + exp(x0_i);\n\n    CCTK_REAL x1_j = x1_beg + ((CCTK_REAL)j + 0.5)*dx1;\n    CCTK_REAL th = -1e300;\n    if(theta_option == 1) {\n       th = th_c + (M_PI - 2.0*th_c)*x1_j + xi*sin(2.0*M_PI*x1_j);\n    } else if (theta_option == 2) {\n       th = M_PI/2.0 * ( 1.0 + (1.0 - xi)*(2.0*x1_j - 1.0) + (xi - 2.0*th_c/M_PI)*pow(2.0*x1_j - 1.0 ,th_n) );\n    } else {\n       printf("Error: theta_option = %d NOT SUPPORTED.",theta_option);\n       exit(1);\n    }\n\n    CCTK_REAL x2_k = x2_beg + ((CCTK_REAL)k + 0.5)*dx2;\n    CCTK_REAL ph = x2_k;\n\n    points_x[which_pt] = rr*sin(th)*cos(ph);\n    points_y[which_pt] = rr*sin(th)*sin(ph);\n    points_z[which_pt] = rr*cos(th);\n    which_pt++;\n  }\n}\n')





get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/output_to_file.h', '\n#include "define_NumInterpFunctions.h"\n\n// output_to_file() starts order and InterpCounter both with the value 1\nvoid output_to_file(CCTK_ARGUMENTS,char gf_name[100],int *order,CCTK_REAL *output_f[1]) {\n  DECLARE_CCTK_ARGUMENTS;\n  DECLARE_CCTK_PARAMETERS;\n\n  char filename[100];\n  sprintf (filename, "%s/interp_sph_grids_MO.dat", out_dir);\n  FILE *file;\n  if(*InterpCounter == 1 && *order==1) {\n    file = fopen (filename,"w");\n  } else {\n    file = fopen (filename,"a+");\n  }\n  if (! file) {\n    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,\n                "interp_sph_grid__ET_thorn: Cannot open output file \'%s\'", filename);\n    exit(1);\n  }\n\n  fwrite(gf_name, 100*sizeof(char), 1, file);\n  fwrite(order, sizeof(CCTK_INT), 1, file);\n\n  fwrite(&N0, sizeof(CCTK_INT), 1, file);\n  fwrite(&R0, sizeof(CCTK_REAL), 1, file);\n  fwrite(&Rin, sizeof(CCTK_REAL), 1, file);\n  fwrite(&Rout, sizeof(CCTK_REAL), 1, file);\n\n  fwrite(&N1, sizeof(CCTK_INT), 1, file);\n  fwrite(&x1_beg, sizeof(CCTK_REAL), 1, file);\n  fwrite(&theta_option, sizeof(CCTK_INT), 1, file);\n  fwrite(&th_c, sizeof(CCTK_REAL), 1, file);\n  fwrite(&xi, sizeof(CCTK_REAL), 1, file);\n  fwrite(&th_n, sizeof(CCTK_INT), 1, file);\n\n  fwrite(&N2, sizeof(CCTK_INT), 1, file);\n  fwrite(&x2_beg, sizeof(CCTK_REAL), 1, file);\n\n  CCTK_REAL magic_number = 1.130814081305130e-21;\n  fwrite(&magic_number, sizeof(CCTK_REAL), 1, file);\n  fwrite(&cctk_iteration, sizeof(CCTK_INT), 1, file);\n  fwrite(&cctk_time, sizeof(CCTK_REAL), 1, file);\n  for(CCTK_INT i=0;i<1;i++) {\n    fwrite(output_f[i], sizeof(CCTK_REAL)*N0*N1*N2, 1, file);\n  }\n\n  fclose(file);\n}\n')





get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/main_function.cc', '\n// Include needed ETK & C library header files:\n#include <stdio.h>\n#include <stdlib.h>\n#include <string.h>\n#include <math.h>\n// Needed for dealing with Cactus/ETK infrastructure\n#include "cctk.h"\n#include "cctk_Arguments.h"\n#include "cctk_Parameters.h"\n// Needed for low-level interpolation functions\n#include "util_Table.h"\n#include "util_String.h"\n\n// Include locally-defined C++ functions:\n#include "Set_up_interp_points_on_sph_grid.h"\n#include "Interpolate_to_sph_grid.h"\n#include "output_to_file.h"\n#include "get_gf_name.h"\n\nvoid Interpolate_to_sph_grid_main_function(CCTK_ARGUMENTS) {\n  DECLARE_CCTK_ARGUMENTS;\n  DECLARE_CCTK_PARAMETERS;\n\n  // Perform interpolation only at iteration == interp_out_iteration:\n  if(cctk_iteration != interp_out_iteration) return;\n\n  // Set up spherically sampled interpolation grid arrays points_x,points_y,points_z:\n  sph_grid_Interpolate_many_pts__set_interp_pts(CCTK_PASS_CTOC);\n\n  // Set up output array:\n  CCTK_REAL *output_f[1];\n  output_f[0] = output_interped;\n  // The name of the input gridfunction is always "interp_sphgrid_MO_ETK::interped_gf":\n  const CCTK_STRING input_array_names[1] = { "interp_sphgrid_MO_ETK::interped_gf" };\n\n  // Perform interpolation!\n  for(int order=1; order <= 4; order *=2) {\n      char gf_name[100];\n      get_gf_name(*InterpCounter,gf_name);\n      printf("Interpolating\\033[1m %s \\033[0m... using interpolation order = %d\\n",gf_name,order);\n      Interpolate_to_sph_grid(cctkGH, N0*N1*N2, order,\n                                 points_x,points_y,points_z, input_array_names, output_f);\n\n      if(CCTK_MyProc(cctkGH)==0) {\n        for(int i=0;i<N0*N1*N2;i++) {\n            if(output_f[0][i] > 1e20) {\n                printf("BAD POINT: %s %d %e %e %e %e\\n",gf_name,i,points_x[i],points_y[i],points_z[i], output_f[0][i]);\n            }\n        }\n        output_to_file(CCTK_PASS_CTOC,gf_name,&order,output_f);\n        printf("Interpolate_to_sph_grid_main_function(): Just output to file at iteration %d\\n",cctk_iteration);\n      } else {\n        printf("Interpolate_to_sph_grid_main_function(): Process !=0 waiting for file output at iteration %d\\n",cctk_iteration);\n      }\n  }\n}\n')





import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import grid as gri               # NRPy+: Functions having to do with numerical grids
import finite_difference as fin  # NRPy+: Finite difference C code generation module
from outputC import lhrh         # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import loop as lp                # NRPy+: Generate C code loops

par.set_parval_from_str("grid::GridFuncMemAccess","ETK")

from collections import namedtuple
gf_interp = namedtuple('gf_interp', 'gf_description')
gf_interp_list = []
gf_interp_list.append(gf_interp("dummy -- used because this is a 1-offset array"))

interped_gf = gri.register_gridfunctions("AUX","interped_gf")

def interp_fileout(which_InterpCounter, expression, filename):
    kernel = fin.FD_outputC("returnstring",lhrh(lhs=gri.gfaccess("out_gfs","interped_gf"),rhs=expression),"outCverbose=False")
    output_type="a"
    if which_InterpCounter == 1:
        output_type="w"

    with open(filename, output_type) as file:
        file.write("if(*InterpCounter == "+str(which_InterpCounter)+") {\n")
        file.write(lp.loop(["i2","i1","i0"],
                           ["cctk_nghostzones[2]","cctk_nghostzones[1]","cctk_nghostzones[0]"],\
                           ["cctk_lsh[2]-cctk_nghostzones[2]",
                            "cctk_lsh[1]-cctk_nghostzones[1]",
                            "cctk_lsh[0]-cctk_nghostzones[0]"],\
                           ["1","1","1"],\
                           ["#pragma omp parallel for","",""],"   ",kernel))
        file.write("}\n")
    # If successful, return incremented which_InterpCounter:
    return which_InterpCounter+1





NRPyoutfilename = os.path.join(Ccodesdir,"src","list_of_functions_to_interpolate.h")

which_InterpCounter = 1





gammaDD = ixp.register_gridfunctions_for_single_rank2("AUX","gammaDD", "sym01")
betaU = ixp.register_gridfunctions_for_single_rank1("AUX","betaU")
alpha = gri.register_gridfunctions("AUX","alpha")

DIM=3

gf_interp_list.append(gf_interp("IGM density primitive"))
rho_b       = gri.register_gridfunctions("AUX","rho_b")
interp_expr = rho_b
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

gf_interp_list.append(gf_interp("IGM pressure primitive"))
P = gri.register_gridfunctions("AUX","P")
interp_expr = P
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)





IGMvU = ixp.register_gridfunctions_for_single_rank1("AUX","IGMvU")
Valenciav = ixp.zerorank1()
for i in range(DIM):
    Valenciav[i] = 1/alpha * (IGMvU[i] + betaU[i])
v_dot_v = sp.sympify(0)
for i in range(DIM):
    for j in range(DIM):
        v_dot_v += gammaDD[i][j]*Valenciav[i]*Valenciav[j]

Gamma_times_ValenciavU = ixp.zerorank1()
for i in range(DIM):
    Gamma_times_ValenciavU[i] = sp.sqrt(1/(1 - v_dot_v))*Valenciav[i]
    gf_interp_list.append(gf_interp("Lorentz factor, times Valencia vU"+str(i)))
    interp_expr = Gamma_times_ValenciavU[i]
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)



BU = ixp.register_gridfunctions_for_single_rank1("AUX","BU")
for i in range(DIM):
    gf_interp_list.append(gf_interp("IGM magnetic field component B"+str(i)))
    interp_expr = BU[i]
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)





betaD = ixp.zerorank1()
for i in range(DIM):
    for j in range(DIM):
        betaD[i] += gammaDD[i][j]*betaU[j]

beta2 = sp.sympify(0)
for i in range(DIM):
    beta2 += betaU[i]*betaD[i]

g4DD = ixp.zerorank2(DIM=4)
g4DD[0][0] = -alpha**2 + beta2
for i in range(DIM):
    g4DD[i+1][0] = g4DD[0][i+1] = betaD[i]
for i in range(DIM):
    for j in range(DIM):
        g4DD[i+1][j+1] = gammaDD[i][j]

for mu in range(4):
    for nu in range(mu,4):
        gf_interp_list.append(gf_interp("4-metric component g4DD"+str(mu)+str(nu)))
        interp_expr = g4DD[mu][nu]
        which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)





betaDdD = ixp.zerorank2()
gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01")
betaU_dD   = ixp.declarerank2("betaU_dD","nosym")
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            # Recall that betaD[i] = gammaDD[i][j]*betaU[j] (Eq. 2.121 in B&S)
            betaDdD[i][k] += gammaDD_dD[i][j][k]*betaU[j] + gammaDD[i][j]*betaU_dD[j][k]

g4DDdD = ixp.zerorank3(DIM=4)
alpha_dD   = ixp.declarerank1("alpha_dD")
for i in range(DIM):
    # Recall that g4DD[0][0] = -alpha^2 + betaU[i]*betaD[i]
    g4DDdD[0][0][i+1] += -2*alpha*alpha_dD[i]
    for j in range(DIM):
        g4DDdD[0][0][i+1] += betaU_dD[j][i]*betaD[j] + betaU[j]*betaDdD[j][i]

for i in range(DIM):
    for j in range(DIM):
        # Recall that g4DD[i][0] = g4DD[0][i] = betaD[i]
        g4DDdD[i+1][0][j+1] = g4DDdD[0][i+1][j+1] = betaDdD[i][j]
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            # Recall that g4DD[i][j] = gammaDD[i][j]
            g4DDdD[i+1][j+1][k+1] = gammaDD_dD[i][j][k]

gammaUU, dummyDET = ixp.symm_matrix_inverter3x3(gammaDD)

g4UU = ixp.zerorank2(DIM=4)
g4UU[0][0] = -1 / alpha**2
for i in range(DIM):
    g4UU[0][i+1] = g4UU[i+1][0] = betaU[i]/alpha**2
for i in range(DIM):
    for j in range(DIM):
        g4UU[i+1][j+1] = gammaUU[i][j] - betaU[i]*betaU[j]/alpha**2





Gamma4UDD = ixp.zerorank3(DIM=4)
for mu in range(4):
    for nu in range(4):
        for delta in range(4):
            for eta in range(4):
                Gamma4UDD[mu][nu][delta] += sp.Rational(1,2)*g4UU[mu][eta]*\
                (g4DDdD[eta][nu][delta] + g4DDdD[eta][delta][nu] - g4DDdD[nu][delta][eta])

for mu in range(4):
    for nu in range(4):
        for delta in range(nu,4):
            gf_interp_list.append(gf_interp("4-Christoffel GammaUDD"+str(mu)+str(nu)+str(delta)))
            interp_expr = Gamma4UDD[mu][nu][delta]
            which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)





get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/construct_function_to_interpolate__store_to_interped_gf.cc', '#include <stdio.h>\n#include <stdlib.h>\n#include "cctk.h"\n#include "cctk_Arguments.h"\n#include "cctk_Parameters.h"\n\n// Set the gridfunction interped_gf, according to the interpolation counter variable interp_counter.\n//    For example, we might interpolate "IllinoisGRMHD::rho_b" if interp_counter==0. The following\n//    function takes care of these\nvoid list_of_functions_to_interpolate(cGH *cctkGH,const CCTK_INT *cctk_lsh,const CCTK_INT *cctk_nghostzones,\n                                     const CCTK_REAL invdx0,const CCTK_REAL invdx1,const CCTK_REAL invdx2,\n                                     const CCTK_INT *InterpCounter,\n                                     const CCTK_REAL *rho_bGF,const CCTK_REAL *PGF,\n                                     const CCTK_REAL *IGMvU0GF,const CCTK_REAL *IGMvU1GF,const CCTK_REAL *IGMvU2GF,\n                                     const CCTK_REAL *BU0GF,const CCTK_REAL *BU1GF,const CCTK_REAL *BU2GF,\n                                     const CCTK_REAL *gammaDD00GF,const CCTK_REAL *gammaDD01GF,const CCTK_REAL *gammaDD02GF,\n                                     const CCTK_REAL *gammaDD11GF,const CCTK_REAL *gammaDD12GF,const CCTK_REAL *gammaDD22GF,\n                                     const CCTK_REAL *betaU0GF,const CCTK_REAL *betaU1GF,const CCTK_REAL *betaU2GF,\n                                     const CCTK_REAL *alphaGF,   CCTK_REAL *interped_gfGF) {\n#include "list_of_functions_to_interpolate.h"\n}\n\nvoid construct_function_to_interpolate__store_to_interped_gf(CCTK_ARGUMENTS) {\n  DECLARE_CCTK_ARGUMENTS;\n  DECLARE_CCTK_PARAMETERS;\n  const CCTK_REAL invdx0 = 1.0 / CCTK_DELTA_SPACE(0);\n  const CCTK_REAL invdx1 = 1.0 / CCTK_DELTA_SPACE(1);\n  const CCTK_REAL invdx2 = 1.0 / CCTK_DELTA_SPACE(2);\n  list_of_functions_to_interpolate(cctkGH,cctk_lsh,cctk_nghostzones,invdx0,invdx1,invdx2,\n                                   InterpCounter,\n                                   rho_b,P,\n                                   vx,vy,vz,\n                                   Bx,By,Bz,\n                                   gxx,gxy,gxz,gyy,gyz,gzz,\n                                   betax,betay,betaz,alp, interped_gf);\n// interped_gf will be interpolated across AMR boundaries, meaning that\n//    it must be prolongated. Only gridfunctions with 3 timelevels stored\n//    may be prolongated (provided time_interpolation_order is set to the\n//    usual value of 2). We should only call this interpolation routine\n//    at iterations in which all gridfunctions are on the same timelevel\n//    (usually a power of 2), which will ensure that the following\n//    "filling of the timelevels" is completely correct.\n#pragma omp parallel for\n    for(int i=0;i<cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];i++) {\n        interped_gf_p[i]   = interped_gf[i];\n        interped_gf_p_p[i] = interped_gf[i];\n    }\n}\n')





with open(os.path.join(Ccodesdir,"src","get_gf_name.h"), "w") as file:
    file.write("void get_gf_name(const int InterpCounter,char gf_name[100]) {\n")
    for i in range(1,which_InterpCounter):
        file.write("    if(InterpCounter=="+str(i)+") { snprintf(gf_name,100,\""+gf_interp_list[i].gf_description+"\"); return; }\n")
    file.write("    printf(\"Error. InterpCounter = %d unsupported. I should not be here.\\n\",InterpCounter); exit(1);\n")
    file.write("}\n")





with open(os.path.join(Ccodesdir,"src","define_NumInterpFunctions.h"), "w") as file:
    file.write("#define NumInterpFunctions "+str(which_InterpCounter)+"\n")




get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/interp_counter.cc', '#include <assert.h>\n#include <stdio.h>\n#include <stdlib.h>\n#include <string.h>\n#include <math.h>\n#include <ctype.h>\n#include "cctk.h"\n#include "cctk_Arguments.h"\n#include "cctk_Parameters.h"\n\n#include "define_NumInterpFunctions.h"\n\nvoid SphGrid_InitializeInterpCounterToZero(CCTK_ARGUMENTS)\n{\n  DECLARE_CCTK_ARGUMENTS;\n  DECLARE_CCTK_PARAMETERS;\n  *InterpCounter = 0;\n\n  if(verbose==2) printf("interp_sphgrid_MO_ETK: Just set InterpCounter to %d\\n",*InterpCounter);\n}\n\nvoid SphGrid_InitializeInterpCounter(CCTK_ARGUMENTS)\n{\n  DECLARE_CCTK_ARGUMENTS;\n  DECLARE_CCTK_PARAMETERS;\n\n  if(cctk_iteration == interp_out_iteration) {\n    *InterpCounter = 1;\n    if(verbose==2) printf("interp_sphgrid_MO_ETK: Just set InterpCounter to %d ; ready to start looping over interpolated gridfunctions!\\n",\n                          *InterpCounter);\n  }\n}\n\n// This function increments InterpCounter if we are at the interp_out_iteration until\n// it hits NumInterpFunctions. At this iteration, InterpCounter is set to zero, which\n// exits the loop.\nvoid SphGrid_IncrementInterpCounter(CCTK_ARGUMENTS)\n{\n  DECLARE_CCTK_ARGUMENTS;\n  DECLARE_CCTK_PARAMETERS;\n\n    if(*InterpCounter == NumInterpFunctions-1) {\n        *InterpCounter = 0;\n        if(verbose==2) printf("interp_sphgrid_MO_ETK: Finished! Just zeroed InterpCounter.\\n");\n    } else {\n        (*InterpCounter)++;\n        if(verbose==2) printf("interp_sphgrid_MO_ETK: Just incremented InterpCounter to %d of %d\\n",*InterpCounter,NumInterpFunctions-1);\n    }\n}\n')






get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/make.code.defn', '# Main make.code.defn file for thorn interp_sphgrid_MO_ETK\n\n# Source files in this directory\nSRCS =  main_function.cc interp_counter.cc construct_function_to_interpolate__store_to_interped_gf.cc\n')





get_ipython().run_cell_magic('writefile', '$Ccodesdir/interface.ccl', '\n# With "implements", we give our thorn its unique name.\nimplements: interp_sphgrid_MO_ETK\n\n# By "inheriting" other thorns, we tell the Toolkit that we\n#   will rely on variables/function that exist within those\n#   functions.\ninherits:   admbase IllinoisGRMHD Grid\n\n# Tell the Toolkit that we want "interped_gf" and "InterpCounter"\n#    and invariants to NOT be visible to other thorns, by using\n#    the keyword "private". Note that declaring these\n#    gridfunctions here *does not* allocate memory for them;\n#    that is done by the schedule.ccl file.\nprivate:\nCCTK_REAL interpolation_gf type=GF timelevels=3 tags=\'Checkpoint="no"\'\n{\n  interped_gf\n} "Gridfunction containing output from interpolation."\n\nint InterpCounterVar type = SCALAR tags=\'checkpoint="no"\'\n{\n  InterpCounter\n} "Counter that keeps track of which function we are interpolating."\n\nCCTK_REAL interp_pointcoords_and_output_arrays TYPE=ARRAY DISTRIB=CONSTANT DIM=1 SIZE=N0*N1*N2 tags=\'checkpoint="no"\'\n{\n  points_x,points_y,points_z,\n  output_interped\n}\n')





get_ipython().run_cell_magic('writefile', '$Ccodesdir/param.ccl', '\n# Output the interpolated data to the IO::out_dir directory:\nshares: IO\nUSES STRING out_dir\n\nrestricted:\n\n########################################\n# BASIC THORN STEERING PARAMETERS\nCCTK_INT interp_out_iteration "Which iteration to interpolate to spherical grids?" STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 960000\n\n## Interpolator information\nCCTK_STRING interpolator_name "Which interpolator to use?" STEERABLE=ALWAYS\n{\n  ".+" :: "Any nonempty string; an unsupported value will throw an error."\n} "Lagrange polynomial interpolation"\n\nCCTK_INT verbose "Set verbosity level: 1=useful info; 2=moderately annoying (though useful for debugging)" STEERABLE=ALWAYS\n{\n  0:2 :: "0 = no output; 1=useful info; 2=moderately annoying (though useful for debugging)"\n} 2\n########################################\n# SPHERICAL COORDINATE SYSTEM PARAMETERS\nCCTK_INT N0 "Number of points in r direction" STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 96\n\nCCTK_INT N1 "Number of points in theta direction" STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 96\n\nCCTK_INT N2 "Number of points in phi direction" STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 96\n\n##########\n# Cartesian position of center of spherical grid (usually center of BH) -- CURRENTLY UNSUPPORTED!\nCCTK_REAL x_center "x-position of center." STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 0.0\n\nCCTK_REAL y_center "y-position of center." STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 0.0\n\nCCTK_REAL z_center "z-position of center." STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 0.0\n\n##########\n# Radial parameters:\nCCTK_REAL R0 "Radial offset: r(x0) = R_0 + exp(x0). Probably should keep it set to zero." STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 0.0\n\nCCTK_REAL Rin "x0 offset: x0 = log(Rin-R0) + (i + 0.5)Dx0." STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 1.08986052555408\n\nCCTK_REAL Rout "Dx0 = log( (Rout-R0) / (Rin-R0) )/N0" STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 80.0\n\n##########\n# Theta parameters:\nCCTK_REAL x1_beg "x1 offset: x1 = x1_beg + (j + 0.5)Dx1. Probably should keep it set to zero." STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 0.0\n\nCCTK_INT theta_option "Which prescription for theta should be used? 1 or 2?" STEERABLE=ALWAYS\n{\n  1:2 :: ""\n} 1\n\nCCTK_REAL th_c "theta_c: Angular cutout size for theta = 0 and pi" STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 0.053407075111026485 # 0.017*pi\n\nCCTK_REAL xi "Amplitude of nonlinear part of the theta distribution." STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 0.25\n\nCCTK_INT th_n "Power of nonlinear part of theta distribution. Only for theta_option=2" STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 9\n\n##########\n# Phi parameters:\nCCTK_REAL x2_beg "x2 offset: x2 = x2_beg + (k + 0.5)Dx2. Probably should keep it set to zero." STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 0.0\n########################################\n')





get_ipython().run_cell_magic('writefile', '$Ccodesdir/schedule.ccl', '\nSTORAGE: interpolation_gf[3]\nSTORAGE: InterpCounterVar\nSTORAGE: interp_pointcoords_and_output_arrays\n\n#############################\nSCHEDULE SphGrid_InitializeInterpCounterToZero AT CCTK_INITIAL\n{\n  LANG: C\n  OPTIONS: GLOBAL\n} "Initialize InterpCounter variable to zero"\n\nSCHEDULE SphGrid_InitializeInterpCounterToZero AT CCTK_POST_RECOVER_VARIABLES\n{\n  LANG: C\n  OPTIONS: GLOBAL\n} "Initialize InterpCounter variable to zero"\n\nSCHEDULE SphGrid_InitializeInterpCounter before SphGrid_InterpGroup AT CCTK_ANALYSIS\n{\n  LANG: C\n  OPTIONS: GLOBAL\n} "Initialize InterpCounter variable"\n##################\n\nSCHEDULE GROUP SphGrid_InterpGroup AT CCTK_ANALYSIS BEFORE CarpetLib_printtimestats BEFORE CarpetLib_printmemstats AFTER Convert_to_HydroBase WHILE interp_sphgrid_MO_ETK::InterpCounter\n{\n} "Perform all spherical interpolations. This group is only actually scheduled at cctk_iteration==interp_out_iteration."\n\nSCHEDULE construct_function_to_interpolate__store_to_interped_gf in SphGrid_InterpGroup before DoSum\n{\n  STORAGE: interpolation_gf[3],InterpCounterVar,interp_pointcoords_and_output_arrays\n  OPTIONS: GLOBAL,LOOP-LOCAL\n  SYNC: interpolation_gf\n  LANG: C\n} "Construct the function to interpolate"\n\nSCHEDULE Interpolate_to_sph_grid_main_function in SphGrid_InterpGroup after construct_function_to_interpolate__store_to_interped_gf\n{\n  OPTIONS: GLOBAL\n  LANG: C\n} "Perform interpolation and output result to file."\n#######\nSCHEDULE SphGrid_IncrementInterpCounter in SphGrid_InterpGroup after Interpolate_to_sph_grid_main_function\n{\n  LANG: C\n  OPTIONS: GLOBAL\n} "Increment InterpCounter variable, or set to zero once loop is complete."\n##################\n')






import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ETK_thorn-Interpolation_to_Spherical_Grids_multi_order")

