#!/usr/bin/env python






import cmdline_helper as cmd  # NRPy+: Multi-platform Python command-line interface
import shutil, os             # Standard Python modules for multiplatform OS-level functions, benchmarking

Ccodesdir = "interp_arbgrid_MO_ETK"
shutil.rmtree(Ccodesdir, ignore_errors=True)
cmd.mkdir(Ccodesdir)
cmd.mkdir(os.path.join(Ccodesdir,"src/"))
cmd.mkdir(os.path.join(Ccodesdir,"src","standalone/"))





get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/Interpolate_to_dest_grid.h', '\nvoid Interpolate_to_dest_grid(const cGH *restrict cctkGH,\n                              const CCTK_INT interp_num_points,\n                              const CCTK_INT interp_order,\n                              char interpolator_name[100],\n                              const CCTK_REAL *restrict point_x_temp,\n                              const CCTK_REAL *restrict point_y_temp,\n                              const CCTK_REAL *restrict point_z_temp,\n                              const CCTK_STRING input_array_names[1],\n                              CCTK_REAL *output_f[1]) {\n  DECLARE_CCTK_PARAMETERS;\n  CCTK_INT ierr;\n  const CCTK_INT NUM_INPUT_ARRAYS=1;\n  const CCTK_INT NUM_OUTPUT_ARRAYS=1;\n\n  CCTK_STRING coord_system = "cart3d";\n\n  // Set up handles\n  const CCTK_INT coord_system_handle = CCTK_CoordSystemHandle(coord_system);\n  if (coord_system_handle < 0) {\n    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,\n               "can\'t get coordinate system handle for coordinate system \\"%s\\"!",\n               coord_system);\n  }\n\n  const CCTK_INT operator_handle = CCTK_InterpHandle(interpolator_name);\n  if (operator_handle < 0)\n    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,\n               "couldn\'t find interpolator \\"%s\\"!",\n               interpolator_name);\n\n  char interp_order_string[10];\n  snprintf(interp_order_string, 10, "order=%d", interp_order);\n  CCTK_STRING interpolator_pars = interp_order_string;\n  CCTK_INT param_table_handle = Util_TableCreateFromString(interpolator_pars);\n  if (param_table_handle < 0) {\n    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,\n               "bad interpolator parameter(s) \\"%s\\"!",\n               interpolator_pars);\n  }\n\n  CCTK_INT operand_indices[NUM_INPUT_ARRAYS]; //NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS];\n  for(int i = 0 ; i < NUM_INPUT_ARRAYS  ; i++) {\n    operand_indices[i] = i;\n  }\n  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,\n                        operand_indices, "operand_indices");\n\n  CCTK_INT operation_codes[NUM_INPUT_ARRAYS];\n  for(int i = 0 ; i < NUM_INPUT_ARRAYS  ; i++) {\n    operation_codes[i] = 0;\n  }\n  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,\n                        operation_codes, "operation_codes");\n\n  const void* interp_coords[3] = { (const void *) point_x_temp,\n                                   (const void *) point_y_temp,\n                                   (const void *) point_z_temp };\n\n  CCTK_INT input_array_indices[NUM_INPUT_ARRAYS];\n  for(int i = 0 ; i < NUM_INPUT_ARRAYS ; i++) {\n    input_array_indices[i] = CCTK_VarIndex(input_array_names[i]);\n    if(input_array_indices[i] < 0) {\n      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,\n                 "COULD NOT FIND VARIABLE \'%s\'.",\n                 input_array_names[i]);\n      exit(1);\n    }\n  }\n\n  CCTK_INT output_array_types[NUM_OUTPUT_ARRAYS];\n  for(int i = 0 ; i < NUM_OUTPUT_ARRAYS ; i++) {\n    output_array_types[i] = CCTK_VARIABLE_REAL;\n  }\n\n  void* output_arrays[NUM_OUTPUT_ARRAYS] = { (void *) output_f[0] };\n\n  // actual interpolation call\n  ierr = CCTK_InterpGridArrays(cctkGH,\n                               3, // number of dimensions\n                               operator_handle,\n                               param_table_handle,\n                               coord_system_handle,\n                               interp_num_points,\n                               CCTK_VARIABLE_REAL,\n                               interp_coords,\n                               NUM_INPUT_ARRAYS, // Number of input arrays\n                               input_array_indices,\n                               NUM_OUTPUT_ARRAYS, // Number of output arrays\n                               output_array_types,\n                               output_arrays);\n\n  if (ierr<0) {\n    CCTK_WARN(1,"interpolation screwed up");\n    Util_TableDestroy(param_table_handle);\n    exit(1);\n  }\n\n  ierr = Util_TableDestroy(param_table_handle);\n  if (ierr != 0) {\n    CCTK_WARN(1,"Could not destroy table");\n    exit(1);\n  }\n}\n')




get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/interpolate_set_of_points_in_file.h', '\n#define ALLOCATE_2D_GENERIC(type,array,ni,nj) type **array=(type **)malloc(ni * sizeof(type *)); \\\n  for(int cc = 0; cc < ni; cc++)                 array[cc]=(type * )malloc(nj * sizeof(type));\n#define FREE_2D_GENERIC(type,array,ni,nj) for(int cc = 0; cc < ni;cc++) free((void *)array[cc]); \\\n  /**/                                                                  free((void *)array);\n#include "output_to_file.h"\n\n// Calls the above function and output_to_file().\nvoid interpolate_set_of_points_in_file(CCTK_ARGUMENTS,char filename_basename[100],char gf_name[100],char interpolator_name[100],int num_interp_orders,int *interp_orders_list) {\n\n  DECLARE_CCTK_ARGUMENTS;\n  DECLARE_CCTK_PARAMETERS; // Needed for x_center,y_center,z_center\n  // Set up output array:\n  // The name of the input gridfunction is always "interp_arbgrid_MO_ETK::interped_gf":\n  const CCTK_STRING input_array_names[1] = { "interp_arbgrid_MO_ETK::interped_gf" };\n  CCTK_REAL *points_dest_grid_x,*points_dest_grid_y,*points_dest_grid_z; // Coordinates of points of destination grid\n  CCTK_REAL **output_f; // Output to be written to dataset, will be filled with NaNs at points out of bounds\n  // For benchmarking purposes:\n  time_t start_timer,end_timer;\n  time(&start_timer); // Resolution of one second...\n  CCTK_REAL time_in_seconds;\n\n  int num_dest_grid_points;\n  if(CCTK_MyProc(cctkGH)==0) {\n    // Step 1: Read list of desired interpolation destination points from file:\n\n    // Step 1.a: Read integer at top of file indicating number of points.\n    int num_dest_grid_pointsx,num_dest_grid_pointsy,num_dest_grid_pointsz;\n    char pointsx_filename[100]; snprintf(pointsx_filename,100,"%s-x.dat",filename_basename);\n    printf("Reading list of x data points from file %s...\\n",pointsx_filename);\n    FILE *pointsx_file = fopen(pointsx_filename, "rb");\n    if(!pointsx_file) { printf("Error: Unable to open %s\\n",pointsx_filename); exit(1); }\n    fread(&num_dest_grid_pointsx, sizeof(int), 1, pointsx_file);\n\n    char pointsy_filename[100]; snprintf(pointsy_filename,100,"%s-y.dat",filename_basename);\n    printf("Reading list of y data points from file %s...\\n",pointsy_filename);\n    FILE *pointsy_file = fopen(pointsy_filename, "rb");\n    if(!pointsy_file) { printf("Error: Unable to open %s\\n",pointsy_filename); exit(1); }\n    fread(&num_dest_grid_pointsy, sizeof(int), 1, pointsy_file);\n\n    char pointsz_filename[100]; snprintf(pointsz_filename,100,"%s-z.dat",filename_basename);\n    printf("Reading list of z data points from file %s...\\n",pointsz_filename);\n    FILE *pointsz_file = fopen(pointsz_filename, "rb");\n    if(!pointsz_file) { printf("Error: Unable to open %s\\n",pointsz_filename); exit(1); }\n    fread(&num_dest_grid_pointsz, sizeof(int), 1, pointsz_file);\n\n    // Step 1.a.i: Sanity check: make sure that num_dest_grid_pointsx == num_dest_grid_pointsy == num_dest_grid_pointsz\n    if(num_dest_grid_pointsx != num_dest_grid_pointsy || num_dest_grid_pointsy != num_dest_grid_pointsz) {\n      printf("Error: Failed sanity check. Number of interpolation points different in %s-{x,y,z}.dat data files!\\n",\n             filename_basename);\n      exit(1);\n    } else {\n      // If sanity check passes:\n      num_dest_grid_points = num_dest_grid_pointsx;\n    } // END sanity check\n\n      // Step 1.b: Allocate memory for destination grids and interpolation output\n    if(num_dest_grid_points <= 0 || num_dest_grid_points > 2000000000) {\n      printf("Error: Failed sanity check. Number of interpolation points was found to be: %d",num_dest_grid_points);\n      exit(1);\n    } // END sanity check\n    points_dest_grid_x = (CCTK_REAL  *)malloc(sizeof(CCTK_REAL)*num_dest_grid_points);\n    points_dest_grid_y = (CCTK_REAL  *)malloc(sizeof(CCTK_REAL)*num_dest_grid_points);\n    points_dest_grid_z = (CCTK_REAL  *)malloc(sizeof(CCTK_REAL)*num_dest_grid_points);\n    output_f = (CCTK_REAL **)malloc(1 * sizeof(CCTK_REAL *));\n    for(int cc = 0; cc < 1; cc++) output_f[cc]=(CCTK_REAL *)malloc(num_dest_grid_points * sizeof(CCTK_REAL));\n\n    // Step 1.c: Store cell-centered points to allocated memory.\n    fread(points_dest_grid_x, sizeof(CCTK_REAL), num_dest_grid_points, pointsx_file);\n    fread(points_dest_grid_y, sizeof(CCTK_REAL), num_dest_grid_points, pointsy_file);\n    fread(points_dest_grid_z, sizeof(CCTK_REAL), num_dest_grid_points, pointsz_file);\n    int magic_numberx; fread(&magic_numberx, sizeof(int), 1, pointsx_file);\n    int magic_numbery; fread(&magic_numbery, sizeof(int), 1, pointsy_file);\n    int magic_numberz; fread(&magic_numberz, sizeof(int), 1, pointsz_file);\n    int correct_magicnum = -349289480;\n    if(magic_numberx != correct_magicnum || magic_numbery != correct_magicnum || magic_numberz != correct_magicnum) {\n      printf("Error: Failed sanity check. Magic numbers in x,y,z data files were: %d %d %d, respectively, but should have been: %d",\n             magic_numberx,magic_numbery,magic_numberz,correct_magicnum);\n      exit(1);\n    }\n    fclose(pointsx_file);\n    fclose(pointsy_file);\n    fclose(pointsz_file);\n    time(&end_timer); time_in_seconds = difftime(end_timer,start_timer); time(&start_timer);\n    printf("Finished in %e seconds.\\n",time_in_seconds);\n\n    // Step 1.d: Apply offset to x,y,z coordinates to ensure they are centered on (x_center,y_center,z_center)\n    //           For example, if a black hole is situated at (x,y,z) = (1,2,3), then we set\n    //           (x_center,y_center,z_center) = (1,2,3) in our ETK parameter file (i.e., with extension .par)\n    //           and if we desire a point at (x_dest,y_dest,z_dest) = (0,0,0) on the *destination* grid,\n    //           this will correspond to point (x_src,y_src,z_src) = (1,2,3) = (x_center,y_center,z_center)\n    //           on the source grid. Thus the translation between source and destination grids is given by\n    //           (x_src,y_src,z_src) = (x_dest+x_center, y_dest+y_center, z_dest+z_center),\n    //           where (x_src,y_src,z_src) = (points_dest_grid_x[i],points_dest_grid_y[i],points_dest_grid_z[i]) for point i.\n    for(int point=0;point<num_dest_grid_points;point++) {\n      points_dest_grid_x[point] += x_center;\n      points_dest_grid_y[point] += y_center;\n      points_dest_grid_z[point] += z_center;\n    }\n\n    // Step 1.e: Look for the points that are out of bounds and set their coordinates to out_of_bounds_interp_xyz.\n    //           At the end, we will replace the interpolation output at these points with NaNs.\n    printf("Looking for points that are out of bounds.\\n");\n    int num_out_of_bounds_points = 0;\n    for(int i=0;i<num_dest_grid_points;i++){\n      if(fabs(points_dest_grid_x[i])>out_of_bounds_interp_xyz ||\n         fabs(points_dest_grid_y[i])>out_of_bounds_interp_xyz ||\n         fabs(points_dest_grid_z[i])>out_of_bounds_interp_xyz) {\n        points_dest_grid_x[i] = out_of_bounds_interp_xyz;\n        points_dest_grid_y[i] = out_of_bounds_interp_xyz;\n        points_dest_grid_z[i] = out_of_bounds_interp_xyz;\n        num_out_of_bounds_points++;\n      }\n    }\n    time(&end_timer); time_in_seconds = difftime(end_timer,start_timer); time(&start_timer);\n    printf("Finished in %e seconds, found %i points out of bounds.\\n", time_in_seconds, num_out_of_bounds_points);\n\n  } // END if(CCTK_MyProc(cctkGH)==0)\n\n\n  // Step 1.f: Looping over interp order as desired, interpolate to destination points & output to file\n  for(int order_i=0; order_i<num_interp_orders; order_i++) {\n    int order = interp_orders_list[order_i];\n    printf("Interpolating\\033[1m %s \\033[0m... using interpolation order = %d\\n",gf_name,order);\n    printf("and %s \\n",interpolator_name);\n    if(CCTK_MyProc(cctkGH)==0) {\n      Interpolate_to_dest_grid(cctkGH, num_dest_grid_points, order, interpolator_name,\n                               points_dest_grid_x,points_dest_grid_y,points_dest_grid_z, input_array_names, output_f);\n\n      // Step 1.f.i: Sanity check -- check for bad point:\n#pragma omp parallel for\n      for(int i=0;i<num_dest_grid_points;i++) {\n        if(output_f[0][i] > 1e20) {\n          printf("BAD POINT: %s %d %e %e %e %e\\n",gf_name,i,points_dest_grid_x[i],points_dest_grid_y[i],points_dest_grid_z[i], output_f[0][i]);\n          exit(1);\n        }\n      }\n      time(&end_timer); time_in_seconds = difftime(end_timer,start_timer); time(&start_timer);\n      printf("Finished in %d seconds. Next: Filling cells out of bounds with NaNs\\n",time_in_seconds);\n      // Step 1.f.ii: Filling cells out of bounds with NaNs:\n      for(int i=0;i<num_dest_grid_points;i++){\n        if(fabs(points_dest_grid_x[i])==out_of_bounds_interp_xyz ||\n           fabs(points_dest_grid_y[i])==out_of_bounds_interp_xyz ||\n           fabs(points_dest_grid_z[i])==out_of_bounds_interp_xyz) {\n          output_f[0][i] = 0.0 / 0.0; // ie NaN\n        }\n      }\n      time(&end_timer); time_in_seconds = difftime(end_timer,start_timer); time(&start_timer);\n      printf("Finished in %e seconds. Next: Interpolate_to_dest_grid_main_function(): Outputting to file at iteration %d\\n",time_in_seconds,cctk_iteration);\n      output_to_file(CCTK_PASS_CTOC,gf_name,&order,&num_dest_grid_points,output_f);\n      time(&end_timer); time_in_seconds = difftime(end_timer,start_timer); time(&start_timer);\n      printf("Finished in %e seconds. Interpolate_to_dest_grid_main_function(): Finished output to file at iteration %d\\n",time_in_seconds,cctk_iteration);\n    } else {\n      // On all MPI processes that are nonzero, only call the interpolation function\n      //    to ensure the MPI calls from the actual interpolation (driven by proc==0) are seen.\n      // Setting num_dest_grid_points to zero results in a segfault on certain (ahem, Frontera)\n      //    systems. So we set num_dest_grid_points = 1 and interpolate to the origin\n      //    only for MPI processes that are nonzero, leaving the heavy lifting to MPI process 0.\n      num_dest_grid_points = 1;\n      points_dest_grid_x = (CCTK_REAL  *)malloc(sizeof(CCTK_REAL)*num_dest_grid_points);\n      points_dest_grid_y = (CCTK_REAL  *)malloc(sizeof(CCTK_REAL)*num_dest_grid_points);\n      points_dest_grid_z = (CCTK_REAL  *)malloc(sizeof(CCTK_REAL)*num_dest_grid_points);\n      output_f = (CCTK_REAL **)malloc(1 * sizeof(CCTK_REAL *));\n      for(int cc = 0; cc < 1; cc++) output_f[cc]=(CCTK_REAL *)malloc(num_dest_grid_points * sizeof(CCTK_REAL));\n\n      Interpolate_to_dest_grid(cctkGH, num_dest_grid_points, order, interpolator_name,\n                               points_dest_grid_x,points_dest_grid_y,points_dest_grid_z, input_array_names, output_f);\n      output_f[0][0] = 0.0;\n    } // END if(CCTK_MyProc(cctkGH)==0)\n  } // END for(int order_i=0; order_i<num_interp_orders; order_i++)\n\n  // Step 1.g: Free memory for destination grids and interpolation output\n  free(points_dest_grid_x);\n  free(points_dest_grid_y);\n  free(points_dest_grid_z);\n  FREE_2D_GENERIC(CCTK_REAL,output_f,1,num_dest_grid_points);\n} // END function\n#undef ALLOCATE_2D_GENERIC\n#undef FREE_2D_GENERIC\n')





get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/output_to_file.h', '\n#include "define_NumInterpFunctions.h"\n\n// output_to_file() starts order and InterpCounter both with the value 1\nvoid output_to_file(CCTK_ARGUMENTS,\n                    char gf_name[100],\n                    CCTK_INT *restrict order,\n                    CCTK_INT *restrict num_interp_points,\n                    CCTK_REAL *restrict output_f[1]) {\n  DECLARE_CCTK_ARGUMENTS;\n  DECLARE_CCTK_PARAMETERS;\n\n  char filename[100];\n  sprintf (filename, "%s/interp_dest_grids_MO.dat", out_dir);\n  FILE *file;\n  if(*InterpCounter == 1 && *order==1) {\n    file = fopen (filename,"w");\n    printf("WRITING to file %s\\n",filename);\n    // Write EOS information to the beginning of the file\n    char eos_info[256];\n    if( enable_nuc_eos ) {\n      // If using nuclear EOS, give the table name\n      sprintf(eos_info,"Nuclear (tabulated) EOS from file %s",nuceos_table_name);\n    }\n    else {\n      // If using hybrid EOS, then give Gamma_th and the number of polytropic pieces\n      sprintf(eos_info,"Hybrid EOS with Gamma_th = %.15e and %d polytropic pieces",Gamma_th,neos);\n    }\n    fwrite(eos_info, 256*sizeof(char), 1, file);\n  }\n  else {\n    file = fopen (filename,"a+");\n    printf("Appending to file %s\\n",filename);\n  }\n  if (! file) {\n    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,\n                "interp_dest_grid__ET_thorn: Cannot open output file \'%s\'", filename);\n    exit(1);\n  }\n\n  fwrite(gf_name, 100*sizeof(char), 1, file);\n  fwrite(order, sizeof(CCTK_INT), 1, file);\n  fwrite(num_interp_points, sizeof(CCTK_INT),1,file);\n\n  CCTK_REAL magic_number = 1.130814081305130e-21;\n  fwrite(&magic_number, sizeof(CCTK_REAL), 1, file);\n  fwrite(&cctk_iteration, sizeof(CCTK_INT), 1, file);\n  fwrite(&cctk_time, sizeof(CCTK_REAL), 1, file);\n  for(CCTK_INT i=0;i<1;i++) {\n    fwrite(output_f[i], sizeof(CCTK_REAL)*(*num_interp_points), 1, file);\n  }\n\n  fclose(file);\n}\n')





get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/main_function.cc', '\n\n// Include needed ETK & C library header files:\n#include <stdio.h>\n#include <stdlib.h>\n#include <string.h>\n#include <math.h>\n#include <time.h> // for benchmarking\n// Needed for dealing with Cactus/ETK infrastructure\n#include "cctk.h"\n#include "cctk_Arguments.h"\n#include "cctk_Parameters.h"\n// Needed for low-level interpolation functions\n#include "util_Table.h"\n#include "util_String.h"\n\n// Include locally-defined C++ functions:\n#include "Interpolate_to_dest_grid.h"\n#include "get_gf_name.h"\n#include "interpolate_set_of_points_in_file.h"\n\nvoid Interpolate_to_dest_grid_main_function(CCTK_ARGUMENTS) {\n  DECLARE_CCTK_ARGUMENTS;\n  DECLARE_CCTK_PARAMETERS;\n\n  // Perform interpolation only at iteration == interp_out_iteration:\n  if(cctk_iteration != interp_out_iteration) return;\n  // Perform interpolation!\n  // Process zero (CCTK_MyProc(cctkGH)==0) is responsible for directing the interpolation.\n  //    All other processes must see the cctk_InterpGridArrays() within Interpolate_to_dest_grid(),\n  //    so that the MPI calls work properly, but these nonzero processes can call\n  //    Interpolate_to_dest_grid() with number of interpolated points set to zero, and\n  //    without needing a malloc().\n  char gf_name[100]; get_gf_name(*InterpCounter,gf_name);\n  char filename_basename[100];\n  char interpolator_name[100];\n\n  // The following "if" statement is destination-code dependent.\n  // In our case, since we are interpolating variables for harm3d,\n  // we interpolate the vector potential to the corners of each cell.\n  // Every other quantity is interpolated to the center of each cell.\n if(strncmp(gf_name,"Unstaggered",11) == 0){\n    sprintf(filename_basename,"corner_points");\n  }\n  else{\n    sprintf(filename_basename,"cell_centered_points");\n  }\n\n  int num_interp_orders,*interp_orders_list;\n  // 4-metric, 4-Christoffels and A_mu only output Hermite interpolation order==3\n  // so we secure continuity of first derivative.\n  if( (strncmp(gf_name,"4-",2) == 0) ||\n      (strncmp(gf_name,"Unstaggered",11) == 0) ||\n     (strncmp(gf_name,"Constraints",11) == 0)) {\n    num_interp_orders = 1;\n    sprintf(interpolator_name, "Hermite polynomial interpolation");\n    interp_orders_list = (int *)malloc(sizeof(int)*num_interp_orders);\n    interp_orders_list[0] = 3;\n  } else {\n    num_interp_orders = 3;\n    sprintf(interpolator_name, "Lagrange polynomial interpolation");\n    interp_orders_list = (int *)malloc(sizeof(int)*num_interp_orders);\n    int count = 0; for(int order=1;order<=4;order*=2) { interp_orders_list[count] = order; count++; }\n  }\n\n  interpolate_set_of_points_in_file(CCTK_PASS_CTOC,filename_basename,gf_name,interpolator_name,num_interp_orders,interp_orders_list);\n  free(interp_orders_list);\n\n  // Now perform interpolation of 4-metric on\n  //   faces (i-1/2,j,k), (i,j-1/2,k), (i,j,k-1/2) and corners (i-1/2,j-1/2,k-1/2)\n  if(strncmp(gf_name,"4-metric",8) == 0) {\n    num_interp_orders = 1;\n    sprintf(interpolator_name, "Hermite polynomial interpolation");\n    interp_orders_list = (int *)malloc(sizeof(int)*num_interp_orders);\n    interp_orders_list[0] = 3;\n\n    char gf_name_new[100];\n\n    sprintf(filename_basename,"faceim_points");\n    snprintf(gf_name_new,100,"faceim (i-1/2,j,k): %s",gf_name);\n    interpolate_set_of_points_in_file(CCTK_PASS_CTOC,filename_basename,gf_name_new,interpolator_name,num_interp_orders,interp_orders_list);\n\n    sprintf(filename_basename,"facejm_points");\n    snprintf(gf_name_new,100,"facejm (i,j-1/2,k): %s",gf_name);\n    interpolate_set_of_points_in_file(CCTK_PASS_CTOC,filename_basename,gf_name_new,interpolator_name,num_interp_orders,interp_orders_list);\n\n    sprintf(filename_basename,"facekm_points");\n    snprintf(gf_name_new,100,"facekm (i,j,k-1/2): %s",gf_name);\n    interpolate_set_of_points_in_file(CCTK_PASS_CTOC,filename_basename,gf_name_new,interpolator_name,num_interp_orders,interp_orders_list);\n\n    sprintf(filename_basename,"corner_points");\n    snprintf(gf_name_new,100,"cornr (i-1/2,j-1/2,k-1/2): %s",gf_name);\n    interpolate_set_of_points_in_file(CCTK_PASS_CTOC,filename_basename,gf_name_new,interpolator_name,num_interp_orders,interp_orders_list);\n  } // END if(strncmp(gf_name,"4-metric",8) == 0)\n} // END function\n')





get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/standalone/standalone_C_code_genpoints.c', '\n// Part P1: Import needed header files\n#include "stdio.h"\n#include "stdlib.h"\n#include "math.h"\n\nconst double xyzmin = -1000.0;\nconst double xyzmax =  1000.0;\n\nvoid write_to_xyz_files(int num_interp_points, char filename_basename[100]) {\n  char filenamex[100],filenamey[100],filenamez[100];\n  snprintf(filenamex,100,"%s-x.dat",filename_basename);\n  snprintf(filenamey,100,"%s-y.dat",filename_basename);\n  snprintf(filenamez,100,"%s-z.dat",filename_basename);\n  FILE *filex = fopen(filenamex,"wb");\n  FILE *filey = fopen(filenamey,"wb");\n  FILE *filez = fopen(filenamez,"wb");\n\n  // Write file headers:\n  fwrite(&num_interp_points, sizeof(int), 1, filex);\n  fwrite(&num_interp_points, sizeof(int), 1, filey);\n  fwrite(&num_interp_points, sizeof(int), 1, filez);\n\n  // Write guts of file:\n  for(int ii=0;ii<num_interp_points;ii++) {\n    double rngx = xyzmin + (xyzmax - xyzmin)*drand48(); // drand48() returns between 0.0 & 1.0\n    double rngy = xyzmin + (xyzmax - xyzmin)*drand48();\n    double rngz = xyzmin + (xyzmax - xyzmin)*drand48();\n    fwrite(&rngx, sizeof(double), 1, filex);\n    fwrite(&rngy, sizeof(double), 1, filey);\n    fwrite(&rngz, sizeof(double), 1, filez);\n  }\n\n  // Write magic number as file footers:\n  int magic_number = -349289480;\n  fwrite(&magic_number, sizeof(int), 1, filex);\n  fwrite(&magic_number, sizeof(int), 1, filey);\n  fwrite(&magic_number, sizeof(int), 1, filez);\n\n  // Close files.\n  fclose(filex);\n  fclose(filey);\n  fclose(filez);\n}\n\nint main(int argc, const char *argv[]) {\n\n  // Step 0a: Read command-line input, error out if nonconformant\n  if(argc != 2 || atoi(argv[1]) < 1) {\n    printf("Error: Expected one command-line argument: ./standalone_C_code_genpoints [num_interp_points],\\n");\n    exit(1);\n  }\n\n  const int num_interp_points = atoi(argv[1]);\n\n  char filename_basename[100];\n  sprintf(filename_basename,"cell_centered_points");\n  write_to_xyz_files(num_interp_points, filename_basename);\n  sprintf(filename_basename,"faceim_points");\n  write_to_xyz_files(num_interp_points, filename_basename);\n  sprintf(filename_basename,"facejm_points");\n  write_to_xyz_files(num_interp_points, filename_basename);\n  sprintf(filename_basename,"facekm_points");\n  write_to_xyz_files(num_interp_points, filename_basename);\n  sprintf(filename_basename,"corner_points");\n  write_to_xyz_files(num_interp_points, filename_basename);\n\n  return 0;\n}\n')





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
        # Write the file header, which includes #define's for GAMMA_SPEED_LIMIT and TINYDOUBLE:
        with open(filename, "w") as file:
            file.write("// Parameters needed to ensure velocity computations are robust:\n")
            file.write("#define GAMMA_SPEED_LIMIT 20\n")
            file.write("#define TINYDOUBLE        1e-100\n\n")
    compute_xx0xx1xx2 = ""
    if "SPHERICAL" in gf_interp_list[which_InterpCounter].gf_description:
        compute_xx0xx1xx2 = """
// ONLY NEEDED/USED IF CONVERTING TO SPHERICAL BASIS:
const double Cartx = x[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] - x_center;
const double Carty = y[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] - y_center;
const double Cartz = z[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] - z_center;

const double xx0 = sqrt(Cartx*Cartx + Carty*Carty + Cartz*Cartz);
const double xx1 = acos(Cartz/xx0);
const double xx2 = atan2(Carty,Cartx);\n
"""
    with open(filename, output_type) as file:
        if( "temperature" in str(interp_expr) or
            "Y_e" in str(interp_expr) or
            "eps" in str(interp_expr) or
            "entropy"in str(interp_expr) ):
            file.write("if(*InterpCounter == "+str(which_InterpCounter)+" && enable_nuc_eos) {\n")
        else:
            file.write("if(*InterpCounter == "+str(which_InterpCounter)+") {\n")
        file.write("    // Interpolating: "+gf_interp_list[which_InterpCounter].gf_description+"\n")
        file.write(lp.loop(["i2","i1","i0"],
                           ["cctk_nghostzones[2]","cctk_nghostzones[1]","cctk_nghostzones[0]"],\
                           ["cctk_lsh[2]-cctk_nghostzones[2]",
                            "cctk_lsh[1]-cctk_nghostzones[1]",
                            "cctk_lsh[0]-cctk_nghostzones[0]"],\
                           ["1","1","1"],\
                           ["#pragma omp parallel for","",""],"   ",
                           compute_xx0xx1xx2+kernel))
        file.write("}\n")
    # If successful, return incremented which_InterpCounter:
    return which_InterpCounter+1





NRPyoutfilename = os.path.join(Ccodesdir,"src","list_of_functions_to_interpolate.h")

which_InterpCounter = 1






gammaDD = ixp.register_gridfunctions_for_single_rank2("AUX","gammaDD", "sym01")
alpha = gri.register_gridfunctions("AUX","alpha")
betaU = ixp.register_gridfunctions_for_single_rank1("AUX","betaU")

beta_offsetU0,beta_offsetU1,beta_offsetU2 = par.Cparameters("REAL","modulenamedoesntmatter",
                                                            ["beta_offsetU0","beta_offsetU1","beta_offsetU2"],
                                                            [0.0,0.0,0.0])
betaU[0] += beta_offsetU0
betaU[1] += beta_offsetU1
betaU[2] += beta_offsetU2

gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01")
betaU_dD   = ixp.declarerank2("betaU_dD","nosym")
alpha_dD   = ixp.declarerank1("alpha_dD")

DIM=3

IGMvU = ixp.register_gridfunctions_for_single_rank1("AUX","IGMvU")
BU    = ixp.register_gridfunctions_for_single_rank1("AUX","BU")

gf_interp_list.append(gf_interp("IGM density primitive"))
rho_b       = gri.register_gridfunctions("AUX","rho_b")
interp_expr = rho_b
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

gf_interp_list.append(gf_interp("IGM pressure primitive"))
P = gri.register_gridfunctions("AUX","P")
interp_expr = P
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)





import GRHD.equations as Ge
Ge.u4U_in_terms_of_vU__rescale_vU_by_applying_speed_limit(alpha, betaU, gammaDD, IGMvU)


u4Uzero = Ge.u4U_ito_vU[0]
gf_interp_list.append(gf_interp("u^0: zero (time) component of 4-velocity"))
interp_expr = u4Uzero
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

Gamma_times_ValenciavU = ixp.zerorank1()
Gamma      = alpha*Ge.u4U_ito_vU[0]
ValenciavU = ixp.zerorank1()
for i in range(DIM):
    ValenciavU[i] = (1/alpha * (Ge.rescaledvU[i] + betaU[i])) # simplify?

for i in range(DIM):
    Gamma_times_ValenciavU[i] = Gamma*ValenciavU[i]





import reference_metric as rfm    # NRPy+: Reference metric support

par.set_parval_from_str("reference_metric::CoordSystem","Spherical")
rfm.reference_metric()

Jac_dUCart_dDrfmUD,Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

Sph_basis_Gamma_times_ValenciavU = ixp.zerorank1()
for i in range(DIM):
    for l in range(DIM):
        Sph_basis_Gamma_times_ValenciavU[i] += Jac_dUrfm_dDCartUD[i][l] * Gamma_times_ValenciavU[l]

for i in range(DIM):
    gf_interp_list.append(gf_interp("*SPHERICAL BASIS* Lorentz factor, times Valencia vU"+str(i)))
    interp_expr = Sph_basis_Gamma_times_ValenciavU[i]
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)




for i in range(DIM):
    gf_interp_list.append(gf_interp("(speed-limited) Valencia 3-velocity vU"+str(i)))
    interp_expr = ValenciavU[i]
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

for i in range(DIM):
    if i==0:
        gf_interp_list.append(gf_interp("Local grid resolution dx=dy=dz"))
        invdx0 = sp.symbols('invdx0', real=True)
        interp_expr = 1/invdx0
    else:
        gf_interp_list.append(gf_interp("(speed-limited) IGM 3-velocity vU"+str(i)+" = u^i divided by u^0"))
        interp_expr = Ge.u4U_ito_vU[i]
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)



for i in range(DIM):
    gf_interp_list.append(gf_interp("IGM magnetic field component B"+str(i)))
    interp_expr = BU[i]
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)





get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/unstagger_A_fields.cc', '\n// Include needed ETK & C library header files:\n#include <stdio.h>\n#include <stdlib.h>\n#include <string.h>\n#include <math.h>\n// Needed for dealing with Cactus/ETK infrastructure\n#include "cctk.h"\n#include "cctk_Arguments.h"\n#include "cctk_Parameters.h"\n\nvoid unstagger_A_fields(CCTK_ARGUMENTS) {\n    DECLARE_CCTK_ARGUMENTS;\n    DECLARE_CCTK_PARAMETERS;\n\n    // Set Ai_unstaggered = Ai and exit the function if A fields are unstaggered already.\n    if(A_fields_are_staggered == 0) {\n#pragma omp parallel for\n        for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {\n            int index=CCTK_GFINDEX3D(cctkGH,i,j,k);\n            Ax_unstaggered[index] = Ax[index];\n            Ay_unstaggered[index] = Ay[index];\n            Az_unstaggered[index] = Az[index];\n        }\n        return;\n    }\n    printf("Unstaggering A fields on grid with dx = %e!\\n",CCTK_DELTA_SPACE(0));\n    // If A fields are staggered (IllinoisGRMHD-style), then unstagger them:\n    // First unstagger A_x, which is defined at (i, j+1/2, k+1/2). Unstaggering\n    //   is as simple as A_x(i,j,k) = 1/4 * (A_x(i,j-1/2,k-1/2)+A_x(i,j-1/2,k+1/2)+A_x(i,j+1/2,k-1/2)+A_x(i,j+1/2,k+1/2))\n#pragma omp parallel for\n    for(int k=1;k<cctk_lsh[2];k++) for(int j=1;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {\n        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);\n        Ax_unstaggered[index] = 0.25*(Ax[CCTK_GFINDEX3D(cctkGH,i,j,k)]     + Ax[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] +\n                                      Ax[CCTK_GFINDEX3D(cctkGH,i,j-1,k-1)] + Ax[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);\n    }\n#pragma omp parallel for\n    for(int k=1;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=1;i<cctk_lsh[0];i++) {\n        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);\n        Ay_unstaggered[index] = 0.25*(Ay[CCTK_GFINDEX3D(cctkGH,i,j,k)]     + Ay[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] +\n                                      Ay[CCTK_GFINDEX3D(cctkGH,i-1,j,k-1)] + Ay[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);\n    }\n#pragma omp parallel for\n    for(int k=0;k<cctk_lsh[2];k++) for(int j=1;j<cctk_lsh[1];j++) for(int i=1;i<cctk_lsh[0];i++) {\n        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);\n        Az_unstaggered[index] = 0.25*(Az[CCTK_GFINDEX3D(cctkGH,i,j,k)]     + Az[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] +\n                                      Az[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k)] + Az[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);\n    }\n}\n')





Ax_unstaggered = gri.register_gridfunctions("AUX","Ax_unstaggered")
Ay_unstaggered = gri.register_gridfunctions("AUX","Ay_unstaggered")
Az_unstaggered = gri.register_gridfunctions("AUX","Az_unstaggered")

gf_interp_list.append(gf_interp("Unstaggered vector potential component Ax"))
interp_expr = Ax_unstaggered
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

gf_interp_list.append(gf_interp("Unstaggered vector potential component Ay"))
interp_expr = Ay_unstaggered
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

gf_interp_list.append(gf_interp("Unstaggered vector potential component Az"))
interp_expr = Az_unstaggered
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

for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            # Recall that betaD[i] = gammaDD[i][j]*betaU[j] (Eq. 2.121 in B&S)
            betaDdD[i][k] += gammaDD_dD[i][j][k]*betaU[j] + gammaDD[i][j]*betaU_dD[j][k]

g4DDdD = ixp.zerorank3(DIM=4)
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











































Weylgfs = ["Psi0r","Psi0i","Psi1r","Psi1i","Psi2r","Psi2i","Psi3r","Psi3i","Psi4r","Psi4i",
                                    "curvIr","curvIi","curvJr","curvJi"]
Psi0r,Psi0i,Psi1r,Psi1i,Psi2r,Psi2i,Psi3r,Psi3i,Psi4r,Psi4i,curvIr,curvIi,curvJr,curvJi = \
    gri.register_gridfunctions("AUX",Weylgfs)
count = 0
for gf in [Psi0r,Psi0i,Psi1r,Psi1i,Psi2r,Psi2i,Psi3r,Psi3i,Psi4r,Psi4i,curvIr,curvIi,curvJr,curvJi]:
    gf_interp_list.append(gf_interp("4-Weyl scalar or invariant "+Weylgfs[count]))
    interp_expr = gf
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
    count = count + 1





gf_interp_list.append(gf_interp("Constraints - Hamiltonian"))
H = gri.register_gridfunctions("AUX","InterpH")
interp_expr = H
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

MU = ixp.register_gridfunctions_for_single_rank1("AUX","InterpMU")


for i in range(DIM):
    gf_interp_list.append(gf_interp("Constraints - Momentum M^"+chr(ord('x')+i)))
    interp_expr = MU[i]
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)





gf_interp_list.append(gf_interp("Temperature primitive"))
temperature = gri.register_gridfunctions("AUX","temperature")
interp_expr = temperature
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

gf_interp_list.append(gf_interp("Electron fraction primitive"))
Ye = gri.register_gridfunctions("AUX","Y_e")
interp_expr = Ye
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

gf_interp_list.append(gf_interp("Specific internal energy primitive"))
eps = gri.register_gridfunctions("AUX","eps")
interp_expr = eps
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

gf_interp_list.append(gf_interp("Entropy primitive"))
entropy = gri.register_gridfunctions("AUX","entropy")
interp_expr = entropy
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)





get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/construct_function_to_interpolate__store_to_interped_gf.cc', '#include <stdio.h>\n#include <stdlib.h>\n#include "cctk.h"\n#include "cctk_Arguments.h"\n#include "cctk_Parameters.h"\n\n// Set the gridfunction interped_gf, according to the interpolation counter variable interp_counter.\n//    For example, we might interpolate "IllinoisGRMHD::rho_b" if interp_counter==0. The following\n//    function takes care of these\nvoid list_of_functions_to_interpolate(cGH *restrict cctkGH,\n                                      const CCTK_INT *restrict cctk_lsh,\n                                      const CCTK_INT *restrict cctk_nghostzones,\n                                      const CCTK_REAL *restrict x,const CCTK_REAL *restrict y,const CCTK_REAL *restrict z,\n                                      const CCTK_REAL invdx0,const CCTK_REAL invdx1,const CCTK_REAL invdx2,\n                                      const CCTK_INT *restrict InterpCounter,\n                                      const CCTK_REAL *restrict rho_bGF,\n                                      const CCTK_REAL *restrict PGF,\n                                      const CCTK_REAL *restrict temperatureGF,\n                                      const CCTK_REAL *restrict Y_eGF,\n                                      const CCTK_REAL *restrict epsGF,\n                                      const CCTK_REAL *restrict entropyGF,\n                                      const CCTK_REAL *restrict IGMvU0GF,const CCTK_REAL *restrict IGMvU1GF,const CCTK_REAL *restrict IGMvU2GF,\n                                      const CCTK_REAL *restrict BU0GF,const CCTK_REAL *restrict BU1GF,const CCTK_REAL *restrict BU2GF,\n                                      const CCTK_REAL *restrict gammaDD00GF,const CCTK_REAL *restrict gammaDD01GF,const CCTK_REAL *restrict gammaDD02GF,\n                                      const CCTK_REAL *restrict gammaDD11GF,const CCTK_REAL *restrict gammaDD12GF,const CCTK_REAL *restrict gammaDD22GF,\n                                      const CCTK_REAL *restrict betaU0GF,const CCTK_REAL *restrict betaU1GF,const CCTK_REAL *restrict betaU2GF,\n                                      const CCTK_REAL *restrict alphaGF,\n                                      CCTK_REAL *restrict interped_gfGF,\n                                      CCTK_REAL *restrict Ax_unstaggeredGF,CCTK_REAL *restrict Ay_unstaggeredGF,CCTK_REAL *restrict Az_unstaggeredGF,\n                                      const CCTK_REAL *restrict Psi0rGF,const CCTK_REAL *restrict Psi0iGF,\n                                      const CCTK_REAL *restrict Psi1rGF,const CCTK_REAL *restrict Psi1iGF,\n                                      const CCTK_REAL *restrict Psi2rGF,const CCTK_REAL *restrict Psi2iGF,\n                                      const CCTK_REAL *restrict Psi3rGF,const CCTK_REAL *restrict Psi3iGF,\n                                      const CCTK_REAL *restrict Psi4rGF,const CCTK_REAL *restrict Psi4iGF,\n                                      const CCTK_REAL *restrict curvIrGF,const CCTK_REAL *restrict curvIiGF,\n                                      const CCTK_REAL *restrict curvJrGF,const CCTK_REAL *restrict curvJiGF,\n                                      const CCTK_REAL *restrict InterpHGF,\n                                      const CCTK_REAL *restrict InterpMU0GF,\n                                      const CCTK_REAL *restrict InterpMU1GF,\n                                      const CCTK_REAL *restrict InterpMU2GF) {\n  DECLARE_CCTK_PARAMETERS;\n#include "list_of_functions_to_interpolate.h"\n}\n\nvoid construct_function_to_interpolate__store_to_interped_gf(CCTK_ARGUMENTS) {\n  DECLARE_CCTK_ARGUMENTS;\n  DECLARE_CCTK_PARAMETERS;\n  printf("Called construct_function_to_interpolate__store_to_interped_gf() on grid with dx = %e!\\n",CCTK_DELTA_SPACE(0));\n  const CCTK_REAL invdx0 = 1.0 / CCTK_DELTA_SPACE(0);\n  const CCTK_REAL invdx1 = 1.0 / CCTK_DELTA_SPACE(1);\n  const CCTK_REAL invdx2 = 1.0 / CCTK_DELTA_SPACE(2);\n\n  // .------------------------------------.\n  // | Constraint equations gridfunctions |\n  // .------------------------------------.\n  // We get the pointer using the CCTK_VarDataPtr() function because\n  // if the gridfunction does not exist in this version of IGM then\n  // we will get compilation errors.\n  int timelevel = 0;\n\n  // Hamiltonian constraint gridfunction\n  CCTK_REAL *InterpHGF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,HVarString));\n  if( !InterpHGF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "Hamiltonian constraint - Couldn\'t get data pointer of input array variable \'%s\'", HVarString);\n\n  // Momentum constraint gridfunction (x-direction)\n  CCTK_REAL *InterpMU0GF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,M0VarString));\n  if( !InterpMU0GF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "Momentum constraint (x-direction) - Couldn\'t get data pointer of input array variable \'%s\'", M0VarString);\n\n  // Momentum constraint gridfunction (y-direction)\n  CCTK_REAL *InterpMU1GF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,M1VarString));\n  if( !InterpMU1GF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "Momentum constraint (y-direction) - Couldn\'t get data pointer of input array variable \'%s\'", M1VarString);\n\n  // Momentum constraint gridfunction (z-direction)\n  CCTK_REAL *InterpMU2GF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,M2VarString));\n  if( !InterpMU2GF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "Momentum constraint (z-direction) - Couldn\'t get data pointer of input array variable \'%s\'", M2VarString);\n\n  // .-------------------------------------------------------.\n  // | Additional hydro quantities for nuclear/tabulated EOS |\n  // .-------------------------------------------------------.\n  CCTK_REAL *epsGF, *temperatureGF, *Y_eGF, *entropyGF;\n  if( enable_nuc_eos ) {\n    // Temperature gridfunction\n    temperatureGF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,temperatureVarString));\n    if( !temperatureGF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "temperature - Couldn\'t get data pointer of input array variable \'%s\'", temperatureVarString);\n\n    // Electron fraction gridfunction\n    Y_eGF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,Y_eVarString));\n    if( !Y_eGF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "Y_e - Couldn\'t get data pointer of input array variable \'%s\'", Y_eVarString);\n\n    // Specific internal energy gridfunction\n    epsGF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,epsVarString));\n    if( !epsGF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "eps - Couldn\'t get data pointer of input array variable \'%s\'", epsVarString);\n\n    // Entropy gridfunction\n    entropyGF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,entropyVarString));\n    if( !entropyGF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "entropy - Couldn\'t get data pointer of input array variable \'%s\'", entropyVarString);\n  }\n  else {\n    // If we are not using nuclear EOS, then set the gridfunction pointers to NULL.\n    temperatureGF = NULL;\n    Y_eGF         = NULL;\n    epsGF         = NULL;\n    entropyGF     = NULL;\n  }\n\n  list_of_functions_to_interpolate(cctkGH,cctk_lsh,cctk_nghostzones,\n                                   x,y,z,\n                                   invdx0,invdx1,invdx2,\n                                   InterpCounter,\n                                   rho_b,P,temperatureGF,Y_eGF,epsGF,entropyGF,\n                                   vx,vy,vz,\n                                   Bx,By,Bz,\n                                   gxx,gxy,gxz,gyy,gyz,gzz,\n                                   betax,betay,betaz,alp, interped_gf,\n                                   Ax_unstaggered,Ay_unstaggered,Az_unstaggered,\n                                   Psi0r,Psi0i,Psi1r,Psi1i,Psi2r,Psi2i,Psi3r,Psi3i,Psi4r,Psi4i,\n                                   curvIr,curvIi,curvJr,curvJi,\n                                   InterpHGF,InterpMU0GF,InterpMU1GF,InterpMU2GF);\n\n  // interped_gf will be interpolated across AMR boundaries, meaning that\n  //    it must be prolongated. Only gridfunctions with 3 timelevels stored\n  //    may be prolongated (provided time_interpolation_order is set to the\n  //    usual value of 2). We should only call this interpolation routine\n  //    at iterations in which all gridfunctions are on the same timelevel\n  //    (usually a power of 2), which will ensure that the following\n  //    "filling of the timelevels" is completely correct.\n#pragma omp parallel for\n  for(int i=0;i<cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];i++) {\n    interped_gf_p[i]   = interped_gf[i];\n    interped_gf_p_p[i] = interped_gf[i];\n  }\n}\n')





with open(os.path.join(Ccodesdir,"src","get_gf_name.h"), "w") as file:
    file.write("void get_gf_name(const int InterpCounter,char gf_name[100]) {\n")
    for i in range(1,which_InterpCounter):
        file.write("    if(InterpCounter=="+str(i)+") { snprintf(gf_name,100,\""+gf_interp_list[i].gf_description+"\"); return; }\n")
    file.write("    printf(\"Error. InterpCounter = %d unsupported. I should not be here.\\n\",InterpCounter); exit(1);\n")
    file.write("}\n")





with open(os.path.join(Ccodesdir,"src","define_NumInterpFunctions.h"), "w") as file:
    file.write("#define NumInterpFunctions ("+str(which_InterpCounter)+"-4*(1-enable_nuc_eos))\n")




get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/interp_counter.cc', '#include <assert.h>\n#include <stdio.h>\n#include <stdlib.h>\n#include <string.h>\n#include <math.h>\n#include <ctype.h>\n#include "cctk.h"\n#include "cctk_Arguments.h"\n#include "cctk_Parameters.h"\n\n#include "define_NumInterpFunctions.h"\n\nvoid ArbGrid_InitializeInterpCounterToZero(CCTK_ARGUMENTS)\n{\n  DECLARE_CCTK_ARGUMENTS;\n  DECLARE_CCTK_PARAMETERS;\n  *InterpCounter = 0;\n\n  if(verbose==2) printf("interp_arbgrid_MO_ETK: Just set InterpCounter to %d\\n",*InterpCounter);\n}\n\nvoid ArbGrid_InitializeInterpCounter(CCTK_ARGUMENTS)\n{\n  DECLARE_CCTK_ARGUMENTS;\n  DECLARE_CCTK_PARAMETERS;\n\n  if(cctk_iteration == interp_out_iteration) {\n    *InterpCounter = 1;\n    if(verbose==2) printf("interp_arbgrid_MO_ETK: Just set InterpCounter to %d ; ready to start looping over interpolated gridfunctions!\\n",\n                          *InterpCounter);\n  }\n}\n\n// This function increments InterpCounter if we are at the interp_out_iteration until\n// it hits NumInterpFunctions. At this iteration, InterpCounter is set to zero, which\n// exits the loop.\nvoid ArbGrid_IncrementInterpCounter(CCTK_ARGUMENTS)\n{\n  DECLARE_CCTK_ARGUMENTS;\n  DECLARE_CCTK_PARAMETERS;\n\n  if(*InterpCounter == NumInterpFunctions-1) {\n    *InterpCounter = 0;\n    if(verbose==2) printf("interp_arbgrid_MO_ETK: Finished! Just zeroed InterpCounter.\\n");\n  } else {\n    (*InterpCounter)++;\n    if(verbose==2) printf("interp_arbgrid_MO_ETK: Just incremented InterpCounter to %d of %d\\n",*InterpCounter,NumInterpFunctions-1);\n  }\n}\n')








































get_ipython().run_cell_magic('writefile', '$Ccodesdir/src/make.code.defn', '# Main make.code.defn file for thorn interp_arbgrid_MO_ETK\n\n# Source files in this directory\nSRCS =  main_function.cc unstagger_A_fields.cc interp_counter.cc \\\n        construct_function_to_interpolate__store_to_interped_gf.cc # FM_validation.cc # <- For FishboneMoncriefID validation\n')





get_ipython().run_cell_magic('writefile', '$Ccodesdir/interface.ccl', '\n# With "implements", we give our thorn its unique name.\nimplements: interp_arbgrid_MO_ETK\n\n# By "inheriting" other thorns, we tell the Toolkit that we\n#   will rely on variables/function that exist within those\n#   functions.\ninherits:   admbase IllinoisGRMHD Grid\ninherits:   WeylScal4 # Needed for Weyl scalars psi4, psi3, psi..., and Weyl invariants I & J.\n# For FM ID comparisons:\n# inherits:   FishboneMoncriefID\n\n# Tell the Toolkit that we want "interped_gf" and "InterpCounter"\n#    and invariants to NOT be visible to other thorns, by using\n#    the keyword "private". Note that declaring these\n#    gridfunctions here *does not* allocate memory for them;\n#    that is done by the schedule.ccl file.\nprivate:\nCCTK_REAL interpolation_gf type=GF timelevels=3 tags=\'Checkpoint="no"\'\n{\n  interped_gf\n} "Gridfunction containing output from interpolation."\n\nCCTK_REAL unstaggered_A_fields type=GF timelevels=3 tags=\'Checkpoint="no"\'\n{\n  Ax_unstaggered,Ay_unstaggered,Az_unstaggered\n} "Unstaggered A-field components."\n\n\nint InterpCounterVar type = SCALAR tags=\'checkpoint="no"\'\n{\n  InterpCounter\n} "Counter that keeps track of which function we are interpolating."\n')





get_ipython().run_cell_magic('writefile', '$Ccodesdir/param.ccl', '\n# Output the interpolated data to the IO::out_dir directory:\nshares: IO\nUSES STRING out_dir\n\n# These parameters are used to output Hybrid EOS information\nshares: IllinoisGRMHD\nUSES CCTK_INT neos\nUSES CCTK_REAL Gamma_th\n\n# These parameters are used to output nuclear (tabulated) EOS information\nshares: EOS_Omni\nUSES CCTK_STRING nuceos_table_name\n\n# For FM ID comparisons:\n# shares: FishboneMoncriefID\n# USES KEYWORD M\n# USES KEYWORD a\n# USES KEYWORD r_in\n# USES KEYWORD r_at_max_density\n\nrestricted:\n\n########################################\n# BASIC THORN STEERING PARAMETERS\nCCTK_INT interp_out_iteration "Which iteration to interpolate to destination grids?" STEERABLE=ALWAYS\n{\n  0:* :: ""\n} 960000\n\n## Interpolator information\nCCTK_INT verbose "Set verbosity level: 1=useful info; 2=moderately annoying (though useful for debugging)" STEERABLE=ALWAYS\n{\n  0:2 :: "0 = no output; 1=useful info; 2=moderately annoying (though useful for debugging)"\n} 2\n\nCCTK_INT A_fields_are_staggered "Are A fields staggered? 1 = yes; 0 = no. Default to yes." STEERABLE=ALWAYS\n{\n  0:1 :: ""\n} 1\n\n##########\n# Cartesian position of center of output grid (usually center of BH).\nCCTK_REAL x_center "x-position of center." STEERABLE=ALWAYS\n{\n  *:* :: ""\n} 0.0\n\nCCTK_REAL y_center "y-position of center." STEERABLE=ALWAYS\n{\n  *:* :: ""\n} 0.0\n\nCCTK_REAL z_center "z-position of center." STEERABLE=ALWAYS\n{\n  *:* :: ""\n} 0.0\n\n##########\n# Shift offset:\nCCTK_REAL beta_offsetU0 "Offset to betax, to account for coordinate drift in x direction." STEERABLE=ALWAYS\n{\n  *:* :: ""\n} 0.0\n\nCCTK_REAL beta_offsetU1 "Offset to betay, to account for coordinate drift in y direction." STEERABLE=ALWAYS\n{\n  *:* :: ""\n} 0.0\n\nCCTK_REAL beta_offsetU2 "Offset to betaz, to account for coordinate drift in z direction." STEERABLE=ALWAYS\n{\n  *:* :: ""\n} 0.0\n\nCCTK_REAL out_of_bounds_interp_xyz "Do not interpolate points with fabs(xyz) > out_of_bounds_interp_xyz, where xyz are centered at xyz_center (usually center of BH). Fill dataset with NaN instead." STEERABLE=ALWAYS\n{\n  0:* :: "Any positive number"\n} 1E100\n\n###########\n# Nuclear EOS options\nCCTK_INT enable_nuc_eos "Use nuclear (tabulated/microphysics) equation of state? 1 = yes; 0 = no. Default to no." STEERABLE=ALWAYS\n{\n  0:1 :: ""\n} 0\n\nCCTK_STRING temperatureVarString "Temperature GF name. Defaults to IllinoisGRMHD\'s." STEERABLE=ALWAYS\n{\n  "IllinoisGRMHD::igm_temperature" :: "IllinoisGRMHD\'s temperature gridfunction name"\n  "HydroBase::temperature"         :: "HydroBase\'s temperature gridfunction name"\n  ".+"                             :: "Or use you can use your own thorn\'s temperature gridfunction name"\n} "IllinoisGRMHD::igm_temperature"\n\nCCTK_STRING Y_eVarString "Electron fraction GF name. Defaults to IllinoisGRMHD\'s." STEERABLE=ALWAYS\n{\n  "IllinoisGRMHD::igm_Ye" :: "IllinoisGRMHD\'s electron fraction gridfunction name"\n  "HydroBase::Y_e"        :: "HydroBase\'s  electron fraction gridfunction name"\n  ".+"                    :: "Or use you can use your own thorn\'s  electron fraction gridfunction name"\n} "IllinoisGRMHD::igm_Ye"\n\nCCTK_STRING epsVarString "Specific internal energy GF name. Defaults to IllinoisGRMHD\'s." STEERABLE=ALWAYS\n{\n  "IllinoisGRMHD::igm_eps" :: "IllinoisGRMHD\'s specific internal energy gridfunction name"\n  "HydroBase::eps"         :: "HydroBase\'s specific internal energy gridfunction name"\n  ".+"                     :: "Or use you can use your own thorn\'s specific internal energy gridfunction name"\n} "IllinoisGRMHD::igm_eps"\n\nCCTK_STRING entropyVarString "Entropy GF name. Defaults to IllinoisGRMHD\'s." STEERABLE=ALWAYS\n{\n  "IllinoisGRMHD::igm_entropy" :: "IllinoisGRMHD\'s entropy gridfunction name"\n  "HydroBase::entropy"         :: "HydroBase\'s entropy gridfunction name"\n  ".+"                         :: "Or use you can use your own thorn\'s entropy gridfunction name"\n} "IllinoisGRMHD::igm_entropy"\n\nCCTK_STRING HVarString "Hamiltonian constraint GF name. Defaults to Baika\'s." STEERABLE=ALWAYS\n{\n  "Baikal::HGF" :: "Baikal\'s Hamiltonian constraint gridfunction name"\n  "ML_BSSN::H"  :: "ML_BSSN\'s Hamiltonian constraint gridfunction name"\n  ".+"                         :: "Or use you can use your own thorn\'s entropy gridfunction name"\n} "Baikal::HGF"\n\nCCTK_STRING M0VarString "Momentum constraint (x-direction) Defaults to Baikal\'s." STEERABLE=ALWAYS\n{\n  "Baikal::MU0GF" :: "Baikal\'s Momentum constraint (x-direction) gridfunction name"\n  "ML_BSSN::M1"   :: "ML_BSSN\'s Momentum constraint (x-direction) gridfunction name"\n  ".+"                         :: "Or use you can use your own thorn\'s entropy gridfunction name"\n} "Baikal::MU0GF"\n\nCCTK_STRING M1VarString "Momentum constraint (y-direction) GF name. Defaults to Baikal\'s." STEERABLE=ALWAYS\n{\n  "Baikal::MU1GF" :: "Baikal\'s Momentum constraint (y-direction) gridfunction name"\n  "ML_BSSN::M2"   :: "ML_BSSN\'s Momentum constraint (y-direction) gridfunction name"\n  ".+"                         :: "Or use you can use your own thorn\'s entropy gridfunction name"\n} "Baikal::MU1GF"\n\nCCTK_STRING M2VarString "Momentum constraint (z-direction) GF name. Defaults to Baikal\'s." STEERABLE=ALWAYS\n{\n  "Baikal::MU2GF" :: "Baikal\'s Momentum constraint (z-direction) gridfunction name"\n  "ML_BSSN::M3"   :: "ML_BSSN\'s Momentum constraint (z-direction) gridfunction name"\n  ".+"                         :: "Or use you can use your own thorn\'s entropy gridfunction name"\n} "Baikal::MU2GF"\n')





get_ipython().run_cell_magic('writefile', '$Ccodesdir/schedule.ccl', '\nSTORAGE: interpolation_gf[3]\nSTORAGE: unstaggered_A_fields[3]\nSTORAGE: InterpCounterVar\n# STORAGE: interp_pointcoords_and_output_arrays\n\n#############################\nSCHEDULE ArbGrid_InitializeInterpCounterToZero AT CCTK_INITIAL\n{\n  LANG: C\n  OPTIONS: GLOBAL\n} "Initialize InterpCounter variable to zero"\n\nSCHEDULE ArbGrid_InitializeInterpCounterToZero AT CCTK_POST_RECOVER_VARIABLES\n{\n  LANG: C\n  OPTIONS: GLOBAL\n} "Initialize InterpCounter variable to zero"\n\nSCHEDULE ArbGrid_InitializeInterpCounter before ArbGrid_InterpGroup AT CCTK_ANALYSIS\n{\n  LANG: C\n  OPTIONS: GLOBAL\n} "Initialize InterpCounter variable"\n##################\n\nSCHEDULE GROUP ArbGrid_InterpGroup AT CCTK_ANALYSIS BEFORE CarpetLib_printtimestats BEFORE CarpetLib_printmemstats AFTER Convert_to_HydroBase WHILE interp_arbgrid_MO_ETK::InterpCounter\n{\n} "Perform all interpolations. This group is only actually scheduled at cctk_iteration==interp_out_iteration."\n\nSCHEDULE unstagger_A_fields in ArbGrid_InterpGroup before construct_function_to_interpolate__store_to_interped_gf\n{\n  STORAGE: unstaggered_A_fields[3]\n  OPTIONS: GLOBAL,LOOP-LOCAL\n  SYNC:    unstaggered_A_fields\n  LANG: C\n} "Unstagger A fields."\n\nSCHEDULE construct_function_to_interpolate__store_to_interped_gf in ArbGrid_InterpGroup before DoSum\n{\n  STORAGE: interpolation_gf[3],InterpCounterVar\n  OPTIONS: GLOBAL,LOOP-LOCAL\n  SYNC:    interpolation_gf\n  LANG: C\n} "Construct the function to interpolate"\n\nSCHEDULE Interpolate_to_dest_grid_main_function in ArbGrid_InterpGroup after construct_function_to_interpolate__store_to_interped_gf\n{\n  OPTIONS: GLOBAL\n  LANG: C\n} "Perform interpolation and output result to file."\n\n# For FishboneMoncriefID validation only.\n# SCHEDULE FM_validation in ArbGrid_InterpGroup after Interpolate_to_dest_grid_main_function\n# {\n#   OPTIONS: GLOBAL\n#   LANG: C\n# } "Perform interpolation and output result to 2D ASCII file."\n\n#######\nSCHEDULE ArbGrid_IncrementInterpCounter in ArbGrid_InterpGroup after Interpolate_to_dest_grid_main_function\n{\n  LANG: C\n  OPTIONS: GLOBAL\n} "Increment InterpCounter variable, or set to zero once loop is complete."\n##################\n')






import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ETK_thorn-Interpolation_to_Arbitrary_Grids_multi_order")

