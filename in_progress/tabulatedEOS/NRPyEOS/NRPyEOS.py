# (c) 2022, Leo Werneck
#
# NRPyEOS.py
#
# This file contains functions to generate the C code
# for interpolating the hydrodynamic quantities in the
# equation of state tables.

# Step 1: Initialize NRPy+/Python modules
import os, sys                           # Standard Python modules for multiplatform OS-level functions, benchmarking
from collections import namedtuple       # Standard Python: Enable namedtuple data type
sys.path.append(os.path.join("..","..")) # Add NRPy+'s base directory to Python's path
import outputC as outC                   # NRPy+: Core C code output module

# Step 2: Identify EOS table quantities
# Step 2.a: Create the EOS named tuple
eos_tuple = namedtuple("eos_tuple","n desc key var")

# Step 2.b: Create tuples for all EOS quantities
P      = eos_tuple( 0, "Pressure"                          ,"NRPyEOS_press_key"  ,"P"     )
eps    = eos_tuple( 1, "Energy"                            ,"NRPyEOS_eps_key"    ,"eps"   )
S      = eos_tuple( 2, "Entropy"                           ,"NRPyEOS_entropy_key","S"     )
munu   = eos_tuple( 3, "Neutrino chemical potential"       ,"NRPyEOS_munu_key"   ,"munu"  )
cs2    = eos_tuple( 4, "Soundspeed"                        ,"NRPyEOS_cs2_key"    ,"cs2"   )
depsdT = eos_tuple( 5, "Derivative of eps w.r.t T"         ,"NRPyEOS_depsdT_key" ,"depsdT")
dPdrho = eos_tuple( 6, "Derivative of P w.r.t rho"         ,"NRPyEOS_dPdrho_key" ,"dPdrho")
dPdeps = eos_tuple( 7, "Derivative of P w.r.t eps"         ,"NRPyEOS_dPdeps_key" ,"dPdT"  )
muhat  = eos_tuple( 8, "mu_n - mu_p"                       ,"NRPyEOS_muhat_key"  ,"muhat" )
mu_e   = eos_tuple( 9, "Electron chemical potential"       ,"NRPyEOS_mu_e_key"   ,"mu_e"  )
mu_p   = eos_tuple(10, "Proton chemical potential"         ,"NRPyEOS_mu_p_key"   ,"mu_p"  )
mu_n   = eos_tuple(11, "Neutron chemical potential"        ,"NRPyEOS_mu_n_key"   ,"mu_n"  )
X_a    = eos_tuple(12, "Alpha particle mass fraction"      ,"NRPyEOS_X_a_key"    ,"X_a"   )
X_h    = eos_tuple(13, "Heavy nuclei mass fraction"        ,"NRPyEOS_X_h_key"    ,"X_h"   )
X_n    = eos_tuple(14, "Neutron mass fraction"             ,"NRPyEOS_X_n_key"    ,"X_n"   )
X_p    = eos_tuple(15, "Proton mass fraction"              ,"NRPyEOS_X_p_key"    ,"X_p"   )
Abar   = eos_tuple(16, "Avg. mass number of heavy nuclei"  ,"NRPyEOS_Abar_key"   ,"Abar"  )
Zbar   = eos_tuple(17, "Avg. charge number of heavy nuclei","NRPyEOS_Zbar_key"   ,"Zbar"  )
Gamma  = eos_tuple(18, "Adiabatic index"                   ,"NRPyEOS_Gamma_key"  ,"Gamma" )

# Step 3: Helper functions
# Step 3.a: Set function name
def func_name(eos_params,auxvar_name):
    N_params = len(eos_params)
    name = "NRPyEOS"
    for i in range(N_params):
        if auxvar_name == "T" and i == len(eos_params)-1 and N_params > 1:
            name += "_and"
        name += "_"+eos_params[i].var.replace("_","")
    if auxvar_name != "T":
        name += "_and_T"
    name += "_from_rho_Ye_"+auxvar_name
    return name

# Step 3.b: Determine identation of function parameters
def param_indentation(c_type,name):
    indent = "  " # Parenthesis and space between type and name
    for i in range(len(c_type)+len(name)):
        indent += " "
    return indent

# Step 3.c: Set function parameters
def func_params(c_type,name,auxvar_name,eos_params,unknown_T=False):
    indent   = param_indentation(c_type,name)
    params   = "const NRPyEOS_params *restrict eos_params,\n"
    params  += indent+"const double rho,\n"
    params  += indent+"const double Y_e,\n"
    params  += indent+"const double "+auxvar_name+",\n"
    N_params = len(eos_params)
    for i in range(N_params):
        if i == N_params-1:
            params += indent+"double *restrict "+eos_params[i].var
        else:
            params += indent+"double *restrict "+eos_params[i].var+",\n"
    if unknown_T:
        params += ",\n"+indent+"double *restrict T"
    return params

# Step 3.d: Set function body
def func_body(name,eos_params,auxvar):
    N_params   = len(eos_params)
    indent     = "  "
    body       = ""
    body      += indent+"// Step 1: Set EOS table keys\n"
    body      += indent+"const int keys["+str(N_params)+"] = {"
    for i in range(N_params):
        if i == N_params-1:
            body += eos_params[i].key
        else:
            body += eos_params[i].key+","
    body      += "};\n\n"
    body      += indent+"// Step 2: Declare EOS error report struct\n"
    body      += indent+"NRPyEOS_error_report report;\n\n"
    body      += indent+"// Step 3: Declare output array\n"
    body      += indent+"double outvars["+str(N_params)+"];\n\n"
    body      += indent+"// Step 4: Perform the interpolation\n"
    if auxvar == "T":
        body  += indent+"NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, "+str(N_params)+",rho,Y_e,T, keys,outvars, &report );\n\n"
    else:
        body  += indent+"const double root_finding_precision = 1e-10;"
        body  += """
  NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( eos_params, """+str(N_params)+""",root_finding_precision,
                                                               rho,Y_e,"""+auxvar.var+","+auxvar.key+""", keys,outvars, T, &report );\n\n"""
    body      += indent+"// Step 5: Check for errors"
    body      += r"""
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside """+name+""". Error message: %s (key = %d)",report.message,report.error_key);
  }\n\n"""
    body      += indent+"// Step 6: Update output variables\n"
    for i in range(N_params):
        body  += indent+"*"+eos_params[i].var+" = outvars["+str(i)+"];\n"
    return body

# Step 3.e: Functions for which the temperature is known
def Cfunc_known_T(eos_params_in):
    eos_params = sorted(eos_params_in)
    includes   = ["NRPy_basic_defines.h","NRPy_function_prototypes.h"]
    desc       = "(c) 2022 Leo Werneck"
    c_type     = "void"
    name       = func_name(eos_params,"T")
    params     = func_params(c_type,name,"T",eos_params)
    body       = func_body(name,eos_params,"T")
    outC.add_to_Cfunction_dict(includes=includes,desc=desc,c_type=c_type,name=name,
                               params=params,body=body,enableCparameters=False)

# Step 3.f: Functions for which the temperature is unknown
def Cfunc_unknown_T(auxvar,eos_params_in):
    eos_params = sorted(eos_params_in)
    includes   = ["NRPy_basic_defines.h","NRPy_function_prototypes.h"]
    desc       = "(c) 2022 Leo Werneck"
    c_type     = "void"
    name       = func_name(eos_params,auxvar.var)
    params     = func_params(c_type,name,auxvar.var,eos_params,unknown_T=True)
    body       = func_body(name,eos_params,auxvar)
    outC.add_to_Cfunction_dict(includes=includes,desc=desc,c_type=c_type,name=name,
                               params=params,body=body,enableCparameters=False)


# Step 4: Supplemental dictionary for NRPy_basic_defines.h
def add_NRPyEOS_header_to_supplementary_dict(supplementary_dict):
    supplementary_dict["NRPyEOS"] = r"""
#include <stdbool.h>
#include <hdf5.h>
#define H5_USE_16_API 1

// EOS struct
typedef struct _NRPyEOS_params_ {

  // Number of points
  int nrho;
  int ntemp;
  int nye;

  // Table arrays
  double *restrict alltables;
  double *restrict epstable;
  double *restrict logrho;
  double *restrict logtemp;
  double *restrict yes;

  // Minimum and maximum values of
  // rho, Ye, and T
  double eos_rhomax , eos_rhomin;
  double eos_tempmin, eos_tempmax;
  double eos_yemin  , eos_yemax;

  // Auxiliary variables
  double energy_shift;
  double temp0, temp1;
  double dlintemp, dlintempi;
  double drholintempi;
  double dlintempyei;
  double drholintempyei;
  double dtemp, dtempi;
  double drho, drhoi;
  double dye, dyei;
  double drhotempi;
  double drhoyei;
  double dtempyei;
  double drhotempyei;

} NRPyEOS_params;

// Table keys
#define NRPyEOS_press_key    0
#define NRPyEOS_eps_key      1
#define NRPyEOS_entropy_key  2
#define NRPyEOS_munu_key     3
#define NRPyEOS_cs2_key      4
#define NRPyEOS_depsdT_key   5
#define NRPyEOS_dPdrho_key   6
#define NRPyEOS_dPdeps_key   7
#define NRPyEOS_muhat_key    8
#define NRPyEOS_mu_e_key     9
#define NRPyEOS_mu_p_key    10
#define NRPyEOS_mu_n_key    11
#define NRPyEOS_X_a_key     12
#define NRPyEOS_X_h_key     13
#define NRPyEOS_X_n_key     14
#define NRPyEOS_X_p_key     15
#define NRPyEOS_Abar_key    16
#define NRPyEOS_Zbar_key    17
#define NRPyEOS_Gamma_key   18
#define NRPyEOS_ntablekeys  19

// Unit conversion
#define LENGTHGF 6.77269222552442e-06
#define TIMEGF 2.03040204956746e05
#define RHOGF 1.61887093132742e-18
#define PRESSGF 1.80123683248503e-39
#define EPSGF 1.11265005605362e-21
#define INVRHOGF 6.17714470405638e17
#define INVEPSGF 8.98755178736818e20
#define INVPRESSGF 5.55174079257738e38

// Name of the variables. This is only used to print
// information about the keys during startup
static const char table_var_names[NRPyEOS_ntablekeys][10] = {
  "logpress","logenergy","entropy","munu","cs2","dedt",
  "dpdrhoe", "dpderho", "muhat", "mu_e", "mu_p", "mu_n",
  "Xa","Xh","Xn","Xp","Abar","Zbar","Gamma"
};

// Error handling struct
typedef struct _NRPyEOS_error_report_ {
  bool error;
  int error_key;
  char message[512];
} NRPyEOS_error_report;
"""

# Step 5: Interpolation helper functions
def gen_Cheader_interpolation_helpers(Ccodesdir):
    with open(os.path.join(Ccodesdir,"NRPyEOS_tabulated_helpers.h"),"w") as file:
        file.write(r"""
/*
 * (c) 2022 Leo Werneck
 *
 * This file contains modified functions from the original
 * helpers.hh file from the Zelmani eosdrivercxx repository.
 * Source: https://bitbucket.org/zelmani/eosdrivercxx
 */

//------------------------------------------
static inline __attribute__((always_inline))
int NRPyEOS_checkbounds(const NRPyEOS_params *restrict eos_params,
                        const double xrho,
                        const double xtemp,
                        const double xye) {

  // keyerr codes:
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 103 -- temp too high (if keytemp = 1)
  // 104 -- temp too low (if keytemp = 1)
  // 105 -- rho too high
  // 106 -- rho too low

  if(xrho > eos_params->eos_rhomax) {
    return 105;
  }
  if(xrho < eos_params->eos_rhomin) {
    return 106;
  }
  if(xye > eos_params->eos_yemax) {
    return 101;
  }
  if(xye < eos_params->eos_yemin) {
    // this is probably not pure and should be removed
    fprintf(stderr,"xye: %15.6E eos_yemin: %15.6E\n",xye,eos_params->eos_yemin);
    return 102;
  }
  if(xtemp > eos_params->eos_tempmax) {
    return 103;
  }
  if(xtemp < eos_params->eos_tempmin) {
    return 104;
  }
  return 0;
}
//------------------------------------------
static inline __attribute__((always_inline))
int NRPyEOS_checkbounds_kt0_noTcheck(const NRPyEOS_params *restrict eos_params,
                                     const double xrho,
                                     const double xye) {

  // keyerr codes:
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 105 -- rho too high
  // 106 -- rho too low

  if(xrho > eos_params->eos_rhomax) {
    return 105;
  }
  if(xrho < eos_params->eos_rhomin) {
    return 106;
  }
  if(xye > eos_params->eos_yemax) {
    return 101;
  }
  if(xye < eos_params->eos_yemin) {
    return 102;
  }
  return 0;
}
//------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_get_interp_spots(const NRPyEOS_params *restrict eos_params,
                              const double x,
                              const double y,
                              const double z,
                              double *restrict delx,
                              double *restrict dely,
                              double *restrict delz,
                              int *restrict idx) {

  int ix = 1 + (int)( (x - eos_params->logrho[0]  - 1.0e-10) * eos_params->drhoi );
  int iy = 1 + (int)( (y - eos_params->logtemp[0] - 1.0e-10) * eos_params->dtempi );
  int iz = 1 + (int)( (z - eos_params->yes[0]     - 1.0e-10) * eos_params->dyei );

  ix = MAX( 1, MIN( ix, eos_params->nrho -1 ) );
  iy = MAX( 1, MIN( iy, eos_params->ntemp-1 ) );
  iz = MAX( 1, MIN( iz, eos_params->nye  -1 ) );

  idx[0] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*(iy     + eos_params->ntemp*iz    ));
  idx[1] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*(iy     + eos_params->ntemp*iz    ));
  idx[2] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*((iy-1) + eos_params->ntemp*iz    ));
  idx[3] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*(iy     + eos_params->ntemp*(iz-1)));
  idx[4] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*((iy-1) + eos_params->ntemp*iz    ));
  idx[5] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*(iy     + eos_params->ntemp*(iz-1)));
  idx[6] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*((iy-1) + eos_params->ntemp*(iz-1)));
  idx[7] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*((iy-1) + eos_params->ntemp*(iz-1)));

  // set up aux vars for interpolation
  *delx = eos_params->logrho[ix]  - x;
  *dely = eos_params->logtemp[iy] - y;
  *delz = eos_params->yes[iz]     - z;

}
//------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_get_interp_spots_linT_low(const NRPyEOS_params *restrict eos_params,
                                       const double x,
                                       const double y,
                                       const double z,
                                       double *restrict delx,
                                       double *restrict dely,
                                       double *restrict delz,
                                       int *restrict idx) {

  int ix = 1 + (int)( (x - eos_params->logrho[0] - 1.0e-10) * eos_params->drhoi );
  int iy = 1;
  int iz = 1 + (int)( (z - eos_params->yes[0]    - 1.0e-10) * eos_params->dyei );

  ix = MAX( 1, MIN( ix, eos_params->nrho-1 ) );
  iz = MAX( 1, MIN( iz, eos_params->nye -1 ) );

  idx[0] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*(iy     + eos_params->ntemp*iz));
  idx[1] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*(iy     + eos_params->ntemp*iz));
  idx[2] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*((iy-1) + eos_params->ntemp*iz));
  idx[3] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*(iy     + eos_params->ntemp*(iz-1)));
  idx[4] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*((iy-1) + eos_params->ntemp*iz));
  idx[5] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*(iy     + eos_params->ntemp*(iz-1)));
  idx[6] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*((iy-1) + eos_params->ntemp*(iz-1)));
  idx[7] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*((iy-1) + eos_params->ntemp*(iz-1)));

  // set up aux vars for interpolation
  *delx = eos_params->logrho[ix] - x;
  *dely = eos_params->temp1      - y;
  *delz = eos_params->yes[iz]    - z;

}
//------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_get_interp_spots_linT_low_eps(const NRPyEOS_params *restrict eos_params,
                                           const double x,
                                           const double y,
                                           const double z,
                                           double *restrict delx,
                                           double *restrict dely,
                                           double *restrict delz,
                                           int *restrict idx) {

  int ix = 1 + (int)( (x - eos_params->logrho[0] - 1.0e-10) * eos_params->drhoi );
  int iy = 1;
  int iz = 1 + (int)( (z - eos_params->yes[0]    - 1.0e-10) * eos_params->dyei );

  ix = MAX( 1, MIN( ix, eos_params->nrho-1 ) );
  iz = MAX( 1, MIN( iz, eos_params->nye -1 ) );

  idx[0] = (ix     + eos_params->nrho*(iy     + eos_params->ntemp*iz));
  idx[1] = ((ix-1) + eos_params->nrho*(iy     + eos_params->ntemp*iz));
  idx[2] = (ix     + eos_params->nrho*((iy-1) + eos_params->ntemp*iz));
  idx[3] = (ix     + eos_params->nrho*(iy     + eos_params->ntemp*(iz-1)));
  idx[4] = ((ix-1) + eos_params->nrho*((iy-1) + eos_params->ntemp*iz));
  idx[5] = ((ix-1) + eos_params->nrho*(iy     + eos_params->ntemp*(iz-1)));
  idx[6] = (ix     + eos_params->nrho*((iy-1) + eos_params->ntemp*(iz-1)));
  idx[7] = ((ix-1) + eos_params->nrho*((iy-1) + eos_params->ntemp*(iz-1)));

  // set up aux vars for interpolation
  *delx = eos_params->logrho[ix] - x;
  *dely = eos_params->temp1      - y;
  *delz = eos_params->yes[iz]    - z;

}
//------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_linterp_one(const NRPyEOS_params *restrict eos_params,
                         const int *restrict idx,
                         const double delx,
                         const double dely,
                         const double delz,
                         double *restrict f,
                         const int iv) {

  // helper variables
  double fh[8], a[8];

  fh[0] = eos_params->alltables[iv+idx[0]];
  fh[1] = eos_params->alltables[iv+idx[1]];
  fh[2] = eos_params->alltables[iv+idx[2]];
  fh[3] = eos_params->alltables[iv+idx[3]];
  fh[4] = eos_params->alltables[iv+idx[4]];
  fh[5] = eos_params->alltables[iv+idx[5]];
  fh[6] = eos_params->alltables[iv+idx[6]];
  fh[7] = eos_params->alltables[iv+idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = eos_params->drhoi       * ( fh[1] - fh[0] );
  a[2] = eos_params->dtempi      * ( fh[2] - fh[0] );
  a[3] = eos_params->dyei        * ( fh[3] - fh[0] );
  a[4] = eos_params->drhotempi   * ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = eos_params->drhoyei     * ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = eos_params->dtempyei    * ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = eos_params->drhotempyei * ( fh[7] - fh[0] + fh[1] + fh[2] +
                                     fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0]
     + a[1] * delx
     + a[2] * dely
     + a[3] * delz
     + a[4] * delx * dely
     + a[5] * delx * delz
     + a[6] * dely * delz
     + a[7] * delx * dely * delz;

}
//------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_linterp_one_linT_low(const NRPyEOS_params *restrict eos_params,
                                  const int *restrict idx,
                                  const double delx,
                                  const double dely,
                                  const double delz,
                                  double *restrict f,
                                  const int iv) {

  // helper variables
  double fh[8], a[8];

  fh[0] = eos_params->alltables[iv+idx[0]];
  fh[1] = eos_params->alltables[iv+idx[1]];
  fh[2] = eos_params->alltables[iv+idx[2]];
  fh[3] = eos_params->alltables[iv+idx[3]];
  fh[4] = eos_params->alltables[iv+idx[4]];
  fh[5] = eos_params->alltables[iv+idx[5]];
  fh[6] = eos_params->alltables[iv+idx[6]];
  fh[7] = eos_params->alltables[iv+idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = eos_params->drhoi          * ( fh[1] - fh[0] );
  a[2] = eos_params->dlintempi      * ( fh[2] - fh[0] );
  a[3] = eos_params->dyei           * ( fh[3] - fh[0] );
  a[4] = eos_params->drholintempi   * ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = eos_params->drhoyei        * ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = eos_params->dlintempyei    * ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = eos_params->drholintempyei * ( fh[7] - fh[0] + fh[1] + fh[2] +
                                        fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0]
     + a[1] * delx
     + a[2] * dely
     + a[3] * delz
     + a[4] * delx * dely
     + a[5] * delx * delz
     + a[6] * dely * delz
     + a[7] * delx * dely * delz;

}
//------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_linterp_one_linT_low_eps(const NRPyEOS_params *restrict eos_params,
                                      const int *restrict idx,
                                      const double delx,
                                      const double dely,
                                      const double delz,
                                      double *restrict f) {

  // helper variables
  double fh[8], a[8];

  fh[0] = eos_params->epstable[idx[0]];
  fh[1] = eos_params->epstable[idx[1]];
  fh[2] = eos_params->epstable[idx[2]];
  fh[3] = eos_params->epstable[idx[3]];
  fh[4] = eos_params->epstable[idx[4]];
  fh[5] = eos_params->epstable[idx[5]];
  fh[6] = eos_params->epstable[idx[6]];
  fh[7] = eos_params->epstable[idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = eos_params->drhoi          * ( fh[1] - fh[0] );
  a[2] = eos_params->dlintempi      * ( fh[2] - fh[0] );
  a[3] = eos_params->dyei           * ( fh[3] - fh[0] );
  a[4] = eos_params->drholintempi   * ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = eos_params->drhoyei        * ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = eos_params->dlintempyei    * ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = eos_params->drholintempyei * ( fh[7] - fh[0] + fh[1] + fh[2] +
                                        fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0]
     + a[1] * delx
     + a[2] * dely
     + a[3] * delz
     + a[4] * delx * dely
     + a[5] * delx * delz
     + a[6] * dely * delz
     + a[7] * delx * dely * delz;

}
//------------------------------------------
static inline __attribute__((always_inline))
double NRPyEOS_linterp2D(const double *restrict xs,
                         const double *restrict ys,
                         const double *restrict fs,
                         const double x,
                         const double y) {

  //  2     3
  //
  //  0     1
  //
  // first interpolate in x between 0 and 1, 2 and 3
  // then interpolate in y
  // assume rectangular grid

  double dxi = 1./(xs[1]-xs[0]);
  double dyi = 1./(ys[1]-ys[0]); // x*1./y uses faster instructions than x/y
  double t1 = (fs[1]-fs[0])*dxi * (x - xs[0]) + fs[0];
  double t2 = (fs[3]-fs[2])*dxi * (x - xs[0]) + fs[2];

  return (t2 - t1)*dyi * (y-ys[0]) + t1;
}
//------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_bisection(const NRPyEOS_params *restrict eos_params,
                       const double lr,
                       const double lt0,
                       const double ye,
                       const double leps0,
                       const double prec,
                       double *restrict ltout,
                       const int iv,
                       int *restrict keyerrt) {
  // iv is the index of the variable we do the bisection on

  int bcount = 0;
  int maxbcount = 80;
  int itmax = 50;

  const double dlt0p = log(1.1);
  const double dlt0m = log(0.9);
  const double dltp  = log(1.2);
  const double dltm  = log(0.8);

  double leps0_prec = fabs(leps0*prec);

  // temporary local vars
  double lt, lt1, lt2;
  double ltmin = eos_params->logtemp[0];
  double ltmax = eos_params->logtemp[eos_params->ntemp-1];
  double f1,f2,fmid,dlt,ltmid;
  double f1a = 0.0;
  double f2a = 0.0;
  double delx,dely,delz;
  int idx[8];

  // LSMOD (Modification made by Lorenzo Sala)
  // LSMOD: The following lines calculate eps in
  //        f2a = eps(rho,Tmin, Ye) and f1a = eps(rho,Tmax,Ye)
  NRPyEOS_get_interp_spots(eos_params,lr,ltmax,ye,&delx,&dely,&delz,idx);
  NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&f1a,iv);
  NRPyEOS_get_interp_spots(eos_params,lr,ltmin,ye,&delx,&dely,&delz,idx);
  NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&f2a,iv);

  // prepare
  // check if your energy is actually tabulated at this rho and ye.
  // f2a is the energy evaluated at ltmin, so it is the minimum energy tabulated
  // at this rho ad ye.
  // If leps0 <= f2a, then ltout is likely to be the minimum temperature tabulated.
  if(leps0 <= f2a) { // + 1.0E-6
    *ltout = ltmin;
    return;
  }

  /* // If leps0 >= f1a, then ltout is likely to be the maximum temperature tabulated.
     if(leps0 >= f1a) { // + 1.0E-6
     *ltout = ltmax;
     return;
     } */

  // otherwise, proceed finding extrema for applying bisection method.
  lt = lt0;
  lt1 = MIN(lt0 + dlt0p,ltmax);
  lt2 = MAX(lt0 + dlt0m,ltmin);

  NRPyEOS_get_interp_spots(eos_params,lr,lt1,ye,&delx,&dely,&delz,idx);
  NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&f1a,iv);

  NRPyEOS_get_interp_spots(eos_params,lr,lt2,ye,&delx,&dely,&delz,idx);
  NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&f2a,iv);

  f1=f1a-leps0;
  f2=f2a-leps0;

  // iterate until we bracket the right eps, but enforce
  // dE/dt > 0, so eps(lt1) > eps(lt2)
  while(f1*f2 >= 0.0) {
    lt1 = MIN(lt1 + dltp,ltmax);
    lt2 = MAX(lt2 + dltm,ltmin);
    NRPyEOS_get_interp_spots(eos_params,lr,lt1,ye,&delx,&dely,&delz,idx);
    NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&f1a,iv);

    NRPyEOS_get_interp_spots(eos_params,lr,lt2,ye,&delx,&dely,&delz,idx);
    NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&f2a,iv);

    f1=f1a-leps0;
    f2=f2a-leps0;

#if DEBUG
    fprintf(stderr,"bisection bracketing it %d, f1: %15.6E, f2: %15.6E, lt1: %15.6E, lt2: %15.6E, f1a: %18.11E, f2a: %18.11E leps0: %18.11E\n",
            bcount,f1,f2,lt1,lt2,f1a,f2a,leps0);
#endif

    bcount++;
    if(bcount >= maxbcount) {
#if DEBUG
      fprintf(stderr,"bcount out of range it %d, lr: %15.6E, lt1: %15.6E, lt2: %15.6E, f1a: %18.11E, f2a: %18.11E leps0: %18.11E, ye: %15.6E\n",
              bcount,lr,lt1,lt2,f1a,f2a,leps0,ye);
#endif
      *keyerrt = 667;
      return;
    }
  } // while

  if(f1 < 0.0) {
    lt = lt1;
    dlt = lt2 - lt1;
  } else {
    lt = lt2;
    dlt = lt1 - lt2;
  }

#if DEBUG
  fprintf(stderr,"bisection step 2 it -1, fmid: %15.6E ltmid: %15.6E dlt: %15.6E\n",
          f2,lt,dlt);
  fprintf(stderr,"ltmax: %15.6E\n",ltmax);
#endif

  int it;
  for(it=0;it<itmax;it++) {
    dlt = dlt * 0.5;
    ltmid = lt + dlt;
    NRPyEOS_get_interp_spots(eos_params,lr,ltmid,ye,&delx,&dely,&delz,idx);
    NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&f2a,iv);

    fmid=f2a-leps0;
    if(fmid <= 0.0) lt=ltmid;
#if DEBUG
    fprintf(stderr,"bisection step 2 it %d, fmid: %15.6E f2a: %15.6E lt: %15.6E ltmid: %15.6E dlt: %15.6E\n",
            it,fmid,f2a,lt,ltmid,dlt);
#endif

    if(fabs(leps0-f2a) <= leps0_prec) {
      *ltout = ltmid;
      return;
    }
  } // for it = 0

  *keyerrt = 667;
  return;
} // bisection
  //------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_findtemp_from_any( const NRPyEOS_params *restrict eos_params,
                                const int tablevar_key,
                                const double lr,
                                const double lt0,
                                const double ye,
                                const double tablevar_in,
                                const double prec,
                                double *restrict ltout,
                                int *keyerrt ) {

  // local variables
  const int itmax = 200; // use at most 10 iterations, then go to bisection
  double dtablevardlti; // 1 / derivative dlogeps/dlogT
  double ldt;
  double tablevar; // temp vars for eps
  double ltn; // temp vars for temperature
  const double ltmax = eos_params->logtemp[eos_params->ntemp-1]; // max temp
  const double ltmin = eos_params->logtemp[0]; // min temp
  int it = 0;

  // setting up some vars
  *keyerrt  = 0;
  double lt = lt0;

  // step 1: do we already have the right temperature
  int idx[8];
  double delx,dely,delz;
  NRPyEOS_get_interp_spots(eos_params,lr,lt,ye,&delx,&dely,&delz,idx);
  NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&tablevar,tablevar_key);

  // TODO: profile this to see which outcome is more likely
  if(fabs(tablevar-tablevar_in) < prec*fabs(tablevar_in)) {
    *ltout = lt0;
    return;
  }

  double oerr = 1.0e90;
  double fac  = 1.0;
  const int irho = MIN(MAX(1 + (int)(( lr - eos_params->logrho[0] - 1.0e-12) * eos_params->drhoi),1),eos_params->nrho-1);
  const int iye  = MIN(MAX(1 + (int)(( ye - eos_params->yes[0]    - 1.0e-12) * eos_params->dyei ),1),eos_params->nye -1);

  /* ******* if temp low for high density, switch directly to bisection.
     Verifying Newton-Raphson result evaluating the derivative.
     The variable shouldgotobisection will be modified accordingly
     to the value of derivative of eps wrt temp ******* */
  bool shouldgotobisection = false; // LSMOD
  while(it < itmax && shouldgotobisection == false) {
    it++;

    // step 2: check if the two bounding values of the temperature
    //         give eps values that enclose the new eps.
    const int itemp = MIN(MAX(1 + (int)(( lt - eos_params->logtemp[0] - 1.0e-12) * eos_params->dtempi),1),eos_params->ntemp-1);

    double tablevart1, tablevart2;
    // lower temperature
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = tablevar_key + NRPyEOS_ntablekeys*(irho-1 + eos_params->nrho*((itemp-1) + eos_params->ntemp*(iye-1)));
      fs[0]   = eos_params->alltables[ifs];
      // point 1
      ifs     = tablevar_key + NRPyEOS_ntablekeys*(irho   + eos_params->nrho*((itemp-1) + eos_params->ntemp*(iye-1)));
      fs[1]   = eos_params->alltables[ifs];
      // point 2
      ifs     = tablevar_key + NRPyEOS_ntablekeys*(irho-1 + eos_params->nrho*((itemp-1) + eos_params->ntemp*(iye)));
      fs[2]   = eos_params->alltables[ifs];
      // point 3
      ifs     = tablevar_key + NRPyEOS_ntablekeys*(irho   + eos_params->nrho*((itemp-1) + eos_params->ntemp*(iye)));
      fs[3]   = eos_params->alltables[ifs];

      tablevart1 = NRPyEOS_linterp2D(&eos_params->logrho[irho-1],&eos_params->yes[iye-1], fs, lr, ye);
    }
    // upper temperature
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = tablevar_key + NRPyEOS_ntablekeys*(irho-1 + eos_params->nrho*((itemp) + eos_params->ntemp*(iye-1)));
      fs[0]   = eos_params->alltables[ifs];
      // point 1
      ifs     = tablevar_key + NRPyEOS_ntablekeys*(irho   + eos_params->nrho*((itemp) + eos_params->ntemp*(iye-1)));
      fs[1]   = eos_params->alltables[ifs];
      // point 2
      ifs     = tablevar_key + NRPyEOS_ntablekeys*(irho-1 + eos_params->nrho*((itemp) + eos_params->ntemp*(iye)));
      fs[2]   = eos_params->alltables[ifs];
      // point 3
      ifs     = tablevar_key + NRPyEOS_ntablekeys*(irho   + eos_params->nrho*((itemp) + eos_params->ntemp*(iye)));
      fs[3]   = eos_params->alltables[ifs];

      tablevart2 = NRPyEOS_linterp2D(&eos_params->logrho[irho-1],&eos_params->yes[iye-1], fs, lr, ye);
    }

    // Check if we are already bracketing the input internal
    // energy. If so, interpolate for new T.
    if((tablevar_in - tablevart1) * (tablevar_in - tablevart2) <= 0.) {

      *ltout = (eos_params->logtemp[itemp]-eos_params->logtemp[itemp-1]) / (tablevart2 - tablevart1) *
        (tablevar_in - tablevart1) + eos_params->logtemp[itemp-1];

      return;
    }

    // well, then do a Newton-Raphson step
    // first, guess the derivative
    dtablevardlti = (eos_params->logtemp[itemp]-eos_params->logtemp[itemp-1])/(tablevart2-tablevart1);
    ldt = -(tablevar - tablevar_in) * dtablevardlti * fac;

    //LSMOD: too large a dlt means that the energy dependence on the temperature
    //       is weak ==> We'd better try bisection.
    //       Factor 1/12.0 come from tests by LSMOD
    //       This is done in order to limit the "velocity" of T variation
    //       given by Newton-Raphson.
    if(ldt > (ltmax-ltmin) / 12.0 ) shouldgotobisection = true;

    ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
    lt = ltn;

    NRPyEOS_get_interp_spots(eos_params,lr,lt,ye,&delx,&dely,&delz,idx);
    NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&tablevar,tablevar_key);

    // drive the thing into the right direction
    double err = fabs(tablevar-tablevar_in);
    if(oerr < err) fac *= 0.9;
    oerr = err;

    if(err < prec*fabs(tablevar_in)) {
      *ltout = lt;
      return;
    }

  } // while(it < itmax)

    // try bisection
  NRPyEOS_bisection(eos_params,lr,lt0,ye,tablevar_in,prec,ltout,tablevar_key,keyerrt);

  return;
}
""")

# Step 6: General wrapper functions for interpolation
# Step 6.a: Wrapper for when T is known
def Cfunc_general_wrapper_known_T():
    includes   = ["NRPy_basic_defines.h","NRPy_function_prototypes.h","NRPyEOS_tabulated_helpers.h"]
    desc       = "(c) 2022 Leo Werneck"
    c_type     = "void"
    name       = "NRPyEOS_from_rho_Ye_T_interpolate_n_quantities"
    indent     = param_indentation(c_type,name)
    params     = "const NRPyEOS_params *restrict eos_params,\n"
    params    += indent+"const int n,\n"
    params    += indent+"const double rho,\n"
    params    += indent+"const double Y_e,\n"
    params    += indent+"const double T,\n"
    params    += indent+"const int *restrict tablevars_keys,\n"
    params    += indent+"double *restrict tablevars,\n"
    params    += indent+"NRPyEOS_error_report *restrict report"
    body       = r"""
  // This function will interpolate n table quantities from
  // (rho,Ye,T). It replaces EOS_Omni calls with keytemp = 1
  if( n > NRPyEOS_ntablekeys ) {
    fprintf(stderr,"(NRPyEOS) from_rho_Ye_T_interpolate_n_quantities: number of quantities exceed maximum allowed: %d > %d. ABORTING.",
            n,NRPyEOS_ntablekeys);
  }

  // Start by assuming no errors
  report->error = false;

  // Check table bounds for input variables
  report->error_key = NRPyEOS_checkbounds(eos_params,rho,T,Y_e);
  if( report->error_key != 0 ) {
    // This should never happen, because we enforce
    // limits before calling this function
    sprintf(report->message,"from_rho_Ye_T_interpolate_n_quantities: problem with checkbounds");
    report->error = true;
    return;
  }

  // Get interpolation spots
  int idx[8];
  double delx,dely,delz;
  const double lr = log(rho);
  const double lt = log(T);
  NRPyEOS_get_interp_spots(eos_params,lr,lt,Y_e,&delx,&dely,&delz,idx);

  for(int i=0;i<n;i++) {
    // Now perform the interpolations
    int key = tablevars_keys[i];
    double tablevar_out;
    NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&tablevar_out,key);

    // We have the result, but we must convert appropriately.
    // The only edge cases are P and eps, for which we obtain
    // log(P) and log(eps+eps0). We must check for them here
    if( key == NRPyEOS_press_key ) {
      tablevar_out = exp(tablevar_out);
    }
    else if( key == NRPyEOS_eps_key ) {
      tablevar_out = exp(tablevar_out) - eos_params->energy_shift;
    }

    // Then update tablevars
    tablevars[i] = tablevar_out;
  }
"""
    outC.add_to_Cfunction_dict(includes=includes,desc=desc,c_type=c_type,name=name,
                               params=params,body=body,enableCparameters=False)

# Step 6.b: Wrapper for when T is unknown
def Cfunc_general_wrapper_unknown_T():
    includes   = ["NRPy_basic_defines.h","NRPy_function_prototypes.h","NRPyEOS_tabulated_helpers.h"]
    desc       = "(c) 2022 Leo Werneck"
    c_type     = "void"
    name       = "NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities"
    indent     = param_indentation(c_type,name)
    params     = "const NRPyEOS_params *restrict eos_params,\n"
    params    += indent+"const int n,\n"
    params    += indent+"const double prec,\n"
    params    += indent+"const double rho,\n"
    params    += indent+"const double Y_e,\n"
    params    += indent+"const double tablevar_in,\n"
    params    += indent+"const int tablevar_in_key,\n"
    params    += indent+"const int *restrict tablevars_keys,\n"
    params    += indent+"double *restrict tablevars,\n"
    params    += indent+"double *restrict T,\n"
    params    += indent+"NRPyEOS_error_report *restrict report"
    body       = r"""
  // This function will interpolate n table quantities from
  // (rho,Ye,aux). It replaces EOS_Omni calls with keytemp != 1
  if( n > NRPyEOS_ntablekeys ) {
    fprintf(stderr,"(NRPyEOS) NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities: number of quantities exceed maximum allowed: %d > %d. ABORTING.",
            n,NRPyEOS_ntablekeys);
  }

  // Check table bounds for input variables
  report->error_key = NRPyEOS_checkbounds_kt0_noTcheck(eos_params,rho,Y_e);
  if( report->error_key != 0 ) {
    // This should never happen, because we enforce
    // limits before calling this function
    sprintf(report->message,"NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities: problem with checkbounds_kt0_noTcheck");
    report->error = true;
    return;
  }

  // First step is to recover the temperature. The variable
  // tablevar_in is the one used in the temperature recovery.
  // For example, if tablevar_in = eps, then we recover T
  // using (rho,Ye,eps).
  double aux = tablevar_in;

  if( tablevar_in_key == NRPyEOS_press_key ) {
    // If aux = P, then we need log(P).
    aux = log(aux);
  }
  else if( tablevar_in_key == NRPyEOS_eps_key ) {
    // If aux = eps, then we need log(eps+eps0).
    // Compute eps+eps0
    aux += eos_params->energy_shift;
    // At this point, aux *must* be positive. If not, error out.
    if( aux < 0.0 ) {
      fprintf(stderr,"(NRPyEOS) NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities: found eps+energy_shift < 0.0 (%e). ABORTING.",
              aux);
    }
    // Compute log(eps+eps0)
    aux = log(aux);
  }

  // Now compute the temperature
  const double lr  = log(rho);
  const double lt0 = log(*T);
  double lt        = 0.0;
  int keyerr=0;
  NRPyEOS_findtemp_from_any(eos_params,tablevar_in_key,lr,lt0,Y_e,aux,prec,&lt,&keyerr);

  // Now set the temperature
  *T = exp(lt);

  // Then interpolate the quantities we want from (rho,Ye,T)
  int anyerr=0;
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(eos_params,n,rho,Y_e,*T,tablevars_keys,tablevars,report);
  report->error_key = keyerr;
  report->error     = anyerr;
"""
    outC.add_to_Cfunction_dict(includes=includes,desc=desc,c_type=c_type,name=name,
                               params=params,body=body,enableCparameters=False)

# Step 7: EOS table reader & memory management
# Step 7.a: Table reader and memory allocation
def Cfunc_read_table_set_EOS_params():
    includes   = ["NRPy_basic_defines.h","NRPy_function_prototypes.h"]
    desc       = "(c) 2022 Leo Werneck"
    c_type     = "void"
    name       = "NRPyEOS_readtable_set_EOS_params"
    params     = "const char *nuceos_table_name, NRPyEOS_params *restrict eos_params"
    prefunc    = r"""
// mini NoMPI
#ifdef HAVE_CAPABILITY_MPI
#include <mpi.h>
#define BCAST(buffer, size) MPI_Bcast(buffer, size, MPI_BYTE, my_reader_process, MPI_COMM_WORLD)
#else
#define BCAST(buffer, size) do { /* do nothing */ } while(0)
#endif

// If on the IO proc (doIO == True) actually perform HDF5 IO, catch possible
// HDF5 errors
#define HDF5_DO_IO(fn_call)                                              \
  {                                                                      \
    int _error_code = fn_call;                                           \
    if (_error_code < 0) {                                               \
      fprintf(stderr,"(NRPyEOS) HDF5 call '%s' returned error code %d",  \
                  #fn_call, _error_code);                                \
    }                                                                    \
  }
"""
    body       = r"""
  fprintf(stderr,"(NRPyEOS) *******************************\n");
  fprintf(stderr,"(NRPyEOS) Reading EOS table from file:\n");
  fprintf(stderr,"(NRPyEOS) %s\n",nuceos_table_name);
  fprintf(stderr,"(NRPyEOS) *******************************\n");

  hid_t file;
  HDF5_DO_IO(file = H5Fopen(nuceos_table_name, H5F_ACC_RDONLY, H5P_DEFAULT));

// Use these two defines to easily read in a lot of variables in the same way
// The first reads in one variable of a given type completely
#define READ_BCAST_EOS_HDF5(NAME,VAR,TYPE,MEM,NELEMS)                   \
  do {                                                                  \
    hid_t dataset;                                                      \
    HDF5_DO_IO(dataset = H5Dopen(file, NAME, H5P_DEFAULT));             \
    HDF5_DO_IO(H5Dread(dataset, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR)); \
    BCAST (VAR, sizeof(*(VAR))*(NELEMS));                               \
    HDF5_DO_IO(H5Dclose(dataset));                                      \
  } while (0)
// The second reads a given variable into a hyperslab of the alltables_temp array
#define READ_BCAST_EOSTABLE_HDF5(NAME,OFF,DIMS)                          \
  do {                                                                   \
    READ_BCAST_EOS_HDF5(NAME,&alltables_temp[(OFF)*(DIMS)[1]],H5T_NATIVE_DOUBLE,H5S_ALL,(DIMS)[1]); \
  } while (0)

  // Read size of tables
  READ_BCAST_EOS_HDF5("pointsrho",  &eos_params->nrho,  H5T_NATIVE_INT, H5S_ALL, 1);
  READ_BCAST_EOS_HDF5("pointstemp", &eos_params->ntemp, H5T_NATIVE_INT, H5S_ALL, 1);
  READ_BCAST_EOS_HDF5("pointsye",   &eos_params->nye,   H5T_NATIVE_INT, H5S_ALL, 1);

  // Allocate memory for tables
  double* alltables_temp;
  if (!(alltables_temp = (double*)malloc(eos_params->nrho * eos_params->ntemp * eos_params->nye * NRPyEOS_ntablekeys * sizeof(double)))) {
    fprintf(stderr,"(NRPyEOS) Cannot allocate memory for EOS table");
  }
  if (!(eos_params->logrho = (double*)malloc(eos_params->nrho * sizeof(double)))) {
    fprintf(stderr,"(NRPyEOS) Cannot allocate memory for EOS table");
  }
  if (!(eos_params->logtemp = (double*)malloc(eos_params->ntemp * sizeof(double)))) {
    fprintf(stderr,"(NRPyEOS) Cannot allocate memory for EOS table");
  }
  if (!(eos_params->yes = (double*)malloc(eos_params->nye * sizeof(double)))) {
    fprintf(stderr,"(NRPyEOS) Cannot allocate memory for EOS table");
  }

  // Prepare HDF5 to read hyperslabs into alltables_temp
  hsize_t table_dims[2] = {NRPyEOS_ntablekeys, (hsize_t)eos_params->nrho * eos_params->ntemp * eos_params->nye};
  hid_t mem3 =  H5Screate_simple(2, table_dims, NULL);

  // Read alltables_temp
  READ_BCAST_EOSTABLE_HDF5("logpress",  0, table_dims);
  READ_BCAST_EOSTABLE_HDF5("logenergy", 1, table_dims);
  READ_BCAST_EOSTABLE_HDF5("entropy",   2, table_dims);
  READ_BCAST_EOSTABLE_HDF5("munu",      3, table_dims);
  READ_BCAST_EOSTABLE_HDF5("cs2",       4, table_dims);
  READ_BCAST_EOSTABLE_HDF5("dedt",      5, table_dims);
  READ_BCAST_EOSTABLE_HDF5("dpdrhoe",   6, table_dims);
  READ_BCAST_EOSTABLE_HDF5("dpderho",   7, table_dims);
  // chemical potentials
  READ_BCAST_EOSTABLE_HDF5("muhat",     8, table_dims);
  READ_BCAST_EOSTABLE_HDF5("mu_e",      9, table_dims);
  READ_BCAST_EOSTABLE_HDF5("mu_p",     10, table_dims);
  READ_BCAST_EOSTABLE_HDF5("mu_n",     11, table_dims);
  // compositions
  READ_BCAST_EOSTABLE_HDF5("Xa",       12, table_dims);
  READ_BCAST_EOSTABLE_HDF5("Xh",       13, table_dims);
  READ_BCAST_EOSTABLE_HDF5("Xn",       14, table_dims);
  READ_BCAST_EOSTABLE_HDF5("Xp",       15, table_dims);
  // average nucleus
  READ_BCAST_EOSTABLE_HDF5("Abar",     16, table_dims);
  READ_BCAST_EOSTABLE_HDF5("Zbar",     17, table_dims);
  // Gamma
  READ_BCAST_EOSTABLE_HDF5("gamma",    18, table_dims);

  // Read additional tables and variables
  READ_BCAST_EOS_HDF5("logrho",       eos_params->logrho,        H5T_NATIVE_DOUBLE, H5S_ALL, eos_params->nrho);
  READ_BCAST_EOS_HDF5("logtemp",      eos_params->logtemp,       H5T_NATIVE_DOUBLE, H5S_ALL, eos_params->ntemp);
  READ_BCAST_EOS_HDF5("ye",           eos_params->yes,           H5T_NATIVE_DOUBLE, H5S_ALL, eos_params->nye);
  READ_BCAST_EOS_HDF5("energy_shift", &eos_params->energy_shift, H5T_NATIVE_DOUBLE, H5S_ALL, 1);

  HDF5_DO_IO(H5Sclose(mem3));
  HDF5_DO_IO(H5Fclose(file));

  // change ordering of alltables array so that
  // the table kind is the fastest changing index
  if (!(eos_params->alltables = (double*)malloc(eos_params->nrho * eos_params->ntemp * eos_params->nye * NRPyEOS_ntablekeys
                                                * sizeof(double)))) {
    fprintf(stderr,"(NRPyEOS) Cannot allocate memory for EOS table");
  }
  for(int iv = 0;iv<NRPyEOS_ntablekeys;iv++)
    for(int k = 0; k<eos_params->nye;k++)
      for(int j = 0; j<eos_params->ntemp; j++)
    for(int i = 0; i<eos_params->nrho; i++) {
      int indold = i + eos_params->nrho*(j + eos_params->ntemp*(k + eos_params->nye*iv));
      int indnew = iv + NRPyEOS_ntablekeys*(i + eos_params->nrho*(j + eos_params->ntemp*k));
      eos_params->alltables[indnew] = alltables_temp[indold];
    }

  // free memory of temporary array
  free(alltables_temp);

  // convert units, convert logs to natural log
  // The latter is great, because exp() is way faster than pow()
  // pressure
  eos_params->energy_shift = eos_params->energy_shift * EPSGF;
  for(int i=0;i<eos_params->nrho;i++) {
    // rewrite:
    //logrho[i] = log(pow(10.0,logrho[i]) * RHOGF);
    // by using log(a^b*c) = b*log(a)+log(c)
    eos_params->logrho[i] = eos_params->logrho[i] * log(10.) + log(RHOGF);
  }

  for(int i=0;i<eos_params->ntemp;i++) {
    //logtemp[i] = log(pow(10.0,logtemp[i]));
    eos_params->logtemp[i] = eos_params->logtemp[i]*log(10.0);
  }

  // allocate epstable; a linear-scale eps table
  // that allows us to extrapolate to negative eps
  if (!(eos_params->epstable = (double*)malloc(eos_params->nrho * eos_params->ntemp * eos_params->nye
                                               * sizeof(double)))) {
    fprintf(stderr,"(NRPyEOS) Cannot allocate memory for eps table\n");
  }

  // convert units
  for(int i=0;i<eos_params->nrho*eos_params->ntemp*eos_params->nye;i++) {

    { // pressure
      int idx = 0 + NRPyEOS_ntablekeys*i;
      eos_params->alltables[idx] = eos_params->alltables[idx] * log(10.0) + log(PRESSGF);
    }

    { // eps
      int idx = 1 + NRPyEOS_ntablekeys*i;
      eos_params->alltables[idx] = eos_params->alltables[idx] * log(10.0) + log(EPSGF);
      eos_params->epstable[i] = exp(eos_params->alltables[idx]);
    }

    { // cs2
      int idx = 4 + NRPyEOS_ntablekeys*i;
      eos_params->alltables[idx] *= LENGTHGF*LENGTHGF/TIMEGF/TIMEGF;
    }

    { // dedT
      int idx = 5 + NRPyEOS_ntablekeys*i;
      eos_params->alltables[idx] *= EPSGF;
    }

    { // dpdrhoe
      int idx = 6 + NRPyEOS_ntablekeys*i;
      eos_params->alltables[idx] *= PRESSGF/RHOGF;
    }

    { // dpderho
      int idx = 7 + NRPyEOS_ntablekeys*i;
      eos_params->alltables[idx] *= PRESSGF/EPSGF;
    }

  }

  eos_params->temp0 = exp(eos_params->logtemp[0]);
  eos_params->temp1 = exp(eos_params->logtemp[1]);

  // set up some vars
  eos_params->dtemp  = (eos_params->logtemp[eos_params->ntemp-1] - eos_params->logtemp[0]) / (1.0*(eos_params->ntemp-1));
  eos_params->dtempi = 1.0/eos_params->dtemp;

  eos_params->dlintemp = eos_params->temp1-eos_params->temp0;
  eos_params->dlintempi = 1.0/eos_params->dlintemp;

  eos_params->drho  = (eos_params->logrho[eos_params->nrho-1] - eos_params->logrho[0]) / (1.0*(eos_params->nrho-1));
  eos_params->drhoi = 1.0/eos_params->drho;

  eos_params->dye  = (eos_params->yes[eos_params->nye-1] - eos_params->yes[0]) / (1.0*(eos_params->nye-1));
  eos_params->dyei = 1.0/eos_params->dye;

  eos_params->drhotempi      = eos_params->drhoi     * eos_params->dtempi;
  eos_params->drholintempi   = eos_params->drhoi     * eos_params->dlintempi;
  eos_params->drhoyei        = eos_params->drhoi     * eos_params->dyei;
  eos_params->dtempyei       = eos_params->dtempi    * eos_params->dyei;
  eos_params->dlintempyei    = eos_params->dlintempi * eos_params->dyei;
  eos_params->drhotempyei    = eos_params->drhoi     * eos_params->dtempi    * eos_params->dyei;
  eos_params->drholintempyei = eos_params->drhoi     * eos_params->dlintempi * eos_params->dyei;

  eos_params->eos_rhomax = exp(eos_params->logrho[eos_params->nrho-1]);
  eos_params->eos_rhomin = exp(eos_params->logrho[0]);

  eos_params->eos_tempmax = exp(eos_params->logtemp[eos_params->ntemp-1]);
  eos_params->eos_tempmin = exp(eos_params->logtemp[0]);

  eos_params->eos_yemax = eos_params->yes[eos_params->nye-1];
  eos_params->eos_yemin = eos_params->yes[0];
"""
    outC.add_to_Cfunction_dict(includes=includes,desc=desc,c_type=c_type,name=name,
                               params=params,body=body,enableCparameters=False,prefunc=prefunc)

# Step 7.b: Memory deallocation
def Cfunc_free_memory():
    includes   = ["NRPy_basic_defines.h"]
    desc       = "(c) 2022 Leo Werneck"
    c_type     = "void"
    name       = "NRPyEOS_free_memory"
    params     = "NRPyEOS_params *restrict eos_params"
    body       = r"""
 fprintf(stderr,"(NRPyEOS) *******************************\n");
 fprintf(stderr,"(NRPyEOS) Freeing up memory.\n");

  // Free memory allocated for the table
  free(eos_params->logrho);
  free(eos_params->logtemp);
  free(eos_params->yes);
  free(eos_params->alltables);
  free(eos_params->epstable);

 fprintf(stderr,"(NRPyEOS) All done!\n");
 fprintf(stderr,"(NRPyEOS) *******************************\n");
"""
    outC.add_to_Cfunction_dict(includes=includes,desc=desc,c_type=c_type,name=name,
                               params=params,body=body,enableCparameters=False)

# Step 7: Add all C functions to the dictionary
def NRPyEOS_generate_interpolators_and_add_all_Cfuncs_to_dict(Ccodesdir,known_T_params_list=None,unknown_T_auxvars_and_params_list=None):
    # Step 7.a: Functions for which the temperature is known
    if known_T_params_list is not None:
        for param_list in known_T_params_list:
            Cfunc_known_T(param_list)

    # Step 7.b: Functions for which the temperature is unknown
    if unknown_T_auxvars_and_params_list is not None:
        for param_list in unknown_T_auxvars_and_params_list:
            Cfunc_unknown_T(param_list[0],param_list[1])

    # Step 7.c: Interpolation helpers
    gen_Cheader_interpolation_helpers(Ccodesdir)

    # Step 7.d: General interpolation wrappers
    Cfunc_general_wrapper_known_T()
    Cfunc_general_wrapper_unknown_T()

    # Step 7.e: Table reader and memory allocation
    Cfunc_read_table_set_EOS_params()

    # Step 7.f: Memory deallocation
    Cfunc_free_memory()
