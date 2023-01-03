# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

from outputC import outCfunction, lhrh, add_to_Cfunction_dict, outC_function_outdir_dict, outC_function_dict, outC_function_prototype_dict, outC_function_master_list, outC_function_element # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import loop as lp                # NRPy+: Generate C code loops
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface

thismodule = __name__

par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER",2)

out_dir = os.path.join("GiRaFFE_standalone_Ccodes")
cmd.mkdir(out_dir)

CoordSystem = "Cartesian"

par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)
rfm.reference_metric() # Create ReU, ReDD needed for rescaling B-L initial data, generating BSSN RHSs, etc.

# Default Kreiss-Oliger dissipation strength
default_KO_strength = 0.1
diss_strength = par.Cparameters("REAL", thismodule, "diss_strength", default_KO_strength)

outCparams = "outCverbose=False,CSE_sorting=none"

import GRHD.equations as GRHD    # NRPy+: Generate general relativistic hydrodynamics equations
import GRFFE.equations as GRFFE  # NRPy+: Generate general relativisitic force-free electrodynamics equations

gammaDD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL","gammaDD","sym01",DIM=3)
betaU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","betaU",DIM=3)
alpha = gri.register_gridfunctions("AUXEVOL","alpha")
AD = ixp.register_gridfunctions_for_single_rank1("EVOL","AD")
BU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","BU")
ValenciavU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","ValenciavU")
psi6Phi = gri.register_gridfunctions("EVOL","psi6Phi")
StildeD = ixp.register_gridfunctions_for_single_rank1("EVOL","StildeD")

# We will pass values of the gridfunction on the cell faces into the function. This requires us
# to declare them as C parameters in NRPy+. We will denote this with the _face infix/suffix.
alpha_face = gri.register_gridfunctions("AUXEVOL","alpha_face")
gamma_faceDD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL","gamma_faceDD","sym01")
beta_faceU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","beta_faceU")

# We'll need some more gridfunctions, now, to represent the reconstructions of BU and ValenciavU
# on the right and left faces
Valenciav_rU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_rU",DIM=3)
B_rU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","B_rU",DIM=3)
Valenciav_lU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_lU",DIM=3)
B_lU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","B_lU",DIM=3)

ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Stilde_flux_HLLED")

ixp.register_gridfunctions_for_single_rank1("AUXEVOL","PhievolParenU",DIM=3)
gri.register_gridfunctions("AUXEVOL","AevolParen")

# Declare this symbol
sqrt4pi = par.Cparameters("REAL",thismodule,"sqrt4pi","sqrt(4.0*M_PI)")

def add_to_Cfunction_dict__AD_gauge_term_psi6Phi_flux_term(includes=None):
    GRHD.compute_sqrtgammaDET(gammaDD)
    GRFFE.compute_AD_source_term_operand_for_FD(GRHD.sqrtgammaDET,betaU,alpha,psi6Phi,AD)
    GRFFE.compute_psi6Phi_rhs_flux_term_operand(gammaDD,GRHD.sqrtgammaDET,betaU,alpha,AD,psi6Phi)

    parens_to_print = [
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","AevolParen"),rhs=GRFFE.AevolParen),
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","PhievolParenU0"),rhs=GRFFE.PhievolParenU[0]),
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","PhievolParenU1"),rhs=GRFFE.PhievolParenU[1]),
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","PhievolParenU2"),rhs=GRFFE.PhievolParenU[2]),
                      ]

    desc = "Calculate quantities to be finite-differenced for the GRFFE RHSs"
    name = "calculate_AD_gauge_term_psi6Phi_flux_term_for_RHSs"
    params = "const paramstruct *restrict params,const REAL *restrict in_gfs,REAL *restrict auxevol_gfs"
    body = fin.FD_outputC("returnstring",parens_to_print,params=outCparams)
    loopopts = "AllPoints"
    rel_path_to_Cparams=os.path.join("../")
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        body=body, loopopts=loopopts)
#     return pickle_NRPy_env()

def add_to_Cfunction_dict__AD_gauge_term_psi6Phi_fin_diff(includes=None):
    xi_damping = par.Cparameters("REAL",thismodule,"xi_damping",0.1)
    GRFFE.compute_psi6Phi_rhs_damping_term(alpha,psi6Phi,xi_damping)

    AevolParen_dD = ixp.declarerank1("AevolParen_dD",DIM=3)
    PhievolParenU_dD = ixp.declarerank2("PhievolParenU_dD","nosym",DIM=3)

    A_rhsD = ixp.zerorank1()
    psi6Phi_rhs = GRFFE.psi6Phi_damping

    for i in range(3):
        A_rhsD[i] += -AevolParen_dD[i]
        psi6Phi_rhs += -PhievolParenU_dD[i][i]

    # Add Kreiss-Oliger dissipation to the GRFFE RHSs:
    # psi6Phi_dKOD = ixp.declarerank1("psi6Phi_dKOD")
    # AD_dKOD    = ixp.declarerank2("AD_dKOD","nosym")
    # for i in range(3):
    #     psi6Phi_rhs += diss_strength*psi6Phi_dKOD[i]*rfm.ReU[i] # ReU[i] = 1/scalefactor_orthog_funcform[i]
    #     for j in range(3):
    #         A_rhsD[j] += diss_strength*AD_dKOD[j][i]*rfm.ReU[i] # ReU[i] = 1/scalefactor_orthog_funcform[i]

    RHSs_to_print = [
                     lhrh(lhs=gri.gfaccess("rhs_gfs","AD0"),rhs=A_rhsD[0]),
                     lhrh(lhs=gri.gfaccess("rhs_gfs","AD1"),rhs=A_rhsD[1]),
                     lhrh(lhs=gri.gfaccess("rhs_gfs","AD2"),rhs=A_rhsD[2]),
                     lhrh(lhs=gri.gfaccess("rhs_gfs","psi6Phi"),rhs=psi6Phi_rhs),
                    ]

    desc = "Calculate AD gauge term and psi6Phi RHSs"
    name = "calculate_AD_gauge_psi6Phi_RHSs"
    params ="const paramstruct *params,const REAL *in_gfs,const REAL *auxevol_gfs,REAL *rhs_gfs"
    body     = fin.FD_outputC("returnstring",RHSs_to_print,params=outCparams)
    loopopts ="InteriorPoints"
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        body=body, loopopts=loopopts)
    outC_function_dict[name] = outC_function_dict[name].replace("= NGHOSTS","= NGHOSTS_A2B").replace("NGHOSTS+Nxx0","Nxx_plus_2NGHOSTS0-NGHOSTS_A2B").replace("NGHOSTS+Nxx1","Nxx_plus_2NGHOSTS1-NGHOSTS_A2B").replace("NGHOSTS+Nxx2","Nxx_plus_2NGHOSTS2-NGHOSTS_A2B")
    # Note the above .replace() functions. These serve to expand the loop range into the ghostzones, since
    # the second-order FD needs fewer than some other algorithms we use do.

import GiRaFFE_NRPy.GiRaFFE_NRPy_C2P_P2C as C2P_P2C
def add_to_Cfunction_dict__cons_to_prims(StildeD,BU,gammaDD,betaU,alpha,includes=None):
    C2P_P2C.GiRaFFE_NRPy_C2P(StildeD,BU,gammaDD,betaU,alpha)

    values_to_print = [
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD0"),rhs=C2P_P2C.outStildeD[0]),
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD1"),rhs=C2P_P2C.outStildeD[1]),
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD2"),rhs=C2P_P2C.outStildeD[2]),
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","ValenciavU0"),rhs=C2P_P2C.ValenciavU[0]),
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","ValenciavU1"),rhs=C2P_P2C.ValenciavU[1]),
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","ValenciavU2"),rhs=C2P_P2C.ValenciavU[2])
                      ]

    desc = "Apply fixes to \tilde{S}_i and recompute the velocity to match with current sheet prescription."
    name = "GiRaFFE_NRPy_cons_to_prims"
    params   ="const paramstruct *params,REAL *xx[3],REAL *auxevol_gfs,REAL *in_gfs"
    body     = fin.FD_outputC("returnstring",values_to_print,params=outCparams)
    loopopts ="AllPoints,Read_xxs"
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        body=body, loopopts=loopopts)

# TINYDOUBLE = par.Cparameters("REAL",thismodule,"TINYDOUBLE",1e-100)

def add_to_Cfunction_dict__prims_to_cons(gammaDD,betaU,alpha,  ValenciavU,BU, sqrt4pi,
                                         includes=None):
    C2P_P2C.GiRaFFE_NRPy_P2C(gammaDD,betaU,alpha,  ValenciavU,BU, sqrt4pi)

    values_to_print = [
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD0"),rhs=C2P_P2C.StildeD[0]),
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD1"),rhs=C2P_P2C.StildeD[1]),
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD2"),rhs=C2P_P2C.StildeD[2]),
                      ]

    desc = "Recompute StildeD after current sheet fix to Valencia 3-velocity to ensure consistency between conservative & primitive variables."
    name = "GiRaFFE_NRPy_prims_to_cons"
    params   ="const paramstruct *params,REAL *auxevol_gfs,REAL *in_gfs"
    body     = fin.FD_outputC("returnstring",values_to_print,params=outCparams)
    loopopts ="AllPoints"
    rel_path_to_Cparams=os.path.join("../")
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        body=body, loopopts=loopopts)

main_evolution_prototype = "void GiRaFFE_NRPy_RHSs(const paramstruct *restrict params,REAL *restrict auxevol_gfs,const REAL *restrict in_gfs,REAL *restrict rhs_gfs);"
main_evolution_func = """#include "NRPy_basic_defines.h"
#include "GiRaFFE_basic_defines.h"
#include "NRPy_function_prototypes.h"
void GiRaFFE_NRPy_RHSs(const paramstruct *restrict params,REAL *restrict auxevol_gfs,const REAL *restrict in_gfs,REAL *restrict rhs_gfs) {
#include "set_Cparameters.h"
    // First thing's first: initialize the RHSs to zero!
#pragma omp parallel for
    for(int ii=0;ii<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;ii++) {
        rhs_gfs[ii] = 0.0;
    }
    // Next calculate the easier source terms that don't require flux directions
    // This will also reset the RHSs for each gf at each new timestep.
    calculate_AD_gauge_term_psi6Phi_flux_term_for_RHSs(params,in_gfs,auxevol_gfs);
    calculate_AD_gauge_psi6Phi_RHSs(params,in_gfs,auxevol_gfs,rhs_gfs);

    // Now, we set up a bunch of structs of pointers to properly guide the PPM algorithm.
    // They also count the number of ghostzones available.
    gf_and_gz_struct in_prims[NUM_RECONSTRUCT_GFS], out_prims_r[NUM_RECONSTRUCT_GFS], out_prims_l[NUM_RECONSTRUCT_GFS];
    int which_prims_to_reconstruct[NUM_RECONSTRUCT_GFS],num_prims_to_reconstruct;
    const int Nxxp2NG012 = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

    REAL *temporary = auxevol_gfs + Nxxp2NG012*AEVOLPARENGF; //We're not using this anymore
    // This sets pointers to the portion of auxevol_gfs containing the relevant gridfunction.
    int ww=0;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAVU0GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU0GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU0GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAVU1GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU1GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU1GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAVU2GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU2GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU2GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BU0GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*B_RU0GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*B_LU0GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BU1GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*B_RU1GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*B_LU1GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BU2GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*B_RU2GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*B_LU2GF;
    ww++;

    // Prims are defined AT ALL GRIDPOINTS, so we set the # of ghostzones to zero:
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { in_prims[i].gz_lo[j]=0; in_prims[i].gz_hi[j]=0; }
    // Left/right variables are not yet defined, yet we set the # of gz's to zero by default:
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { out_prims_r[i].gz_lo[j]=0; out_prims_r[i].gz_hi[j]=0; }
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { out_prims_l[i].gz_lo[j]=0; out_prims_l[i].gz_hi[j]=0; }

    ww=0;
    which_prims_to_reconstruct[ww]=VX; ww++;
    which_prims_to_reconstruct[ww]=VY; ww++;
    which_prims_to_reconstruct[ww]=VZ; ww++;
    which_prims_to_reconstruct[ww]=BX; ww++;
    which_prims_to_reconstruct[ww]=BY; ww++;
    which_prims_to_reconstruct[ww]=BZ; ww++;
    num_prims_to_reconstruct=ww;

    // In each direction, perform the PPM reconstruction procedure.
    // Then, add the fluxes to the RHS as appropriate.
    for(int flux_dirn=0;flux_dirn<3;flux_dirn++) {
        // In each direction, interpolate the metric gfs (gamma,beta,alpha) to cell faces.
        interpolate_metric_gfs_to_cell_faces(params,auxevol_gfs,flux_dirn+1);
        // Then, reconstruct the primitive variables on the cell faces.
        // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE_NRPy.c"
        reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,
                                                which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
        // For example, if flux_dirn==0, then at gamma_faceDD00(i,j,k) represents gamma_{xx}
        // at (i-1/2,j,k), Valenciav_lU0(i,j,k) is the x-component of the velocity at (i-1/2-epsilon,j,k),
        // and Valenciav_rU0(i,j,k) is the x-component of the velocity at (i-1/2+epsilon,j,k).

        if(flux_dirn==0) {
            // Next, we calculate the source term for StildeD. Again, this also resets the rhs_gfs array at
            // each new timestep.
            calculate_StildeD0_source_term(params,auxevol_gfs,rhs_gfs);
            // Now, compute the electric field on each face of a cell and add it to the RHSs as appropriate
            //calculate_E_field_D0_right(params,auxevol_gfs,rhs_gfs);
            //calculate_E_field_D0_left(params,auxevol_gfs,rhs_gfs);
            // Finally, we calculate the flux of StildeD and add the appropriate finite-differences
            // to the RHSs.
            calculate_Stilde_flux_D0(params,auxevol_gfs,rhs_gfs);
        }
        else if(flux_dirn==1) {
            calculate_StildeD1_source_term(params,auxevol_gfs,rhs_gfs);
            //calculate_E_field_D1_right(params,auxevol_gfs,rhs_gfs);
            //calculate_E_field_D1_left(params,auxevol_gfs,rhs_gfs);
            calculate_Stilde_flux_D1(params,auxevol_gfs,rhs_gfs);
        }
        else {
            calculate_StildeD2_source_term(params,auxevol_gfs,rhs_gfs);
            //calculate_E_field_D2_right(params,auxevol_gfs,rhs_gfs);
            //calculate_E_field_D2_left(params,auxevol_gfs,rhs_gfs);
            calculate_Stilde_flux_D2(params,auxevol_gfs,rhs_gfs);
        }
        calculate_Stilde_rhsD(flux_dirn+1,params,auxevol_gfs,rhs_gfs);
        for(int count=0;count<=1;count++) {
            // This function is written to be general, using notation that matches the forward permutation added to AD2,
            // i.e., [F_HLL^x(B^y)]_z corresponding to flux_dirn=0, count=1.
            // The SIGN parameter is necessary because
            // -E_z(x_i,y_j,z_k) = 0.25 ( [F_HLL^x(B^y)]_z(i+1/2,j,k)+[F_HLL^x(B^y)]_z(i-1/2,j,k)
            //                           -[F_HLL^y(B^x)]_z(i,j+1/2,k)-[F_HLL^y(B^x)]_z(i,j-1/2,k) )
            // Note the negative signs on the reversed permutation terms!

            // By cyclically permuting with flux_dirn, we
            // get contributions to the other components, and by incrementing count, we get the backward permutations:
            // Let's suppose flux_dirn = 0. Then we will need to update Ay (count=0) and Az (count=1):
            //     flux_dirn=count=0 -> AD0GF+(flux_dirn+1+count)%3 = AD0GF + (0+1+0)%3=AD1GF <- Updating Ay!
            //        (flux_dirn)%3 = (0)%3 = 0               Vx
            //        (flux_dirn-count+2)%3 = (0-0+2)%3 = 2   Vz .  Inputs Vx, Vz -> SIGN = -1 ; 2.0*((REAL)count)-1.0=-1 check!
            //     flux_dirn=0,count=1 -> AD0GF+(flux_dirn+1+count)%3 = AD0GF + (0+1+1)%3=AD2GF <- Updating Az!
            //        (flux_dirn)%3 = (0)%3 = 0               Vx
            //        (flux_dirn-count+2)%3 = (0-1+2)%3 = 1   Vy .  Inputs Vx, Vy -> SIGN = +1 ; 2.0*((REAL)count)-1.0=2-1=+1 check!
            // Let's suppose flux_dirn = 1. Then we will need to update Az (count=0) and Ax (count=1):
            //     flux_dirn=1,count=0 -> AD0GF+(flux_dirn+1+count)%3 = AD0GF + (1+1+0)%3=AD2GF <- Updating Az!
            //        (flux_dirn)%3 = (1)%3 = 1               Vy
            //        (flux_dirn-count+2)%3 = (1-0+2)%3 = 0   Vx .  Inputs Vy, Vx -> SIGN = -1 ; 2.0*((REAL)count)-1.0=-1 check!
            //     flux_dirn=count=1 -> AD0GF+(flux_dirn+1+count)%3 = AD0GF + (1+1+1)%3=AD0GF <- Updating Ax!
            //        (flux_dirn)%3 = (1)%3 = 1               Vy
            //        (flux_dirn-count+2)%3 = (1-1+2)%3 = 2   Vz .  Inputs Vy, Vz -> SIGN = +1 ; 2.0*((REAL)count)-1.0=2-1=+1 check!
            // Let's suppose flux_dirn = 2. Then we will need to update Ax (count=0) and Ay (count=1):
            //     flux_dirn=2,count=0 -> AD0GF+(flux_dirn+1+count)%3 = AD0GF + (2+1+0)%3=AD0GF <- Updating Ax!
            //        (flux_dirn)%3 = (2)%3 = 2               Vz
            //        (flux_dirn-count+2)%3 = (2-0+2)%3 = 1   Vy .  Inputs Vz, Vy -> SIGN = -1 ; 2.0*((REAL)count)-1.0=-1 check!
            //     flux_dirn=2,count=1 -> AD0GF+(flux_dirn+1+count)%3 = AD0GF + (2+1+1)%3=AD1GF <- Updating Ay!
            //        (flux_dirn)%3 = (2)%3 = 2               Vz
            //        (flux_dirn-count+2)%3 = (2-1+2)%3 = 0   Vx .  Inputs Vz, Vx -> SIGN = +1 ; 2.0*((REAL)count)-1.0=2-1=+1 check!
            calculate_E_field_flat_all_in_one(params,
              &auxevol_gfs[IDX4ptS(VALENCIAV_RU0GF+(flux_dirn)%3, 0)],&auxevol_gfs[IDX4ptS(VALENCIAV_RU0GF+(flux_dirn-count+2)%3, 0)],
              &auxevol_gfs[IDX4ptS(VALENCIAV_LU0GF+(flux_dirn)%3, 0)],&auxevol_gfs[IDX4ptS(VALENCIAV_LU0GF+(flux_dirn-count+2)%3, 0)],
              &auxevol_gfs[IDX4ptS(B_RU0GF        +(flux_dirn)%3, 0)],&auxevol_gfs[IDX4ptS(B_RU0GF        +(flux_dirn-count+2)%3, 0)],
              &auxevol_gfs[IDX4ptS(B_LU0GF        +(flux_dirn)%3, 0)],&auxevol_gfs[IDX4ptS(B_LU0GF        +(flux_dirn-count+2)%3, 0)],
              &auxevol_gfs[IDX4ptS(B_RU0GF        +(flux_dirn-count+2)%3, 0)],
              &auxevol_gfs[IDX4ptS(B_LU0GF        +(flux_dirn-count+2)%3, 0)],
              &auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF,0)],&auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF,0)],&auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF,0)],
              &auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF,0)],&auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF,0)],&auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF,0)],
              &auxevol_gfs[IDX4ptS(BETA_FACEU0GF+flux_dirn,0)],&auxevol_gfs[IDX4ptS(BETA_FACEU0GF+(flux_dirn-count+2)%3,0)],&auxevol_gfs[IDX4ptS(ALPHA_FACEGF,0)],
              &rhs_gfs[IDX4ptS(AD0GF+(flux_dirn+1+count)%3,0)], 2.0*((REAL)count)-1.0, flux_dirn);
        }
    }
}"""

post_step_prototype = "void GiRaFFE_NRPy_post_step(const paramstruct *restrict params,REAL *xx[3],REAL *restrict auxevol_gfs,REAL *restrict evol_gfs,const int n);"
post_step_func = """#include "NRPy_basic_defines.h"
#include "GiRaFFE_basic_defines.h"
#include "NRPy_function_prototypes.h"
void GiRaFFE_NRPy_post_step(const paramstruct *restrict params,REAL *xx[3],REAL *restrict auxevol_gfs,REAL *restrict evol_gfs,const int n) {
    // First, apply BCs to AD and psi6Phi. Then calculate BU from AD
    apply_bcs_potential(params,evol_gfs);
    driver_A_to_B(params,evol_gfs,auxevol_gfs);
    //override_BU_with_old_GiRaFFE(params,auxevol_gfs,n);
    // Apply fixes to StildeD, then recompute the velocity at the new timestep.
    // Apply the current sheet prescription to the velocities
    GiRaFFE_NRPy_cons_to_prims(params,xx,auxevol_gfs,evol_gfs);
    // Then, recompute StildeD to be consistent with the new velocities
    //GiRaFFE_NRPy_prims_to_cons(params,auxevol_gfs,evol_gfs);
    // Finally, apply outflow boundary conditions to the velocities.
    apply_bcs_velocity(params,auxevol_gfs);
}"""

def add_to_Cfunction_dict__driver_function():
    outC_function_master_list.append(outC_function_element("empty", "empty", "empty", "empty", "GiRaFFE_NRPy_RHSs", "empty",
                             "empty", "empty", "empty", "empty",
                             "empty", "empty"))
    outC_function_outdir_dict["GiRaFFE_NRPy_RHSs"] = "default"
    outC_function_dict["GiRaFFE_NRPy_RHSs"] = main_evolution_func
    outC_function_prototype_dict["GiRaFFE_NRPy_RHSs"] = main_evolution_prototype

    outC_function_master_list.append(outC_function_element("empty", "empty", "empty", "empty", "GiRaFFE_NRPy_post_step", "empty",
                             "empty", "empty", "empty", "empty",
                             "empty", "empty"))
    outC_function_outdir_dict["GiRaFFE_NRPy_post_step"] = "default"
    outC_function_dict["GiRaFFE_NRPy_post_step"] = post_step_func
    outC_function_prototype_dict["GiRaFFE_NRPy_post_step"] = post_step_prototype