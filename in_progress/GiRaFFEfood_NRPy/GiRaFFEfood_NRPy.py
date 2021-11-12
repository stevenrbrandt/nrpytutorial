# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)
giraffefood_dir_path = os.path.join("GiRaFFEfood_NRPy")
if giraffefood_dir_path not in sys.path:
    sys.path.append(giraffefood_dir_path)

# Step 0: Import the NRPy+ core modules and set the reference metric to Cartesian
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import reference_metric as rfm   # NRPy+: Reference metric support
import GiRaFFEfood_NRPy_Common_Functions as gfcf # Some useful functions for GiRaFFE initial data.

par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()

# We import all ID modules ahead of time so that options can be changed *before* generating the functions.
import GiRaFFEfood_NRPy_Exact_Wald as gfew
import GiRaFFEfood_NRPy_Magnetospheric_Wald as gfmw
import GiRaFFEfood_NRPy_Split_Monopole as gfsm
import GiRaFFEfood_NRPy_Alfven_Wave as gfaw
import GiRaFFEfood_NRPy_Fast_Wave as gffw
import GiRaFFEfood_NRPy_Degen_Alfven_Wave as gfdaw
import GiRaFFEfood_NRPy_Three_Waves as gftw
import GiRaFFEfood_NRPy_FFE_Breakdown as gffb
import GiRaFFEfood_NRPy_Aligned_Rotator as gfar

# Step 1a: Set commonly used parameters.
thismodule = __name__

def GiRaFFEfood_NRPy_generate_initial_data(ID_type = "DegenAlfvenWave", stagger_enable = False,**params):
    global AD, ValenciavU
    if ID_type == "ExactWald":
        AD = gfcf.Axyz_func_spherical(gfew.Ar_EW,gfew.Ath_EW,gfew.Aph_EW,stagger_enable,**params)
        ValenciavU = gfew.ValenciavU_func_EW(**params)
    elif ID_type == "MagnetosphericWald":
        AD = gfcf.Axyz_func_spherical(gfmw.Ar_MW,gfmw.Ath_MW,gfmw.Aph_MW,stagger_enable,**params)
        ValenciavU = gfmw.ValenciavU_func_MW(**params)
    elif ID_type == "SplitMonopole":
        AD = gfcf.Axyz_func_spherical(gfsm.Ar_SM,gfsm.Ath_SM,gfsm.Aph_SM,stagger_enable,**params)
        ValenciavU = gfsm.ValenciavU_func_SM(**params)
    elif ID_type == "AlfvenWave":
        AD = gfcf.Axyz_func_Cartesian(gfaw.Ax_AW,gfaw.Ay_AW,gfaw.Az_AW, stagger_enable, **params)
        ValenciavU = gfaw.ValenciavU_func_AW(**params)
    elif ID_type == "FastWave":
        AD = gfcf.Axyz_func_Cartesian(gffw.Ax_FW,gffw.Ay_FW,gffw.Az_FW, stagger_enable, **params)
        ValenciavU = gffw.ValenciavU_func_FW(**params)
    elif ID_type == "DegenAlfvenWave":
        AD = gfcf.Axyz_func_Cartesian(gfdaw.Ax_DAW,gfdaw.Ay_DAW,gfdaw.Az_DAW, stagger_enable, **params)
        ValenciavU = gfdaw.ValenciavU_func_DAW(**params)
    elif ID_type == "ThreeWaves":
        AD = gfcf.Axyz_func_Cartesian(gftw.Ax_TW,gftw.Ay_TW,gftw.Az_TW, stagger_enable, **params)
        ValenciavU = gftw.ValenciavU_func_TW(**params)
    elif ID_type == "FFE_Breakdown":
        AD = gfcf.Axyz_func_Cartesian(gffb.Ax_FB,gffb.Ay_FB,gffb.Az_FB, stagger_enable, **params)
        ValenciavU = gffb.ValenciavU_func_FB(**params)
    elif ID_type == "AlignedRotator":
        AD = gfcf.Axyz_func_spherical(gfar.Ar_AR,gfar.Ath_AR,gfar.Aph_AR, stagger_enable, **params)
        ValenciavU = gfar.ValenciavU_func_AR(**params)