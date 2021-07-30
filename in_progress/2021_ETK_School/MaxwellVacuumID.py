#!/usr/bin/env python
# coding: utf-8

# <script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
# <script>
#   window.dataLayer = window.dataLayer || [];
#   function gtag(){dataLayer.push(arguments);}
#   gtag('js', new Date());
# 
#   gtag('config', 'UA-59152712-8');
# </script>
# 
# # `MaxwellVacuumID`: An Einstein Toolkit thorn for generating initial data for Maxwell's equations
# 
# ## Authors: Terrence Pierre Jacques, Patrick Nelson, & Zach Etienne
# ### Formatting improvements courtesy Brandon Clark
# 
# ### NRPy+ Source Code for this module: [Maxwell/InitialData.py](../edit/Maxwell/InitialData.py) [\[**tutorial**\]](Tutorial-VacuumMaxwell_InitialData.ipynb) Contructs the SymPy expressions for toroidal dipole field initial data
# 
# ## Introduction:
# In this part of the tutorial, we will construct an Einstein Toolkit (ETK) thorn (module) that will set up *initial data* for two formulations Maxwell's equations. In a [previous tutorial notebook](Tutorial-VacuumMaxwell_InitialData.ipynb), we used NRPy+ to contruct the SymPy expressions for toroidal dipole initial data. This thorn is largely based on and should function similarly to the NRPy+ generated [`IDScalarWaveNRPy`](Tutorial-ETK_thorn-IDScalarWaveNRPy.ipynb) thorn.
# 
# We will construct this thorn in two steps.
# 
# 1. Call on NRPy+ to convert the SymPy expressions for the initial data into one C-code kernel.
# 1. Write the C code and linkages to the Einstein Toolkit infrastructure (i.e., the .ccl files) to complete this Einstein Toolkit module.

# <a id='toc'></a>
# 
# # Table of Contents
# $$\label{toc}$$
# 
# This notebook is organized as follows
# 
# 1. [Step 1](#initializenrpy): Initialize needed Python/NRPy+ modules
# 1. [Step 2](#toroidal_id): NRPy+-generated C code kernels for toroidal dipole field initial data
# 1. [Step 3](#cclfiles): CCL files - Define how this module interacts and interfaces with the wider Einstein Toolkit infrastructure
#     1. [Step 3.a](#paramccl): `param.ccl`: specify free parameters within `MaxwellVacuumID`
#     1. [Step 3.b](#interfaceccl): `interface.ccl`: define needed gridfunctions; provide keywords denoting what this thorn provides and what it should inherit from other thorns
#     1. [Step 3.c](#scheduleccl): `schedule.ccl`:schedule all functions used within `MaxwellVacuumID`, specify data dependencies within said functions, and allocate memory for gridfunctions
# 1. [Step 4](#cdrivers): C driver functions for ETK registration & NRPy+-generated kernels
#     1. [Step 4.a](#etkfunctions): Initial data function
#     1. [Step 4.b](#makecodedefn): `make.code.defn`: List of all C driver functions needed to compile `MaxwellVacuumID`
# 1. [Step 5](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

# <a id='initializenrpy'></a>
# 
# # Step 1: Initialize needed Python/NRPy+ modules \[Back to [top](#toc)\]
# 
# $$\label{initializenrpy}$$

# In[2]:


# Step 1: Import needed core NRPy+ modules
import os,sys
base_dir      = os.getcwd()
nrpy_core_dir = os.path.join(base_dir,"..","..")
sys.path.append(nrpy_core_dir)

from outputC import lhrh         # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import loop as lp                # NRPy+: Generate C code loops
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import cmdline_helper as cmd  # NRPy+: Multi-platform Python command-line interface

# Step 1a: Create directories for the thorn if they don't exist.
# Create directory for MaxwellVacuumID thorn & subdirectories in case they don't exist.
outrootdir = "MaxwellVacuumID/"
cmd.mkdir(os.path.join(outrootdir))
outdir = os.path.join(outrootdir,"src") # Main C code output directory
cmd.mkdir(outdir)

# Step 1b: This is an Einstein Toolkit (ETK) thorn. Here we
#          tell NRPy+ that gridfunction memory access will
#          therefore be in the "ETK" style.
par.set_parval_from_str("grid::GridFuncMemAccess","ETK")


# <a id='toroidal_id'></a>
# 
# # Step 2: Constructing the Einstein Toolkit C-code calling functions that include the C code kernels \[Back to [top](#toc)\]
# $$\label{toroidal_id}$$
# 
# Using sympy, we construct the exact expressions for toroidal dipole field initial data currently supported in NRPy, documented in [Tutorial-VacuumMaxwell_InitialData.ipynb](Tutorial-VacuumMaxwell_InitialData.ipynb). We write the generated C codes into different C files, corresponding to the type of initial data the may want to choose at run time. Note that the code below can be easily extensible to include other types of initial data.

# In[3]:


import Maxwell.InitialData as mwid

# Set coordinate system. ETK only supports cartesian coordinates
CoordSystem     = "Cartesian"
par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)

# set up ID sympy expressions - System I
mwid.InitialData()

# x,y,z = gri.register_gridfunctions("AUX",["x","y","z"])

AIU = ixp.register_gridfunctions_for_single_rank1("EVOL","AU")
EIU = ixp.register_gridfunctions_for_single_rank1("EVOL","EU")
psiI = gri.register_gridfunctions("EVOL","Phi")

# Set which system to use, which are defined in Maxwell/VacuumMaxwell_Flat_Cartesian_ID.py
par.set_parval_from_str("Maxwell.InitialData::System_to_use","System_II")

# set up ID sympy expressions - System II
mwid.InitialData()

# AIIU = ixp.register_gridfunctions_for_single_rank1("EVOL","AIIU")
# EIIU = ixp.register_gridfunctions_for_single_rank1("EVOL","EIIU")
# psiII = gri.register_gridfunctions("EVOL","psiII")
GammaII = gri.register_gridfunctions("EVOL","Gamma")

Maxwell_ID_SymbExpressions = [                 lhrh(lhs=gri.gfaccess("out_gfs","AU0"),rhs=mwid.AidU[0]),                 lhrh(lhs=gri.gfaccess("out_gfs","AU1"),rhs=mwid.AidU[1]),                 lhrh(lhs=gri.gfaccess("out_gfs","AU2"),rhs=mwid.AidU[2]),                 lhrh(lhs=gri.gfaccess("out_gfs","EU0"),rhs=mwid.EidU[0]),                 lhrh(lhs=gri.gfaccess("out_gfs","EU1"),rhs=mwid.EidU[1]),                 lhrh(lhs=gri.gfaccess("out_gfs","EU2"),rhs=mwid.EidU[2]),                 lhrh(lhs=gri.gfaccess("out_gfs","Phi"),rhs=mwid.psi_ID),#                  lhrh(lhs=gri.gfaccess("out_gfs","AIIU0"),rhs=mwid.AidU[0]),\
#                  lhrh(lhs=gri.gfaccess("out_gfs","AIIU1"),rhs=mwid.AidU[1]),\
#                  lhrh(lhs=gri.gfaccess("out_gfs","AIIU2"),rhs=mwid.AidU[2]),\
#                  lhrh(lhs=gri.gfaccess("out_gfs","EIIU0"),rhs=mwid.EidU[0]),\
#                  lhrh(lhs=gri.gfaccess("out_gfs","EIIU1"),rhs=mwid.EidU[1]),\
#                  lhrh(lhs=gri.gfaccess("out_gfs","EIIU2"),rhs=mwid.EidU[2]),\
#                  lhrh(lhs=gri.gfaccess("out_gfs","psiII"),rhs=mwid.psi_ID),\
                 lhrh(lhs=gri.gfaccess("out_gfs","Gamma"),rhs=mwid.Gamma_ID)]
declare_string = """
const double x = xGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
const double y = yGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
const double z = zGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];

"""
Maxwell_ID_CcodeKernel = fin.FD_outputC("returnstring",
                        Maxwell_ID_SymbExpressions,\
                        params="outCverbose=True")

Maxwell_ID_looped = lp.loop(["i2","i1","i0"],["0","0","0"],["cctk_lsh[2]","cctk_lsh[1]","cctk_lsh[0]"],                               ["1","1","1"],["#pragma omp parallel for","",""],"",                                declare_string+Maxwell_ID_CcodeKernel).replace("time","cctk_time")                                                                      .replace("xx0", "x")                                                                      .replace("xx1", "y")                                                                      .replace("xx2", "z")


# Step 4: Write the C code kernel to file.
with open(os.path.join(outdir,"Maxwell_ID.h"), "w") as file:
    file.write(str(Maxwell_ID_looped))


# <a id='cclfiles'></a>
# 
# # Step 3: ETK `ccl` file generation \[Back to [top](#toc)\]
# $$\label{cclfiles}$$
# 
# <a id='paramccl'></a>
# 
# ## Step 3.a: `param.ccl`: specify free parameters within `MaxwellVacuumID` \[Back to [top](#toc)\]
# $$\label{paramccl}$$
# 
# All parameters necessary for the computation of the initial data expressions are registered within NRPy+; we use this information to automatically generate `param.ccl`. NRPy+ also specifies default values for each parameter. 
# 
# More information on `param.ccl` syntax can be found in the [official Einstein Toolkit documentation](https://einsteintoolkit.org/usersguide/UsersGuide.html#x1-184000D2.3).

# In[4]:


def keep_param__return_type(paramtuple):
    keep_param = True # We'll not set some parameters in param.ccl;
                      #   e.g., those that should be #define'd like M_PI.
    typestring = ""
    # Separate thorns within the ETK take care of grid/coordinate parameters;
    #   thus we ignore NRPy+ grid/coordinate parameters:
    if paramtuple.module == "grid" or paramtuple.module == "reference_metric":
        keep_param = False

    partype = paramtuple.type
    if partype == "bool":
        typestring += "BOOLEAN "
    elif partype == "REAL":
        if paramtuple.defaultval != 1e300: # 1e300 is a magic value indicating that the C parameter should be mutable
            typestring += "CCTK_REAL "
        else:
            keep_param = False
    elif partype == "int":
        typestring += "CCTK_INT "
    elif partype == "#define":
        keep_param = False
    elif partype == "char":
        print("Error: parameter "+paramtuple.module+"::"+paramtuple.parname+
              " has unsupported type: \""+ paramtuple.type + "\"")
        sys.exit(1)
    else:
        print("Error: parameter "+paramtuple.module+"::"+paramtuple.parname+
              " has unsupported type: \""+ paramtuple.type + "\"")
        sys.exit(1)
    return keep_param, typestring


with open(os.path.join(outrootdir,"param.ccl"), "w") as file:
    file.write("""
# This param.ccl file was automatically generated by NRPy+.
#   You are advised against modifying it directly; instead
#   modify the Python code that generates it.

shares: grid

USES KEYWORD type

CCTK_KEYWORD initial_data "Type of initial data"
{
  "toroid"      :: "Toroidal Dipole field"
} "toroid"

restricted:

""")

    paramccl_str = ""
    for i in range(len(par.glb_Cparams_list)):
        # keep_param is a boolean indicating whether we should accept or reject
        #    the parameter. singleparstring will contain the string indicating
        #    the variable type.
        keep_param, singleparstring = keep_param__return_type(par.glb_Cparams_list[i])

        if keep_param:
            parname = par.glb_Cparams_list[i].parname
            partype = par.glb_Cparams_list[i].type
            singleparstring += parname + " \""+ parname +" (see NRPy+ for parameter definition)\"\n"
            singleparstring += "{\n"
            if partype != "bool":
                singleparstring += " *:* :: \"All values accepted. NRPy+ does not restrict the allowed ranges of parameters yet.\"\n"
            singleparstring += "} "+str(par.glb_Cparams_list[i].defaultval)+"\n\n"

            paramccl_str += singleparstring
    file.write(paramccl_str)


# <a id='interfaceccl'></a>
# 
# ## Step 3.b: `interface.ccl`: define needed gridfunctions; provide keywords denoting what this thorn provides and what it should inherit from other thorns \[Back to [top](#toc)\]
# $$\label{interfaceccl}$$
# 
# `interface.ccl` declares all gridfunctions and determines how `MaxwellVacuumID` interacts with other Einstein Toolkit thorns.
# 
# The [official Einstein Toolkit documentation](https://einsteintoolkit.org/usersguide/UsersGuide.html#x1-179000D2.2) defines what must/should be included in an `interface.ccl` file. 

# In[5]:


evol_gfs_list    = []
for i in range(len(gri.glb_gridfcs_list)):
    if gri.glb_gridfcs_list[i].gftype == "EVOL":
        evol_gfs_list.append(   gri.glb_gridfcs_list[i].name+"GF")

# NRPy+'s finite-difference code generator assumes gridfunctions
#    are alphabetized; not sorting may result in unnecessary
#    cache misses.
evol_gfs_list.sort()

with open(os.path.join(outrootdir,"interface.ccl"), "w") as file:
    file.write("""
# With "implements", we give our thorn its unique name.
implements: MaxwellVacuumID

# By "inheriting" other thorns, we tell the Toolkit that we
#   will rely on variables/function that exist within those
#   functions.
inherits: MaxwellVacuum grid
""")


# <a id='scheduleccl'></a>
# 
# ## Step 3.c: `schedule.ccl`: schedule all functions used within `MaxwellVacuumID`, specify data dependencies within said functions, and allocate memory for gridfunctions \[Back to [top](#toc)\]
# $$\label{scheduleccl}$$
# 
# Official documentation on constructing ETK `schedule.ccl` files is found [here](https://einsteintoolkit.org/usersguide/UsersGuide.html#x1-187000D2.4). 

# In[6]:


with open(os.path.join(outrootdir,"schedule.ccl"), "w") as file:
    file.write("""
# This schedule.ccl file was automatically generated by NRPy+.
#   You are advised against modifying it directly; instead
#   modify the Python code that generates it.

schedule Maxwell_InitialData at CCTK_INITIAL as Maxwell_InitialData
{
  STORAGE: MaxwellVacuum::evol_variables[3]
  LANG:          C
} "Initial data for Maxwell's equations"

""")


# <a id='cdrivers'></a>
# 
# # Step 4: C driver functions for ETK registration & NRPy+-generated kernels \[Back to [top](#toc)\]
# $$\label{cdrivers}$$
# 
# Now that we have constructed the basic C code kernels and the needed Einstein Toolkit `ccl` files, we next write the driver functions for registering `MaxwellVacuumID` within the Toolkit and the C code kernels. Each of these driver functions is called directly from [`schedule.ccl`](#scheduleccl).

# In[7]:


make_code_defn_list = []
def append_to_make_code_defn_list(filename):
    if filename not in make_code_defn_list:
        make_code_defn_list.append(filename)
    return os.path.join(outdir,filename)


# <a id='etkfunctions'></a>
# 
# ## Step 4.a: Initial data function \[Back to [top](#toc)\]
# $$\label{etkfunctions}$$
# 
# Here we define the initial data function, and how it's to be called in the `schedule.ccl` file by ETK.

# In[8]:


with open(append_to_make_code_defn_list("InitialData.c"),"w") as file:
    file.write("""

#include <math.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

void Maxwell_InitialData(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  const CCTK_REAL *xGF = x;
  const CCTK_REAL *yGF = y;
  const CCTK_REAL *zGF = z;
#include "Maxwell_ID.h"
}
""")


# <a id='makecodedefn'></a>
# 
# ## Step 4.b: `make.code.defn`: List of all C driver functions needed to compile `MaxwellVacuumID` \[Back to [top](#toc)\]
# $$\label{makecodedefn}$$
# 
# When constructing each C code driver function above, we called the `append_to_make_code_defn_list()` function, which built a list of each C code driver file. We'll now add each of those files to the `make.code.defn` file, used by the Einstein Toolkit's build system.

# In[9]:


with open(os.path.join(outdir,"make.code.defn"), "w") as file:
    file.write("""
# Main make.code.defn file for thorn MaxwellVacuumID

# Source files in this directory
SRCS =""")
    filestring = ""
    for i in range(len(make_code_defn_list)):
        filestring += "      "+make_code_defn_list[i]
        if i != len(make_code_defn_list)-1:
            filestring += " \\\n"
        else:
            filestring += "\n"
    file.write(filestring)


# <a id='latex_pdf_output'></a>
# 
# # Step 5: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
# $$\label{latex_pdf_output}$$
# 
# The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
# [Tutorial-ETK_thorn-MaxwellVacuumID.pdf](Tutorial-ETK_thorn-MaxwellVacuumID.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)

# In[9]:


import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ETK_thorn-MaxwellVacuumID")

