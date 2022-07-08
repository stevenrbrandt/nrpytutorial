import grid
import os
from sympy import sympify
from sympy import pi
from sympy.functions.elementary.miscellaneous import cbrt, Max

import NRPy_param_funcs as par
import indexedexp as ixp
from cactusthorn import CactusThorn, loop
from outputC import lhrh

# Current options are Carpet and CarpetX
grid.ET_driver = "CarpetX"

thorn = CactusThorn("CarpetXNRPy", "Z4cNRPy")

fd_order = thorn.declare_param('fd_order', default=4, vmin=2, vmax=8, doc="Finite differencing order")

fd_order = 4
par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", fd_order)

centering='VVV'

# EVOL: evolved grid functions (possibly multiple time levels)
# AUXEVOL: needed for evolution, can be freed after evaluating rhs (1 time level)
# AUX: e.g. RHS (1 time level)
# TMP: not actually a grid function, only temporary tile
import sympy as sp

# Generic helpers

def flatten(lists):
    return sum(lists, [])

def sum1(expr):
    result = sympify(0)
    for i in range(3):
        result += expr(i)
    return result

def sum2(expr):
    result = sympify(0)
    for i in range(3):
        for j in range(3):
            result += expr(i, j)
    return result

def sum2_symm(expr):
    result = sympify(0)
    for i in range(3):
        for j in range(i+1):
            result += (1 if i==j else 2) * expr(i, j)
    return result

def tensor1(expr):
    return [expr(i) for i in range(3)]

def tensor2(expr):
    return [[expr(i, j) for j in range(3)] for i in range(3)]

def tensor3(expr):
    return [[[expr(i, j, k) for k in range(3)] for j in range(3)] for i in range(3)]

# TODO: Specify names `gxx` etc.
def ixnam(i):
    return ["x","y","z"][i]

def namefun(symbol, index, shape, prefix):
    symbol = prefix
    result = [sp.Symbol(symbol + ''.join(ixnam(n) for n in index + [i]))
            if symbol else sp.sympify(0) for i in range(shape[0])]
    return result

d = ["x","y","z"]
for i in range(3):
    for j in range(i,3):
        n1 = f"metric{i}{j}"
        n2 = f"g{d[i]}{d[j]}"
        grid.rename(n1,n2)

# TODO: turn these into a Cactus parameters
kappa1 = thorn.declare_param('kappa1', default=0.02, doc="kappa1")
kappa2 = thorn.declare_param('kappa2', default=0.0, doc="kappa2")
f_mu_L = thorn.declare_param('f_mu_L', default=2.0, doc="f_mu_L")
f_mu_S = thorn.declare_param('f_mu_S', default=1.0, doc="f_mu_S")
eta = thorn.declare_param('eta', default=2.0, doc="eta")
alphaG_floor = thorn.declare_param('alphaG_floor', default=1.0e-10, vmin=0.0, doc="alphaG_floor")
chi_floor = thorn.declare_param('chi_floor', default=1.0e-10, vmin=0.0, doc="chi_floor")

gDD = ixp.register_gridfunctions_for_single_rankN(2,"EXTERNAL", "metric", "sym01", centering=centering, external_module="ADMBase", namefun=lambda sym,ind,shape: namefun(sym,ind,shape,"g"))
KDD = ixp.register_gridfunctions_for_single_rankN(2,"EXTERNAL", "extcurv", "sym01", centering=centering, external_module="ADMBase", namefun=lambda sym,ind,shape: namefun(sym,ind,shape,"k"))
alp = thorn.register_gridfunctions("EXTERNAL", ["alp"], centering=centering, external_module="ADMBase")
dtalp = thorn.register_gridfunctions("EXTERNAL", ["dtalp"], centering=centering, external_module="ADMBase")
betaU = ixp.register_gridfunctions_for_single_rankN(1,"EXTERNAL", "shift", centering=centering, external_module="ADMBase", namefun=lambda sym,ind,shape: namefun(sym,ind,shape,"beta"))
dtbetaU = ixp.register_gridfunctions_for_single_rankN(1,"EXTERNAL", "dtshift", centering=centering, external_module="ADMBase", namefun=lambda sym,ind,shape: namefun(sym,ind,shape,"dtbeta"))

eTtt = thorn.register_gridfunctions("EXTERNAL", ["eTtt"], centering=centering, external_module="TmunuBase")
eTtiD = ixp.register_gridfunctions_for_single_rankN(1,"EXTERNAL", "eTti", centering=centering, external_module="TmunuBase", namefun=lambda sym,ind,shape: namefun(sym,ind,shape,"eTt"))
eTijDD = ixp.register_gridfunctions_for_single_rankN(2,"EXTERNAL", "eTtij", "sym01", centering=centering, external_module="TmunuBase", namefun=lambda sym,ind,shape: namefun(sym,ind,shape,"eT"))

chi = thorn.register_gridfunctions("EVOL", ["chi"], centering=centering)
gammatildeDD = ixp.register_gridfunctions_for_single_rankN(2,"EVOL", "gammatildeDD", "sym01", centering=centering)
Khat = thorn.register_gridfunctions("EVOL", ["Khat"], centering=centering)
AtildeDD = ixp.register_gridfunctions_for_single_rankN(2,"EVOL", "AtildeDD", "sym01", centering=centering)
GammatildeU = ixp.register_gridfunctions_for_single_rankN(1,"EVOL", "GammatildeU", centering=centering)
Theta = thorn.register_gridfunctions("EVOL", ["Theta"], centering=centering)
alphaG = thorn.register_gridfunctions("EVOL", ["alphaG"], centering=centering)
betaGU = ixp.register_gridfunctions_for_single_rankN(1,"EVOL", "betaGU", centering=centering)

chi_rhs = thorn.register_gridfunctions("AUX", ["chi_rhs"], centering=centering)
gammatildeDD_rhs = ixp.register_gridfunctions_for_single_rankN(2,"AUX", "gammatildeDD_rhs", "sym01", centering=centering)
Khat_rhs = thorn.register_gridfunctions("AUX", ["Khat_rhs"], centering=centering)
AtildeDD_rhs = ixp.register_gridfunctions_for_single_rankN(2,"AUX", "AtildeDD_rhs", "sym01", centering=centering)
GammatildeU_rhs = ixp.register_gridfunctions_for_single_rankN(1,"AUX", "GammatildeU_rhs", centering=centering)
Theta_rhs = thorn.register_gridfunctions("AUX", ["Theta_rhs"], centering=centering)
alphaG_rhs = thorn.register_gridfunctions("AUX", ["alphaG_rhs"], centering=centering)
betaGU_rhs = ixp.register_gridfunctions_for_single_rankN(1,"AUX", "betaGU_rhs", centering=centering)

dchiD = ixp.register_gridfunctions_for_single_rankN(1,"TMP", "dchiD", centering=centering)
dgammatildeDDD = ixp.register_gridfunctions_for_single_rankN(3, "TMP", "dgammatildeDDD", "sym01", centering=centering)
dalphaGD = ixp.register_gridfunctions_for_single_rankN(1,"TMP", "dalphaGD", centering=centering)
ddalphaGDD = ixp.register_gridfunctions_for_single_rankN(2, "TMP", "ddalphaGDD", "sym01", centering=centering)
dbetaGUD = ixp.register_gridfunctions_for_single_rankN(2, "TMP", "dbetaGUD", "", centering=centering)

chi_dD = ixp.declarerank1("chi_dD")
gammatildeDD_dD = ixp.declarerank3("gammatildeDD_dD", "sym01")
alphaG_dD = ixp.declarerank1("alphaG_dD")
alphaG_dDD = ixp.declarerank2("alphaG_dDD", "sym01")
betaGU_dD = ixp.declarerank2("betaGU_dD", "")

# Expressions as calculated from the ADM variables

gUU, detg = ixp.symm_matrix_inverter3x3(gDD)
cbrt_detg = cbrt(detg)
trK = sum2(lambda i, j: gUU[i][j] * KDD[i][j])

# Expressions as calculated from the Z4c variables

gammatildeUU, detgammatilde = ixp.symm_matrix_inverter3x3(gammatildeDD)
cbrt_detgammatilde = cbrt(detgammatilde)
trAtilde = sum2(lambda i, j: gammatildeUU[i][j] * AtildeDD[i][j])

g1DD = tensor2(lambda i, j: 1 / chi * gammatildeDD[i][j])
g1UU = tensor2(lambda i, j: chi * gammatildeUU[i][j])
dg1DDD = tensor3(lambda i, j, k: -dchiD[k] / chi**2 * gammatildeDD[i][j] + 1 / chi * dgammatildeDDD[i][k][k])

# Initial conditions

chi1 = 1 / cbrt(detg)
Theta1 = sympify(0)

initial1_eqns = flatten([
    [lhrh(lhs=chi, rhs=chi1)],
    [lhrh(lhs=gammatildeDD[i][j], rhs=chi1 * gDD[i][j]) for i in range(3) for j in range(i+1)],
    [lhrh(lhs=Theta, rhs=Theta1)],
    [lhrh(lhs=Khat, rhs=trK - 2 * Theta1)],
    [lhrh(lhs=AtildeDD[i][j], rhs=chi1 * (KDD[i][j] - trK / 3 * gDD[i][j])) for i in range(3) for j in range(i+1)],
    [lhrh(lhs=alphaG, rhs=alp)],
    [lhrh(lhs=betaGU[i], rhs=betaU[i]) for i in range(3)],
])

thorn.add_func("Z4cNRPy_Initial1",
    body=initial1_eqns,
    where='everywhere',
    schedule_bin="Z4cNRPy_InitialGroup",
    doc="Convert ADM to Z4c variables, part 1",
    centering=centering)

initial2_eqns = flatten([
    flatten([
        flatten([
            [lhrh(lhs=dgammatildeDDD[i][j][k], rhs=gammatildeDD_dD[i][j][k]) for k in range(3)],
            [loop],
        ])
        for i in range(3) for j in range(i+1)
    ]),
    [lhrh(lhs=GammatildeU[i], rhs=sum2(lambda j, k: gammatildeUU[j][k] * dgammatildeDDD[j][k][i])) for i in range(3)],
])

thorn.add_func("Z4cNRPy_Initial2",
    body=initial2_eqns,
    where='interior',
    schedule_bin="Z4cNRPy_InitialGroup AFTER Z4cNRPy_Initial1",
    doc="Convert ADM to Z4c variables, part 2",
    centering=centering)

# Enforce constaints

enforce_eqns = flatten([
    # Enforce floors
    [lhrh(lhs=chi, rhs=Max(chi_floor, chi))],
    [lhrh(lhs=alphaG, rhs=Max(alphaG_floor, alphaG))],
    # Enforce algebraic constraints; see arXiv:1212.2901 [gr-qc]
    [lhrh(lhs=gammatildeDD[i][j], rhs=(1 / cbrt_detgammatilde) * gammatildeDD[i][j])
     for i in range(3) for j in range(i+1)],
    [lhrh(lhs=AtildeDD[i][j], rhs=(AtildeDD[i][j] - trAtilde / 3 * gammatildeDD[i][j]))
     for i in range(3) for j in range(i+1)],
])

thorn.add_func("Z4cNRPy_Enforce",
    body=enforce_eqns,
    where='interior',
    schedule_bin="Z4cNRPy_PostStepGroup",
    doc="Enforce constaints",
    sync="chiGF gammatildeDDGF KhatGF AtildeDDGF GammatildeUGF ThetaGF alphaGGF betaGUGF",
    centering=centering
    )

# Calculate ADM variables

adm_eqns = flatten([
    [lhrh(lhs=gDD[i][j], rhs=g1DD[i][j]) for i in range(3) for j in range(i+1)],
    [lhrh(lhs=KDD[i][j], rhs=1 / chi * (AtildeDD[i][j] + (Khat + 2 * Theta) / 3 * gammatildeDD[i][j]))
     for i in range(3) for j in range(i+1)],
    [lhrh(lhs=alp, rhs=alphaG)],
    [lhrh(lhs=dtalp, rhs=-alphaG * f_mu_L * Khat)],
    [lhrh(lhs=betaU[i], rhs=betaGU[i]) for i in range(3)],
    [lhrh(lhs=dtbetaU[i], rhs=f_mu_S * GammatildeU[i] - eta * betaGU[i]) for i in range(3)],
])

thorn.add_func("Z4cNRPy_ADM",
    body=adm_eqns,
    where='interior',
    schedule_bin="Z4cNRPy_PostStepGroup AFTER Z4cNRPy_Enforce",
    doc="Calculate ADM variables",
    centering=centering)

# Calculate RHS

AtildeUD = tensor2(lambda i, j: sum1(lambda x: gammatildeUU[i][x] * AtildeDD[x][j]))
GammaDDD = tensor3(lambda i, j, k: 1/2 * dg1DDD[i][j][k] + dg1DDD[i][k][j] - dg1DDD[j][k][i])
GammaUDD = tensor3(lambda i, j, k: sum1(lambda x: gUU[i][x] * GammaDDD[x][j][k]))
DDalphaGDD = tensor2(lambda i, j: ddalphaGDD[i][j] - sum1(lambda x: GammaUDD[x][i][j] * dalphaGD[x]))

rho = 1 / alphaG**2 * (eTtt
                       - 2 * sum1(lambda x: betaGU[x] * eTtiD[x])
                       + sum1(lambda x: betaGU[x] * sum1(lambda y: betaGU[y] * eTijDD[x][y])))
SiD = tensor1(lambda i: -1 / alphaG * (eTtiD[i]
                                       - sum1(lambda x:  betaGU[x] * eTijDD[i][x])))
SijDD = tensor2(lambda i, j: eTijDD[i][j])
SijUD = tensor2(lambda i, j: sum1(lambda x: gUU[i][x] * SijDD[x][j]))
traceSij = sum1(lambda x: SijUD[x][x])

rhs_eqns = flatten([
    flatten([
        flatten([
            [lhrh(lhs=dbetaGUD[i][j], rhs=betaGU_dD[i][j]) for j in range(3)],
            [loop],
        ])
        for i in range(3)
    ]),
    [lhrh(lhs=chi_rhs,
          rhs=2 / 3 * chi * (alphaG * (Khat + 2 * Theta) - sum1(lambda x: dbetaGUD[x][x])))],
    [lhrh(lhs=gammatildeDD_rhs[i][j],
          rhs=(-2 * alphaG * AtildeDD[i][j]
               + sum1(lambda x: gammatildeDD[x][i] * dbetaGUD[x][j] + gammatildeDD[x][j] * dbetaGUD[x][i])
               - 2 / 3 * sum1(lambda x: gammatildeDD[i][j] * dbetaGUD[x][x])))
     for i in range(3) for j in range(i+1)],
    [lhrh(lhs=Khat_rhs,
          rhs=(- sum2_symm(lambda x, y: g1UU[x][y] * DDalphaGDD[x][y])
               + alphaG * (sum2_symm(lambda x, y: AtildeUD[x][y] * AtildeUD[y][x])
                           + 1 / 3 * (Khat + 2 * Theta)**2)
               + 4 * pi * alphaG * (traceSij + rho)
               + alphaG * kappa1 * (1 - kappa2) * Theta))],
])

thorn.add_func("Z4cNRPy_RHS",
    body=rhs_eqns,
    where='interior',
    schedule_bin="ODESolvers_RHS",
    doc="Calculate RHS",
    centering=centering)

# Generate thorn

assert "CACTUS_HOME" in os.environ, "Please set the CACTUS_HOME variable to point to your Cactus installation"
cactus_home = os.environ["CACTUS_HOME"]
cactus_sim = os.environ.get("CACTUS_SIM", "sim")
cactus_thornlist = os.environ.get("CACTUS_THORNLIST", None)

schedule_raw = """
# Define schedule groups

SCHEDULE GROUP Z4c_InitialGroup AT initial AFTER ADMBase_PostInitial
{
} "Convert ADM to Z4c variables"

SCHEDULE GROUP Z4cNRPy_PostStepGroup AT initial AFTER Z4cNRPy_InitialGroup BEFORE ADMBase_SetADMVars
{
} "Post-process Z4c variables"

SCHEDULE GROUP Z4cNRPy_PostStepGroup AT postregrid BEFORE ADMBase_SetADMVars
{
} "Post-process Z4c variables"

SCHEDULE GROUP Z4cNRPy_PostStepGroup IN ODESolvers_PostStep BEFORE ADMBase_SetADMVars
{
} "Post-process Z4c variables"
"""

thorn.generate(cactus_home, cactus_config=cactus_sim, cactus_thornlist=cactus_thornlist,schedule_raw=schedule_raw)
