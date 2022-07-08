import grid
import os
from sympy import sympify
from sympy.functions.elementary.miscellaneous import cbrt

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

# TODO: Specify names `gxx` etc.
def ixnam(i):
    return ["x","y","z"][i]

def namefun(symbol, index, shape, prefix):
    symbol = prefix
    result = [sp.Symbol(symbol + ''.join(ixnam(n) for n in index + [i]))
            if symbol else sp.sympify(0) for i in range(shape[0])]
    return result

gDD = ixp.register_gridfunctions_for_single_rankN(2,"EXTERNAL", "metric", "sym01", centering=centering, external_module="ADMBase",namefun=lambda sym,ind,shape: namefun(sym,ind,shape,"g"))

d = ["x","y","z"]
for i in range(3):
    for j in range(i,3):
        n1 = f"metric{i}{j}"
        n2 = f"g{d[i]}{d[j]}"
        grid.rename(n1,n2)

KDD = ixp.register_gridfunctions_for_single_rankN(2,"EXTERNAL", "extcurv", "sym01", centering=centering, external_module="ADMBase",namefun=lambda sym,ind,shape: namefun(sym,ind,shape,"k"))
alp = thorn.register_gridfunctions("EXTERNAL", ["alp"], centering=centering, external_module="ADMBase")
betaU = ixp.register_gridfunctions_for_single_rankN(1,"EXTERNAL", "shift", centering=centering, external_module="ADMBase",namefun=lambda sym,ind,shape: namefun(sym,ind,shape,"beta"))

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

dchi = ixp.register_gridfunctions_for_single_rankN(1,"TMP", "dchi", centering=centering)
dgammatildeDDD = ixp.register_gridfunctions_for_single_rankN(3,"TMP", "dgammatildeDDD", "sym01", centering=centering)

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

gUU, detg = ixp.symm_matrix_inverter3x3(gDD)
trK = sum2(lambda i, j: gUU[i][j] * KDD[i][j])
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

# access a variable with a different centering using interpolation
# looping cell-centered, access vertex-centered, not vice-versa
# all rhs variables should have the same centering
# wave toy with fluxes, fluxes are faces
# schedule something in post-regrid, apply bc's
thorn.add_func("Z4c_initial1",
    body=initial1_eqns,
    where='everywhere',
    schedule_bin="initial AFTER ADMBase_PostInitial",
    doc="Convert ADM to Z4c variables, part 1",
    centering=centering)

chi_dD = ixp.declarerank1("chi_dD")
gammatildeDD_dD = ixp.declarerank3("gammatildeDD_dD", "sym01")

gammatildeUU, detgammatilde = ixp.symm_matrix_inverter3x3(gammatildeDD)

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

thorn.add_func("Z4c_initial2",
    body=initial2_eqns,
    where='interior',
    schedule_bin="initial AFTER Z4c_initial1",
    doc="Convert ADM to Z4c variables, part 2",
    centering=centering,
    sync="ADMBase::metric # added")

assert "CACTUS_HOME" in os.environ, "Please set the CACTUS_HOME variable to point to your Cactus installation"
cactus_home = os.environ["CACTUS_HOME"]
cactus_sim = os.environ.get("CACTUS_SIM", "sim")
cactus_thornlist = os.environ.get("CACTUS_THORNLIST", None)

thorn.generate(cactus_home, cactus_config=cactus_sim, cactus_thornlist=cactus_thornlist,schedule_raw="# RAW")
