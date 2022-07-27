import grid
import os
from sympy import sympify
from sympy import Eq, Ne, Piecewise, Rational
from sympy import pi
from sympy.functions.elementary.miscellaneous import cbrt, Max, sqrt

import NRPy_param_funcs as par
import indexedexp as ixp
from cactusthorn import CactusThorn, loop
from outputC import lhrh, outCparams

outCparams.CSE_enable = "False"

# "CCTK_REAL" or "CCTK_REALVEC"
par.set_parval_from_str("PRECISION", "CCTK_REALVEC")

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
# TILE_TMP: not actually a grid function, only temporary tile
# SCALAR_TMP: not actually a grid function, only a double
import sympy as sp

# Generic helpers

def flatten(lists):
    new_list = []
    for item in lists:
        if type(item) == list:
            new_list += flatten(item)
        else:
            new_list += [item]
    return new_list

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

def tensor4(expr):
    return [[[[expr(i, j, k, l) for l in range(3)] for k in range(3)] for j in range(3)] for i in range(3)]

# TODO: Specify names `gxx` etc.
def ixnam(i):
    return ["x", "y", "z"][i]

def namefun(symbol, index, shape, prefix):
    symbol = prefix
    result = [sp.Symbol(symbol + ''.join(ixnam(n) for n in index + [i]))
            if symbol else sp.sympify(0) for i in range(shape[0])]
    return result

d = ["x", "y", "z"]
for i in range(3):
    for j in range(i,3):
        n1 = f"metric{i}{j}"
        n2 = f"g{d[i]}{d[j]}"
        grid.rename(n1,n2)

# TODO: turn these into a Cactus parameters
set_Theta_zero = thorn.declare_param('set_Theta_zero', default=False, doc="set_Theta_zero")
kappa1 = thorn.declare_param('kappa1', default=0.02, doc="kappa1")
kappa2 = thorn.declare_param('kappa2', default=0.0, doc="kappa2")
f_mu_L = thorn.declare_param('f_mu_L', default=2.0, doc="f_mu_L")
f_mu_S = thorn.declare_param('f_mu_S', default=1.0, doc="f_mu_S")
eta = thorn.declare_param('eta', default=2.0, doc="eta")
alphaG_floor = thorn.declare_param('alphaG_floor', default=1.0e-10, vmin=0.0, doc="alphaG_floor")
chi_floor = thorn.declare_param('chi_floor', default=1.0e-10, vmin=0.0, doc="chi_floor")

gDD = ixp.register_gridfunctions_for_single_rankN(2, "EXTERNAL", "metric", "sym01", centering=centering, external_module="ADMBase", namefun=lambda sym,ind,shape: namefun(sym,ind,shape, "g"))
KDD = ixp.register_gridfunctions_for_single_rankN(2, "EXTERNAL", "extcurv", "sym01", centering=centering, external_module="ADMBase", namefun=lambda sym,ind,shape: namefun(sym,ind,shape, "k"))
alp = thorn.register_gridfunctions("EXTERNAL", ["alp"], centering=centering, external_module="ADMBase")
dtalp = thorn.register_gridfunctions("EXTERNAL", ["dtalp"], centering=centering, external_module="ADMBase")
betaU = ixp.register_gridfunctions_for_single_rankN(1, "EXTERNAL", "shift", centering=centering, external_module="ADMBase", namefun=lambda sym,ind,shape: namefun(sym,ind,shape, "beta"))
dtbetaU = ixp.register_gridfunctions_for_single_rankN(1, "EXTERNAL", "dtshift", centering=centering, external_module="ADMBase", namefun=lambda sym,ind,shape: namefun(sym,ind,shape, "dtbeta"))

eTtt = thorn.register_gridfunctions("EXTERNAL", ["eTtt"], centering=centering, external_module="TmunuBase")
eTtiD = ixp.register_gridfunctions_for_single_rankN(1, "EXTERNAL", "eTti", centering=centering, external_module="TmunuBase", namefun=lambda sym,ind,shape: namefun(sym,ind,shape, "eTt"))
eTijDD = ixp.register_gridfunctions_for_single_rankN(2, "EXTERNAL", "eTtij", "sym01", centering=centering, external_module="TmunuBase", namefun=lambda sym,ind,shape: namefun(sym,ind,shape, "eT"))

chi = thorn.register_gridfunctions("EVOL", ["chi"], centering=centering)
gammatildeDD = ixp.register_gridfunctions_for_single_rankN(2, "EVOL", "gammatildeDD", "sym01", centering=centering)
Khat = thorn.register_gridfunctions("EVOL", ["Khat"], centering=centering)
AtildeDD = ixp.register_gridfunctions_for_single_rankN(2, "EVOL", "AtildeDD", "sym01", centering=centering)
GammatildeU = ixp.register_gridfunctions_for_single_rankN(1, "EVOL", "GammatildeU", centering=centering)
Theta = thorn.register_gridfunctions("EVOL", ["Theta"], centering=centering)
alphaG = thorn.register_gridfunctions("EVOL", ["alphaG"], centering=centering)
betaGU = ixp.register_gridfunctions_for_single_rankN(1, "EVOL", "betaGU", centering=centering)

chi_rhs = thorn.register_gridfunctions("AUX", ["chi_rhs"], centering=centering)
gammatildeDD_rhs = ixp.register_gridfunctions_for_single_rankN(2, "AUX", "gammatildeDD_rhs", "sym01", centering=centering)
Khat_rhs = thorn.register_gridfunctions("AUX", ["Khat_rhs"], centering=centering)
AtildeDD_rhs = ixp.register_gridfunctions_for_single_rankN(2, "AUX", "AtildeDD_rhs", "sym01", centering=centering)
GammatildeU_rhs = ixp.register_gridfunctions_for_single_rankN(1, "AUX", "GammatildeU_rhs", centering=centering)
Theta_rhs = thorn.register_gridfunctions("AUX", ["Theta_rhs"], centering=centering)
alphaG_rhs = thorn.register_gridfunctions("AUX", ["alphaG_rhs"], centering=centering)
betaGU_rhs = ixp.register_gridfunctions_for_single_rankN(1, "AUX", "betaGU_rhs", centering=centering)

ZetatildeCU = ixp.register_gridfunctions_for_single_rankN(1, "AUX", "ZetatildeCU", centering=centering)
HC = thorn.register_gridfunctions("AUX", ["HC"], centering=centering)
MtildeCU = ixp.register_gridfunctions_for_single_rankN(1, "AUX", "MtildeCU", centering=centering)
allC = thorn.register_gridfunctions("AUX", ["allC"], centering=centering)

# Derivatives as calculated by NRPy
chi_dD = ixp.declarerank1("chi_dD")
chi_dDD = ixp.declarerank2("chi_dDD", "sym01")
gammatildeDD_dD = ixp.declarerank3("gammatildeDD_dD", "sym01")
gammatildeDD_dDD = ixp.declarerank4("gammatildeDD_dDD", "sym01_sym23")
Khat_dD = ixp.declarerank1("Khat_dD")
AtildeDD_dD = ixp.declarerank3("AtildeDD_dD", "sym01")
Theta_dD = ixp.declarerank1("Theta_dD")
GammatildeU_dD = ixp.declarerank2("GammatildeU_dD", None)
alphaG_dD = ixp.declarerank1("alphaG_dD")
alphaG_dDD = ixp.declarerank2("alphaG_dDD", "sym01")
betaGU_dD = ixp.declarerank2("betaGU_dD", None)
betaGU_dDD = ixp.declarerank3("betaGU_dDD", "sym12")

# Our tile-local derivatives
dchiD = ixp.register_gridfunctions_for_single_rankN(1, "TILE_TMP", "dchiD", centering=centering)
ddchiDD = ixp.register_gridfunctions_for_single_rankN(2, "TILE_TMP", "ddchiDD", "sym01", centering=centering)
dgammatildeDDD = ixp.register_gridfunctions_for_single_rankN(3, "TILE_TMP", "dgammatildeDDD", "sym01", centering=centering)
ddgammatildeDDDD = ixp.register_gridfunctions_for_single_rankN(4, "TILE_TMP", "ddgammatildeDDDD", "sym01_sym23", centering=centering)
dKhatD = ixp.register_gridfunctions_for_single_rankN(1, "TILE_TMP", "dKhatD", centering=centering)
dAtildeDDD = ixp.register_gridfunctions_for_single_rankN(3, "TILE_TMP", "dAtildeDDD", "sym01", centering=centering)
dThetaD = ixp.register_gridfunctions_for_single_rankN(1, "TILE_TMP", "dThetaD", centering=centering)
dGammatildeUD = ixp.register_gridfunctions_for_single_rankN(2, "TILE_TMP", "dGammatildeUD", None, centering=centering)
dalphaGD = ixp.register_gridfunctions_for_single_rankN(1, "TILE_TMP", "dalphaGD", centering=centering)
ddalphaGDD = ixp.register_gridfunctions_for_single_rankN(2, "TILE_TMP", "ddalphaGDD", "sym01", centering=centering)
dbetaGUD = ixp.register_gridfunctions_for_single_rankN(2, "TILE_TMP", "dbetaGUD", None, centering=centering)
ddbetaGUDD = ixp.register_gridfunctions_for_single_rankN(3, "TILE_TMP", "ddbetaGUDD", "sym12", centering=centering)

# Expressions

gUU_expr, detg_expr = ixp.symm_matrix_inverter3x3(gDD)
gammatildeUU_expr, detgammatilde_expr = ixp.symm_matrix_inverter3x3(gammatildeDD)

# Initial conditions

gUU = ixp.register_gridfunctions_for_single_rankN(2, "SCALAR_TMP", "gUU", "sym01", centering=centering)
gammatildeUU = ixp.register_gridfunctions_for_single_rankN(2, "SCALAR_TMP", "gammatildeUU", "sym01", centering=centering)
dgammatildeUUD = ixp.register_gridfunctions_for_single_rankN(3, "SCALAR_TMP", "dgammatildeUUD", "sym01", centering=centering)
Theta_val = thorn.register_gridfunctions("SCALAR_TMP", ["Theta_val"], centering=centering)
trK = thorn.register_gridfunctions("SCALAR_TMP", ["trK"], centering=centering)
chi_val = thorn.register_gridfunctions("SCALAR_TMP", ["chi_val"], centering=centering)

def Initial():
    initial1_eqns = flatten([
        [lhrh(lhs=chi_val, rhs=1 / cbrt(detg_expr))],
        [lhrh(lhs=gUU[i][j], rhs=gUU_expr[i][j]) for i in range(3) for j in range(i+1)],
        [lhrh(lhs=Theta_val, rhs=sympify(0))],
        [lhrh(lhs=trK, rhs=sum2_symm(lambda x, y: gUU[x][y] * KDD[x][y]))],
        [lhrh(lhs=chi, rhs=chi_val)],
        [lhrh(lhs=gammatildeDD[i][j], rhs=chi_val * gDD[i][j]) for i in range(3) for j in range(i+1)],
        [lhrh(lhs=Theta, rhs=Theta_val)],
        [lhrh(lhs=Khat, rhs=trK - 2 * Theta_val)],
        [lhrh(lhs=AtildeDD[i][j], rhs=chi_val * (KDD[i][j] - trK / 3 * gDD[i][j])) for i in range(3) for j in range(i+1)],
        [lhrh(lhs=alphaG, rhs=alp)],
        [lhrh(lhs=betaGU[i], rhs=betaU[i]) for i in range(3)],
    ])

    thorn.add_func(
        "Z4cNRPy_Initial1",
        body=initial1_eqns,
        where='everywhere',
        schedule_bin="Z4cNRPy_InitialGroup",
        doc="Convert ADM to Z4c variables, part 1",
        centering=centering)
    
    initial2_eqns = flatten([
        [[
            [lhrh(lhs=dgammatildeDDD[i][j][k], rhs=gammatildeDD_dD[i][j][k]) for k in range(3)],
            [loop],
        ] for i in range(3) for j in range(i+1)],
        [lhrh(lhs=gammatildeUU[i][j],
              rhs=gammatildeUU_expr[i][j])
         for i in range(3) for j in range(i+1)],
        [lhrh(lhs=dgammatildeUUD[i][j][k],
              rhs=- sum2_symm(lambda x, y: gammatildeUU[i][x] * gammatildeUU[j][y] * dgammatildeDDD[x][y][k]))
         for i in range(3) for j in range(i+1) for k in range(3)],
        [[lhrh(lhs=GammatildeU[i],
               rhs=- sum1(lambda x: dgammatildeUUD[i][x][x]))]
         for i in range(3)],
    ])
    
    thorn.add_func(
        "Z4cNRPy_Initial2",
        body=initial2_eqns,
        where='interior',
        schedule_bin="Z4cNRPy_InitialGroup AFTER Z4cNRPy_Initial1",
        doc="Convert ADM to Z4c variables, part 2",
        sync="GammatildeUGF",
        centering=centering)
Initial()

# Enforce constraints

detgammatilde = thorn.register_gridfunctions("SCALAR_TMP", ["detgammatilde"], centering=centering)
trAtilde = thorn.register_gridfunctions("SCALAR_TMP", ["trAtilde"], centering=centering)

def Enforce():
    enforce_eqns = flatten([
        # Enforce floors
        [lhrh(lhs=gammatildeUU[i][j], rhs=gammatildeUU_expr[i][j]) for i in range(3) for j in range(i+1)],
        [lhrh(lhs=detgammatilde, rhs=detgammatilde_expr)],
        [lhrh(lhs=trAtilde, rhs=sum2_symm(lambda x, y: gammatildeUU[x][y] * AtildeDD[x][y]))],
        [lhrh(lhs=chi, rhs=Max(chi_floor, chi))],
        [lhrh(lhs=alphaG, rhs=Max(alphaG_floor, alphaG))],
        # Enforce algebraic constraints; see arXiv:1212.2901 [gr-qc]
        [lhrh(lhs=gammatildeDD[i][j], rhs=(1 / cbrt(detgammatilde)) * gammatildeDD[i][j]) for i in range(3) for j in range(i+1)],
        [lhrh(lhs=AtildeDD[i][j], rhs=(AtildeDD[i][j] - trAtilde / 3 * gammatildeDD[i][j])) for i in range(3) for j in range(i+1)],
    ])
    
    thorn.add_func(
        "Z4cNRPy_Enforce",
        body=enforce_eqns,
        where='interior',
        schedule_bin="Z4cNRPy_PostStepGroup",
        doc="Enforce constraints",
        sync="chiGF gammatildeDDGF KhatGF AtildeDDGF GammatildeUGF ThetaGF alphaGGF betaGUGF",
        centering=centering
        )
Enforce()

# Calculate ADM variables
def ADM():
    adm_eqns = flatten([
        [lhrh(lhs=gDD[i][j], rhs=1 / chi * gammatildeDD[i][j]) for i in range(3) for j in range(i+1)],
        [lhrh(lhs=KDD[i][j], rhs=1 / chi * (AtildeDD[i][j] + (Khat + 2 * Theta) / 3 * gammatildeDD[i][j]))
         for i in range(3) for j in range(i+1)],
        [lhrh(lhs=alp, rhs=alphaG)],
        [lhrh(lhs=dtalp, rhs=-alphaG * f_mu_L * Khat)],
        [lhrh(lhs=betaU[i], rhs=betaGU[i]) for i in range(3)],
        [lhrh(lhs=dtbetaU[i], rhs=f_mu_S * GammatildeU[i] - eta * betaGU[i]) for i in range(3)],
    ])
    
    thorn.add_func(
        "Z4cNRPy_ADM",
        body=adm_eqns,
        where='interior',
        schedule_bin="Z4cNRPy_PostStepGroup AFTER Z4cNRPy_Enforce",
        doc="Calculate ADM variables",
        sync="ADMBase::metric ADMBase::curv ADMBase::lapse ADMBase::dtlapse ADMBase::shift ADMBase::dtshift",
        centering=centering)
ADM()

# Calculate RHS

dgDDD = ixp.register_gridfunctions_for_single_rankN(3, "SCALAR_TMP", "dgDDD", "sym01", centering=centering)
GammaDDD = ixp.register_gridfunctions_for_single_rankN(3, "SCALAR_TMP", "GammaDDD", "sym12", centering=centering)
GammaUDD = ixp.register_gridfunctions_for_single_rankN(3, "SCALAR_TMP", "GammaUDD", "sym12", centering=centering)
DDalphaGDD = ixp.register_gridfunctions_for_single_rankN(2, "SCALAR_TMP", "DDalphaGDD", "sym01", centering=centering)
DbetaGUD = ixp.register_gridfunctions_for_single_rankN(2, "SCALAR_TMP", "DbetaGUD", None, centering=centering)
rho = thorn.register_gridfunctions("SCALAR_TMP", ["rho"], centering=centering)
SiD = ixp.register_gridfunctions_for_single_rankN(1, "SCALAR_TMP", "SiD", None, centering=centering)
SijDD = ixp.register_gridfunctions_for_single_rankN(2, "SCALAR_TMP", "SijDD", "sym01", centering=centering)
trS = thorn.register_gridfunctions("SCALAR_TMP", ["trS"], centering=centering)
AtildeUD = ixp.register_gridfunctions_for_single_rankN(2, "SCALAR_TMP", "AtildeUD", None, centering=centering)
AtildeUU = ixp.register_gridfunctions_for_single_rankN(2, "SCALAR_TMP", "AtildeUU", "sym01", centering=centering)
dAtildeUDD = ixp.register_gridfunctions_for_single_rankN(3, "SCALAR_TMP", "dAtildeUDD", None, centering=centering)
dAtildeUUD = ixp.register_gridfunctions_for_single_rankN(3, "SCALAR_TMP", "dAtildeUUD", "sym01", centering=centering)
GammatildeDDD = ixp.register_gridfunctions_for_single_rankN(3, "SCALAR_TMP", "GammatildeDDD", "sym12", centering=centering)
GammatildeUDD = ixp.register_gridfunctions_for_single_rankN(3, "SCALAR_TMP", "GammatildeUDD", "sym12", centering=centering)
GammatildeUDU = ixp.register_gridfunctions_for_single_rankN(3, "SCALAR_TMP", "GammatildeUDU", None, centering=centering)
GammatildedirectU = ixp.register_gridfunctions_for_single_rankN(1, "SCALAR_TMP", "GammatildedirectU", None, centering=centering)
DDchiDD = ixp.register_gridfunctions_for_single_rankN(2, "SCALAR_TMP", "DDchiDD", "sym01", centering=centering)
RchiDD = ixp.register_gridfunctions_for_single_rankN(2, "SCALAR_TMP", "RchiDD", "sym01", centering=centering)
RtildeDD = ixp.register_gridfunctions_for_single_rankN(2, "SCALAR_TMP", "RtildeDD", "sym01", centering=centering)
RDD = ixp.register_gridfunctions_for_single_rankN(2, "SCALAR_TMP", "RDD", "sym01", centering=centering)
Rsc = thorn.register_gridfunctions("SCALAR_TMP", ["Rsc"], centering=centering)
alphaRicciTmunuDD = ixp.register_gridfunctions_for_single_rankN(2, "SCALAR_TMP", "alphaRicciTmunuDD", "sym01", centering=centering)
ZetatildeCvalU = ixp.register_gridfunctions_for_single_rankN(1, "SCALAR_TMP", "ZetatildeCvalU", None, centering=centering)
HCval = thorn.register_gridfunctions("SCALAR_TMP", ["HCval"], centering=centering)
MtildeCvalU = ixp.register_gridfunctions_for_single_rankN(1, "SCALAR_TMP", "MtildeCvalU", None, centering=centering)

def RHS():
    rhs_eqns = flatten([
        # Derivatives
        [
            [lhrh(lhs=dchiD[i], rhs=chi_dD[i]) for i in range(3)],
            [lhrh(lhs=ddchiDD[i][j], rhs=chi_dDD[i][j]) for i in range(3) for j in range(i+1)],
            [loop],
        ],
        [[
            [lhrh(lhs=dgammatildeDDD[i][j][k], rhs=gammatildeDD_dD[i][j][k]) for k in range(3)],
            [lhrh(lhs=ddgammatildeDDDD[i][j][k][l], rhs=gammatildeDD_dDD[i][j][k][l]) for k in range(3) for l in range(k+1)],
            [loop],
        ] for i in range(3) for j in range(i+1)],
        [[
            [lhrh(lhs=dGammatildeUD[i][j], rhs=GammatildeU_dD[i][j]) for j in range(3)],
            [loop],
        ] for i in range(3)],
        [
            [lhrh(lhs=dKhatD[i], rhs=Khat_dD[i]) for i in range(3)],
            [loop],
        ],
        [[
            [lhrh(lhs=dAtildeDDD[i][j][k], rhs=AtildeDD_dD[i][j][k]) for k in range(3)],
            [loop],
        ] for i in range(3) for j in range(i+1)],
        [
            [lhrh(lhs=dThetaD[i], rhs=Theta_dD[i]) for i in range(3)],
            [loop],
        ],
        [
            [lhrh(lhs=dalphaGD[i], rhs=alphaG_dD[i]) for i in range(3)],
            [lhrh(lhs=ddalphaGDD[i][j], rhs=alphaG_dDD[i][j]) for i in range(3) for j in range(i+1)],
            [loop],
        ],
        [[
            [lhrh(lhs=dbetaGUD[i][j], rhs=betaGU_dD[i][j]) for j in range(3)],
            [lhrh(lhs=ddbetaGUDD[i][j][k], rhs=betaGU_dDD[i][j][k]) for j in range(3) for k in range(j+1)],
            [loop],
        ] for i in range(3)],

        # RHS
        [lhrh(lhs=gammatildeUU[i][j],
              rhs=gammatildeUU_expr[i][j])
         for i in range(3) for j in range(i+1)],
        [lhrh(lhs=gUU[i][j],
              rhs=chi * gammatildeUU[i][j])
         for i in range(3) for j in range(i+1)],
        [lhrh(lhs=dgDDD[i][j][k],
              rhs=(- 1 / chi**2 * dchiD[k] * gammatildeDD[i][j]
                   + 1 / chi * dgammatildeDDD[i][j][k]))
         for i in range(3) for j in range(i+1) for k in range(3)],
        [lhrh(lhs=GammaDDD[i][j][k],
              rhs=Rational(1,2) * (dgDDD[i][j][k] + dgDDD[i][k][j] - dgDDD[j][k][i]))
         for i in range(3) for j in range(3) for k in range(j+1)],
        [lhrh(lhs=GammaUDD[i][j][k],
              rhs=sum1(lambda x: gUU[i][x] * GammaDDD[x][j][k]))
         for i in range(3) for j in range(3) for k in range(j+1)],
        [lhrh(lhs=DbetaGUD[i][j],
              rhs=dbetaGUD[i][j] + sum1(lambda x: GammaUDD[i][x][j] * betaGU[x]))
         for i in range(3) for j in range(3)],
        # arXiv:1212.2901 [gr-qc], (1)
        [lhrh(lhs=chi_rhs,
              rhs=Rational(2,3) * chi * (alphaG * (Khat + 2 * Theta) - sum1(lambda x: DbetaGUD[x][x])))],

        # arXiv:1212.2901 [gr-qc], (2)
        [lhrh(lhs=gammatildeDD_rhs[i][j],
              rhs=(-2 * alphaG * AtildeDD[i][j]
                   + sum1(lambda x: betaGU[x] * dgammatildeDDD[i][j][x])
                   + sum1(lambda x: gammatildeDD[x][i] * dbetaGUD[x][j] + gammatildeDD[x][j] * dbetaGUD[x][i])
                   - Rational(2,3) * sum1(lambda x: gammatildeDD[i][j] * dbetaGUD[x][x])))
         for i in range(3) for j in range(i+1)],

        [lhrh(lhs=DDalphaGDD[i][j],
              rhs=ddalphaGDD[i][j] - sum1(lambda x: GammaUDD[x][i][j] * dalphaGD[x]))
         for i in range(3) for j in range(i+1)],
        [lhrh(lhs=AtildeUD[i][j],
              rhs=sum1(lambda x: gammatildeUU[i][x] * AtildeDD[x][j]))
         for i in range(3) for j in range(3)],
        [lhrh(lhs=rho,
              rhs=1 / alphaG**2 * (+ eTtt
                                   - 2 * sum1(lambda x: betaGU[x] * eTtiD[x])
                                   + sum2_symm(lambda x, y: betaGU[x] * betaGU[y] * eTijDD[x][y])))],
        [lhrh(lhs=SiD[i],
              rhs=-1 / alphaG * (eTtiD[i] - sum1(lambda x:  betaGU[x] * eTijDD[i][x])))
         for i in range(3)],
        [lhrh(lhs=SijDD[i][j],
              rhs=eTijDD[i][j])
         for i in range(3) for j in range(i+1)],
        [lhrh(lhs=trS,
              rhs=sum2_symm(lambda x, y: gUU[x][y] * SijDD[x][y]))],
        # arXiv:1212.2901 [gr-qc], (3)
        [lhrh(lhs=Khat_rhs,
              rhs=(- sum2_symm(lambda x, y: gUU[x][y] * DDalphaGDD[x][y])
                   + alphaG * (+ sum2_symm(lambda x, y: AtildeUD[x][y] * AtildeUD[y][x])
                               + Rational(1,3) * (Khat + 2 * Theta)**2)
                   + 4 * pi * alphaG * (trS + rho)
                   + alphaG * kappa1 * (1 - kappa2) * Theta)
                   + sum1(lambda x: betaGU[x] * dKhatD[x]))],

        [lhrh(lhs=GammatildeDDD[i][j][k],
              rhs=Rational(1,2) * (dgammatildeDDD[i][j][k] + dgammatildeDDD[i][k][j] - dgammatildeDDD[j][k][i]))
         for i in range(3) for j in range(3) for k in range(j+1)],
        [lhrh(lhs=GammatildeUDD[i][j][k],
              rhs=sum1(lambda x: gammatildeUU[i][x] * GammatildeDDD[x][j][k]))
         for i in range(3) for j in range(3) for k in range(j+1)],
        [lhrh(lhs=GammatildeUDU[i][j][k],
              rhs=sum1(lambda x: GammatildeUDD[i][j][x] * gammatildeUU[x][k]))
         for i in range(3) for j in range(3) for k in range(3)],
        [lhrh(lhs=GammatildedirectU[i],
              rhs=sum1(lambda x: GammatildeUDU[i][x][x]))
         for i in range(3)],
        [lhrh(lhs=DDchiDD[i][j],
              rhs=ddchiDD[i][j] - sum1(lambda x: GammatildeUDD[x][i][j] * dchiD[x]))
         for i in range(3) for j in range(i+1)],
        # arXiv:1212.2901 [gr-qc], (8)
        [lhrh(lhs=RchiDD[i][j],
              rhs=(+ Rational(1,2) / chi * DDchiDD[i][j]
                   + Rational(1,2) / chi * gammatildeDD[i][j] * sum2_symm(lambda x, y: gammatildeUU[x][y] * DDchiDD[x][y])
                   - Rational(1,4) / chi**2 * dchiD[i] * dchiD[j]
                   - Rational(3,4) / chi**2 * (gammatildeDD[i][j] *
                                               sum2_symm(lambda x, y: gammatildeUU[x][y] * dchiD[x] * dchiD[y]))))
         for i in range(3) for j in range(i+1)],
        # arXiv:1212.2901 [gr-qc], (9)
        [lhrh(lhs=RtildeDD[i][j],
              rhs=(- Rational(1,2) * sum2_symm(lambda x, y: gammatildeUU[x][y] * ddgammatildeDDDD[i][j][x][y])
                   + Rational(1,2) * sum1(lambda x: (+ gammatildeDD[x][i] * dGammatildeUD[x][j]
                                                     + gammatildeDD[x][j] * dGammatildeUD[x][i]))
                   + Rational(1,2) * sum1(lambda x: (+ GammatildedirectU[x] * GammatildeDDD[i][j][x]
                                                     + GammatildedirectU[x] * GammatildeDDD[j][i][x]))
                   + sum2_symm(lambda x, y: (+ GammatildeUDD[x][i][y] * GammatildeUDU[j][x][y]
                                             + GammatildeUDD[x][j][y] * GammatildeUDU[i][x][y]
                                             + GammatildeUDD[x][i][y] * GammatildeUDU[x][j][y]))))
         for i in range(3) for j in range(i+1)],
        # arXiv:1212.2901 [gr-qc], (7)
        [lhrh(lhs=RDD[i][j],
              rhs=RchiDD[i][j] + RtildeDD[i][j])
         for i in range(3) for j in range(i+1)],
        [lhrh(lhs=Rsc,
              rhs=sum2_symm(lambda x, y: gUU[x][y] * RDD[x][y]))],
        [lhrh(lhs=alphaRicciTmunuDD[i][j],
              rhs=- DDalphaGDD[i][j] + alphaG * (RDD[i][j] - 8 * pi * SijDD[i][j]))
         for i in range(3) for j in range(i+1)],
        # arXiv:1212.2901 [gr-qc], (4)
        [lhrh(lhs=AtildeDD_rhs[i][j],
              rhs=(+ chi * (alphaRicciTmunuDD[i][j]
                            - Rational(1,3) * gDD[i][j] * sum2_symm(lambda x, y: gUU[x][y] * alphaRicciTmunuDD[x][y]))
                   + alphaG * (+ (Khat + 2 * Theta) * AtildeDD[i][j]
                               - 2 * sum1(lambda x: AtildeDD[x][i] * AtildeUD[x][j]))
                   + sum1(lambda x: betaGU[x] * dAtildeDDD[i][j][x])
                   + sum1(lambda x: AtildeDD[x][i] * dbetaGUD[x][j] + AtildeDD[x][j] * dbetaGUD[x][i])
                   - Rational(2,3) * AtildeDD[i][j] * sum1(lambda x: dbetaGUD[x][x])))
         for i in range(3) for j in range(i+1)],

        [lhrh(lhs=AtildeUU[i][j],
              rhs=sum1(lambda x: gammatildeUU[i][x] * AtildeUD[j][x]))
         for i in range(3) for j in range(i+1)],
        # arXiv:1212.2901 [gr-qc], (5)
        [lhrh(lhs=GammatildeU_rhs[i],
              rhs=(- 2 * sum1(lambda x: AtildeUU[i][x] * dalphaGD[x])
                   + 2 * alphaG *
                         (+ sum2_symm(lambda x, y: GammatildeUDD[i][x][y] * AtildeUU[x][y])
                          - Rational(3,2) / chi * sum1(lambda x: AtildeUU[i][x] * dchiD[x])
                          - Rational(1,3) * sum1(lambda x: gammatildeUU[i][x] * (2 * dKhatD[x] + dThetaD[x]))
                          - 8 * pi * sum1(lambda x: gammatildeUU[i][x] * SiD[x]))
                   + sum2_symm(lambda x, y: gammatildeUU[x][y] * ddbetaGUDD[i][x][y])
                   + Rational(1,3) * sum2(lambda x, y: gammatildeUU[i][x] * ddbetaGUDD[y][x][y])
                   + sum1(lambda x: betaGU[x] * dGammatildeUD[i][x])
                   - sum1(lambda x: GammatildedirectU[x] * dbetaGUD[i][x])
                   + Rational(2,3) * GammatildedirectU[i] * sum1(lambda x: dbetaGUD[x][x])
                   - 2 * alphaG * kappa1 * (GammatildeU[i] - GammatildedirectU[i])))
         for i in range(3)],

        # arXiv:1212.2901 [gr-qc], (6)
        [lhrh(lhs=Theta_rhs,
              rhs=Piecewise((0, set_Theta_zero),
                            ((+ Rational(1,2) * alphaG * (+ Rsc
                                                          - sum2_symm(lambda x, y: AtildeUU[x][y] * AtildeDD[x][y])
                                                          + Rational(2,3) * (Khat + 2 * Theta)**2)
                              - alphaG * (+ 8 * pi * rho
                                          + kappa1 * (2 + kappa2) * Theta)
                              + sum1(lambda x: betaGU[x] * dThetaD[x])), True)))],

        [lhrh(lhs=alphaG_rhs,
              rhs=-alphaG * f_mu_L * Khat)],

        [lhrh(lhs=betaGU_rhs[i],
              rhs=f_mu_S * GammatildeU[i] - eta * betaGU[i])
         for i in range(3)],
    ])
    
    thorn.add_func(
        "Z4cNRPy_RHS",
        body=rhs_eqns,
        where='interior',
        schedule_bin="ODESolvers_RHS",
        doc="Calculate RHS",
        sync="chi_rhsGF gammatildeDD_rhsGF Khat_rhsGF AtildeDD_rhsGF GammatildeU_rhsGF Theta_rhsGF alphaG_rhsGF betaGU_rhsGF",
        centering=centering)
RHS()

# Evaluate constraints

def Constraints():
    constraints_eqns = flatten([
        # Derivatives
        [
            [lhrh(lhs=dchiD[i], rhs=chi_dD[i]) for i in range(3)],
            [lhrh(lhs=ddchiDD[i][j], rhs=chi_dDD[i][j]) for i in range(3) for j in range(i+1)],
            [loop],
        ],
        [[
            [lhrh(lhs=dgammatildeDDD[i][j][k], rhs=gammatildeDD_dD[i][j][k]) for k in range(3)],
            [lhrh(lhs=ddgammatildeDDDD[i][j][k][l], rhs=gammatildeDD_dDD[i][j][k][l]) for k in range(3) for l in range(k+1)],
            [loop],
        ] for i in range(3) for j in range(i+1)],
        [[
            [lhrh(lhs=dGammatildeUD[i][j], rhs=GammatildeU_dD[i][j]) for j in range(3)],
            [loop],
        ] for i in range(3)],
        [[
            [lhrh(lhs=dAtildeDDD[i][j][k], rhs=AtildeDD_dD[i][j][k]) for k in range(3)],
            [loop],
        ] for i in range(3) for j in range(i+1)],
        [
            [lhrh(lhs=dKhatD[i], rhs=Khat_dD[i]) for i in range(3)],
            [loop],
        ],
        [
            [lhrh(lhs=dThetaD[i], rhs=Theta_dD[i]) for i in range(3)],
            [loop],
        ],

        # Constraints

        # ZetatildeC, eqn. (13)
        [lhrh(lhs=gammatildeUU[i][j],
              rhs=gammatildeUU_expr[i][j])
         for i in range(3) for j in range(i+1)],
        [lhrh(lhs=GammatildeDDD[i][j][k],
              rhs=Rational(1,2) * (dgammatildeDDD[i][j][k] + dgammatildeDDD[i][k][j] - dgammatildeDDD[j][k][i]))
         for i in range(3) for j in range(3) for k in range(j+1)],
        [lhrh(lhs=GammatildeUDD[i][j][k],
              rhs=sum1(lambda x: gammatildeUU[i][x] * GammatildeDDD[x][j][k]))
         for i in range(3) for j in range(3) for k in range(j+1)],
        [lhrh(lhs=GammatildeUDU[i][j][k],
              rhs=sum1(lambda x: GammatildeUDD[i][j][x] * gammatildeUU[x][k]))
         for i in range(3) for j in range(3) for k in range(3)],
        [lhrh(lhs=GammatildedirectU[i],
              rhs=sum1(lambda x: GammatildeUDU[i][x][x]))
         for i in range(3)],
        [lhrh(lhs=ZetatildeCvalU[i],
              rhs=(GammatildeU[i] - GammatildedirectU[i]) / 2)
         for i in range(3)],
        [lhrh(lhs=ZetatildeCU[i],
              rhs=(ZetatildeCvalU[i]))
         for i in range(3)],

        # HC, eqn. (14)
        [lhrh(lhs=gUU[i][j],
              rhs=chi * gammatildeUU[i][j])
         for i in range(3) for j in range(i+1)],
        [lhrh(lhs=dgDDD[i][j][k],
              rhs=(- 1 / chi**2 * dchiD[k] * gammatildeDD[i][j]
                   + 1 / chi * dgammatildeDDD[i][j][k]))
         for i in range(3) for j in range(i+1) for k in range(3)],
        [lhrh(lhs=GammaDDD[i][j][k],
              rhs=Rational(1,2) * (dgDDD[i][j][k] + dgDDD[i][k][j] - dgDDD[j][k][i]))
         for i in range(3) for j in range(3) for k in range(j+1)],
        [lhrh(lhs=GammaUDD[i][j][k],
              rhs=sum1(lambda x: gUU[i][x] * GammaDDD[x][j][k]))
         for i in range(3) for j in range(3) for k in range(j+1)],
        [lhrh(lhs=AtildeUD[i][j],
              rhs=sum1(lambda x: gammatildeUU[i][x] * AtildeDD[x][j]))
         for i in range(3) for j in range(3)],
        [lhrh(lhs=rho,
              rhs=1 / alphaG**2 * (+ eTtt
                                   - 2 * sum1(lambda x: betaGU[x] * eTtiD[x])
                                   + sum2_symm(lambda x, y: betaGU[x] * betaGU[y] * eTijDD[x][y])))],
        [lhrh(lhs=DDchiDD[i][j],
              rhs=ddchiDD[i][j] - sum1(lambda x: GammatildeUDD[x][i][j] * dchiD[x]))
         for i in range(3) for j in range(i+1)],
        # arXiv:1212.2901 [gr-qc], (8)
        [lhrh(lhs=RchiDD[i][j],
              rhs=(+ Rational(1,2) / chi * DDchiDD[i][j]
                   + Rational(1,2) / chi * gammatildeDD[i][j] * sum2_symm(lambda x, y: gammatildeUU[x][y] * DDchiDD[x][y])
                   - Rational(1,4) / chi**2 * dchiD[i] * dchiD[j]
                   - Rational(3,4) / chi**2 * (gammatildeDD[i][j] *
                                               sum2_symm(lambda x, y:  gammatildeUU[x][y] * dchiD[x] * dchiD[y]))))
         for i in range(3) for j in range(i+1)],
        # arXiv:1212.2901 [gr-qc], (9)
        [lhrh(lhs=RtildeDD[i][j],
              rhs=(- Rational(1,2) * sum2_symm(lambda x, y: gammatildeUU[x][y] * ddgammatildeDDDD[i][j][x][y])
                   + Rational(1,2) * sum1(lambda x: (+ gammatildeDD[x][i] * dGammatildeUD[x][j]
                                                     + gammatildeDD[x][j] * dGammatildeUD[x][i]))
                   + Rational(1,2) * sum1(lambda x: (+ GammatildedirectU[x] * GammatildeDDD[i][j][x]
                                                     + GammatildedirectU[x] * GammatildeDDD[j][i][x]))
                   + sum2_symm(lambda x, y: (+ GammatildeUDD[x][i][y] * GammatildeUDU[j][x][y]
                                             + GammatildeUDD[x][j][y] * GammatildeUDU[i][x][y]
                                             + GammatildeUDD[x][i][y] * GammatildeUDU[x][j][y]))))
         for i in range(3) for j in range(i+1)],
        # arXiv:1212.2901 [gr-qc], (7)
        [lhrh(lhs=RDD[i][j],
              rhs=RchiDD[i][j] + RtildeDD[i][j])
         for i in range(3) for j in range(i+1)],
        [lhrh(lhs=Rsc,
              rhs=sum2_symm(lambda x, y: gUU[x][y] * RDD[x][y]))],
        # arXiv:1212.2901 [gr-qc], (14)
        [lhrh(lhs=HCval,
              rhs=(+ Rsc
                   + sum2_symm(lambda x, y: AtildeUD[x][y] * AtildeUD[y][x])
                   - Rational(2,3) * (Khat + 2 * Theta)**2
                   - 16 * pi * rho))],
        [lhrh(lhs=HC,
              rhs=HCval)],

        # MtildeC, eqn. (15)
        [lhrh(lhs=AtildeUU[i][j],
              rhs=sum1(lambda x: gammatildeUU[i][x] * AtildeUD[j][x]))
         for i in range(3) for j in range(i+1)],
        [lhrh(lhs=dgammatildeUUD[i][j][k],
              rhs=- sum2_symm(lambda x, y: gammatildeUU[i][x] * gammatildeUU[j][y] * dgammatildeDDD[x][y][k]))
         for i in range(3) for j in range(i+1) for k in range(3)],
        [lhrh(lhs=dAtildeUDD[i][j][k],
              rhs=(+ sum1(lambda x: dgammatildeUUD[i][x][k] * AtildeDD[x][j])
                   + sum1(lambda x: gammatildeUU[i][x] * dAtildeDDD[x][j][k])))
         for i in range(3) for j in range(3) for k in range(3)],
        [lhrh(lhs=dAtildeUUD[i][j][k],
              rhs=(+ sum1(lambda x: dgammatildeUUD[i][x][k] * AtildeUD[j][x])
                   + sum1(lambda x: gammatildeUU[i][x] * dAtildeUDD[j][x][k])))
         for i in range(3) for j in range(i+1) for k in range(3)],
        [lhrh(lhs=SiD[i],
              rhs=-1 / alphaG * (eTtiD[i] - sum1(lambda x:  betaGU[x] * eTijDD[i][x])))
         for i in range(3)],
        # arXiv:1212.2901 [gr-qc], (15)
        [lhrh(lhs=MtildeCvalU[i],
              rhs=(+ sum1(lambda x: dAtildeUUD[i][x][x])
                   + sum2_symm(lambda x, y: GammatildeUDD[i][x][y] * AtildeUU[x][y])
                   - Rational(2,3) * sum1(lambda x: gammatildeUU[i][x] * (dKhatD[x] + 2 * dThetaD[x]))
                   - Rational(2,3) * sum1(lambda x: AtildeUU[i][x] * dchiD[x] / chi)
                   - 8 * pi * sum1(lambda x: gammatildeUU[i][x] * SiD[x])))
         for i in range(3)],
        [lhrh(lhs=MtildeCU[i],
              rhs=MtildeCvalU[i])
         for i in range(3)],

        # allC, arXiv:1111.2177, eqn. (73)
        [lhrh(lhs=allC,
              rhs=sqrt(Max(0, (+ HCval**2
                               + sum2_symm(lambda x, y: gammatildeDD[x][y] * MtildeCvalU[x] * MtildeCvalU[y])
                               + Theta**2
                               + 2 * sum2_symm(lambda x, y: gammatildeDD[x][y] * ZetatildeCvalU[x] * ZetatildeCvalU[y])))))],
    ])
    
    thorn.add_func(
        "Z4cNRPy_Constraints",
        body=constraints_eqns,
        where='interior',
        schedule_bin="Z4c_AnalysisGroup",
        doc="Evaluate Constraints",
        sync="ZetatildeCUGF HCGF MtildeCUGF allCGF",
        centering=centering)
Constraints()

# Generate thorn

assert "CACTUS_HOME" in os.environ, "Please set the CACTUS_HOME variable to point to your Cactus installation"
cactus_home = os.environ["CACTUS_HOME"]
cactus_sim = os.environ.get("CACTUS_SIM", "sim")
cactus_thornlist = os.environ.get("CACTUS_THORNLIST", None)

schedule_raw = """
# Define schedule groups

SCHEDULE GROUP Z4cNRPy_InitialGroup AT initial AFTER ADMBase_PostInitial
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

SCHEDULE GROUP Z4c_AnalysisGroup AT analysis
{
} "Analyse Z4c variables"
"""

thorn.generate(cactus_home, cactus_config=cactus_sim, cactus_thornlist=cactus_thornlist,schedule_raw=schedule_raw)
