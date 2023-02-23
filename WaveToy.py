import os
import grid
from cactusthorn import CactusThorn, loop
from sympy import sympify, sin, cos, pi
import NRPy_param_funcs as par

# Current options are Carpet and CarpetX
grid.ET_driver = os.environ.get("CACTUS_DRIVER","Carpet")

thorn = CactusThorn("TestOne","WaveToyNRPy")

#FD_order = thorn.declare_param('FD_order',default=2,vmin=2,vmax=8,doc="The finite difference order")

wave_speed = thorn.declare_param('wave_speed',default=1,vmin=.1,vmax=100,doc="The speed of the wave")
x0 = thorn.declare_param('x0',default=0,vmin=-100,vmax=100,doc="The x pos of the wave")
y0 = thorn.declare_param('y0',default=0,vmin=-100,vmax=100,doc="The y pos of the wave")
z0 = thorn.declare_param('z0',default=0,vmin=-100,vmax=100,doc="The z pos of the wave")
zero = thorn.declare_param('zero',default=0,vmin=0,vmax=0,doc="zero")

centering='VVC'

# AUXEVOL needed for the evo, can be freed after evaluating rhs (1 time level)
# AUX uu_rhs (1 time level)
# EVOL evolved gfs (3 time levels)
uu_rhs, vv_rhs = thorn.register_gridfunctions("AUX", ["rhs_uu", "rhs_vv"], centering=centering)
uu, vv = thorn.register_gridfunctions("EVOL", ["uu", "vv"], centering=centering)
ana = thorn.register_gridfunctions("AUX", ["ana"], centering=centering)
cctk_time = par.Cparameters("CCTK_REAL","Cactus",["cctk_time"],0)
if grid.ET_driver == "CarpetX":
    tmp0, tmp1 = thorn.register_gridfunctions("TILE_TMP",["tmp0v","tmp1v"],centering=centering)
x,y,z = thorn.get_xyz()

from outputC import lhrh
import indexedexp as ixp
import NRPy_param_funcs as par

FD_order = 2
par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER",FD_order)

uu_dDD = ixp.declarerank2("uu_dDD","sym01")

if grid.ET_driver == "CarpetX":
    evol_eqns = [
        lhrh(lhs=tmp0, rhs=uu_dDD[0][0]),
        lhrh(lhs=tmp1, rhs=uu_dDD[1][1]),
        loop,
        lhrh(lhs=uu_rhs, rhs=vv),
        lhrh(lhs=vv_rhs, rhs=wave_speed**2*(tmp0 + tmp1))
    ]
else:
    # Version of evolution equations without temporaries
    evol_eqns = [
        lhrh(lhs=uu_rhs, rhs=vv),
        lhrh(lhs=vv_rhs, rhs=wave_speed**2*(uu_dDD[0][0] + uu_dDD[1][1]))
    ]

k = sympify(pi/20)
toff = sympify(pi/2)
sq2 = sympify(2**.5)

init_eqns = [
    lhrh(lhs=vv, rhs=sympify(0)),
    lhrh(lhs=uu, rhs=sin(k*(x))*sin(k*(y)))
]

# What is omega?
# omega / k = wave_speed

anal_eqns = [
    lhrh(lhs=ana, rhs=sin(k*x)*sin(k*y)*sin(cctk_time*sq2*wave_speed*k+toff)-uu)
]

# access a variable with a different centering using interpolation
# looping cell-centered, access vertex-centered, not vice-versa
# all rhs variables should have the same centering
# wave toy with fluxes, fluxes are faces
# schedule something in post-regrid, apply bc's
thorn.add_func("wave_init",
    body=init_eqns,
    where='everywhere',
    schedule_bin='initial',
    doc='Do the wave init',
    centering=centering)

#thorn.add_func("wave_bound",
#    body=bound_eqns,
#    where='boundary',
#    schedule_bin='RHS after wave_evol',
#    doc='Do the b/c',
#    centering=centering)

thorn.add_func("wave_evol",
    body=evol_eqns,
    where='interior',
    schedule_bin='RHS',
    doc='Do the wave evol',
    centering=centering)

thorn.add_func("wave_anal",
    body=anal_eqns,
    where='everywhere',
    schedule_bin='Analysis',
    doc='Check the result',
    centering=centering)

assert "CACTUS_HOME" in os.environ, "Please set the CACTUS_HOME variable to point to your Cactus installation"
cactus_home = os.environ["CACTUS_HOME"]
cactus_sim = os.environ.get("CACTUS_SIM","sim")
cactus_thornlist = os.environ.get("CACTUS_THORNLIST", None)

#cactus_home = "/project/sbrandt/release/Cactus"
thorn.generate(cactus_home,cactus_config=cactus_sim,cactus_thornlist=cactus_thornlist)
