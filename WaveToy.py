import os
import grid
from cactusthorn import CactusThorn
from sympy import sympify, cos

# Current options are Carpet and CarpetX
grid.ET_driver = os.environ.get("CACTUS_DRIVER","Carpet")

thorn = CactusThorn("TestOne","WaveToyNRPy2")

FD_order = thorn.declare_param('FD_order',default=4,vmin=2,vmax=8,doc="The finite difference order")

wave_speed = thorn.declare_param('wave_speed',default=1,vmin=.1,vmax=100,doc="The speed of the wave")
x0 = thorn.declare_param('x0',default=0,vmin=-100,vmax=100,doc="The x pos of the wave")
y0 = thorn.declare_param('y0',default=0,vmin=-100,vmax=100,doc="The y pos of the wave")
z0 = thorn.declare_param('z0',default=0,vmin=-100,vmax=100,doc="The z pos of the wave")
zero = thorn.declare_param('zero',default=0,vmin=0,vmax=0,doc="zero")

# AUXEVOL needed for the evo, can be freed after evaluating rhs (1 time level)
# AUX uu_rhs (1 time level)
# EVOL evolved gfs (3 time levels)
uu_rhs, vv_rhs = thorn.register_gridfunctions("AUX", ["uu_rhs", "vv_rhs"], centering="CCC")
uu, vv = thorn.register_gridfunctions("EVOL", ["uu", "vv"], centering="CCC")
x,y,z = thorn.get_xyz()

from outputC import lhrh
import indexedexp as ixp

uu_dDD = ixp.declarerank2("uu_dDD","sym01")

evol_eqns = [
    lhrh(lhs=uu_rhs, rhs=vv),
    lhrh(lhs=vv_rhs, rhs=wave_speed**2*(uu_dDD[0][0] + uu_dDD[1][1] + uu_dDD[2][2]))
]

init_eqns = [
    lhrh(lhs=vv, rhs=sympify(0)),
    lhrh(lhs=uu, rhs=cos(-(x-x0)**2-(y-y0)**2-(z-z0)**2)),
]

thorn.add_func("wave_init", body=init_eqns, schedule_bin='initial', doc='Do the wave init')
thorn.add_func("wave_evol", body=evol_eqns, schedule_bin='evol', doc='Do the wave evol')

assert "CACTUS_HOME" in os.environ, "Please set the CACTUS_HOME variable to point to your Cactus installation"
cactus_home = os.environ["CACTUS_HOME"]
cactus_sim = os.environ.get("CACTUS_SIM","sim")
cactus_thornlist = os.environ.get("CACTUS_THORNLIST", None)

#cactus_home = "/project/sbrandt/release/Cactus"
thorn.generate(cactus_home,cactus_config=cactus_sim,cactus_thornlist=cactus_thornlist)
