import os
from cactusthorn import CactusThorn

thorn = CactusThorn("TestOne","WaveToyNRPy")

FD_order = thorn.declare_param('FD_order',default=4,vmin=2,vmax=8,doc="The finite difference order")

wave_speed = thorn.declare_param('wave_speed',default=1,vmin=.1,vmax=100,doc="The speed of the wave")

# AUXEVOL needed for the evo, can be freed after evaluating rhs (1 time level)
# AUX uu_rhs (1 time level)
# EVOL evolved gfs (3 time levels)
uu_rhs, vv_rhs = thorn.register_gridfunctions("AUX", ["uu_rhs", "vv_rhs"])
uu, vv = thorn.register_gridfunctions("EVOL", ["uu", "vv"])

from outputC import lhrh
import indexedexp as ixp

uu_dDD = ixp.declarerank2("uu_dDD","sym01")

evol_eqns = [
    lhrh(lhs=uu_rhs, rhs=vv),
    lhrh(lhs=vv_rhs, rhs=(uu_dDD[0][0] + uu_dDD[1][1] + uu_dDD[2][2]))
]

thorn.add_func("wave_evol", body=evol_eqns, schedule_bin='evol', doc='Do the wave evol')

cactus_home = os.path.join(os.environ["HOME"],"cactus","Cactus")
thorn.generate(cactus_home)
