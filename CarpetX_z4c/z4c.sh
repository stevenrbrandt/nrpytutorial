#!/bin/bash

# Sample invocation
# Edit the environ variables below before running.
export CACTUS_HOME=/project/sbrandt/cactusamrex/Cactus
export CACTUS_THORNLIST=${CACTUS_HOME}/../carpetx.th
export CACTUS_DRIVER=CarpetX
export CACTUS_SIM=sim
export SIMFAC_OPTS="--machine db-sing-nv"

time python3 -c 'from CarpetX_z4c.Z4c import main; main()' |& tee z4.out
