#!/bin/bash

# Parameters for "harm_osc.c"
seed=82364
num_traj=50000
#num_traj=15
init=0

M_0=1.0
mu2=1.00
la=0.0
#a=0.10
a=1.0
NT=200
#NT=10

Nskip=5000
#Nskip=5

# GENERATION OF CONFIGURATIONS:

# Run the MAIN program "harm_osc.c"

time ./harm_osc $NT $seed $num_traj $init $a $M_0 $mu2 $la $Nskip 

