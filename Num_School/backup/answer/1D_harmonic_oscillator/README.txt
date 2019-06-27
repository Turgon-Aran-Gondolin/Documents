1)
generator the configurations:
gen_config.cc
compile: g++ gen_config.cc -o gen_config
(if you get error message compiling include<random>, try add option -std=c++11)
run: gen_config < input_gen_config
input_gen_config contains five parameters:
aomega: the value of a*omega, a is the lattice spacing, omega is the angular frequency. More information in the code. 
Nt: number of lattice sites. 
Nconf: number of configuration
delta: step size 
output_path: the directory of the output file, please make sure the directory exists. Otherwise you get the error information "fail to open output file". 

2)
Calculate the expectation value of x^2:
measure_x2.cc
compile: g++ measure_x2.cc -o measure_x2
run: measure_x2 < input_measure
input_measure contains 8 parameters:
aomega: same as above
Nt: same as above
Nconf: same as above
Nskip: number of configurations to be skipped before making each measurement, to avoid autocorrelation between configurations. 
Ndump: skip the first Ndump configurations before starting measurement.
delta: same as above
input_path: the directory of configuration file.
output_path: the directory for output file. The output file contains values of x^2 for each measurements. mean value and error will be printed out on the screen. 

3) 
Calculate the correlation function C(t)= <x(t+t0)x(t0)>
measure_corr.cc
compile: g++ measure_corr.cc -o measure_corr
run: measure_corr < inpupt_measure
the input parameters are the same as measure_x2
  

