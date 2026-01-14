Here we calculate the pore-size of gel structure. 

The code is written for unsheared configuration. 

One has to change it for sheared configuration. 

Below the details are given, for an example one can run the code for restart_data.

code : poresize_v2.c

compilation: gcc poresize_v2.c -lm

exec: ./a.out filename

output: hist_D.dat where the histogram of the pore-diameter will be saved.

example:  ./a.out gel_restart.30000000.data

Getting the histogram one can find average pore size, distribution, etc.
