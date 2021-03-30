# Introduction

The purpose of the simulation was to compare the bias in indirect or mixed direct/indirect treatment effect estimation made by an aggregated-data (AD) based approach vs an Individual Patients Data (IPD)-based approach in a simple 3 nodes network.


# Simulation parameters

The simulation parameters were as follows:

Two configurations :

* Network without a closed loop
* Network with a closed loop

Two treatment effects (log HR):

* -0.2
* -0.5

Two random effects were simulated: one for the trial specific baseline risk and one for the treatement effect. Both with variance :
* 0.1
* 0.01

Three scenarios :

* None: same age distribution and no interaction
* Interaction: interaction in one comparison not the other but same age distribution 
* Both: interaction in one comparison not the other and age distribution differences between comparison

# Simulation file

There is a lot of simulation so the simulation that will usually excess the RAM, so the data are splits in chunk. 

The code is parallelized by furrr

sim_parallel_3a : generate the simulation data (allowing to divide the data in chunks before processing)
sim_parallel_3b : perform the analysis (chunk by chunk)
sim_parrallel_3c : encapsulate all the analysis into a single function so that it can be efficiently run in parallel
fonctions : the fonctions used.
