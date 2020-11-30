# Introduction

The purpose of the simulation was to evaluate AD vs IPD in a simple 3 nodes network.


# Simulation parameters

THe simulation parameters were as follows:

Two configurations :

* Network without a closed loop
* network with a closed loop

Two treatment effects :

* -0.2
* -0.5

Two random components variation :
* 0.1
* 0.01

Three scenarios :

* None
* Interaction
* Both

# Simulation file

sim_parallel_3a : generate the data allowing to divide the data in chunks before processing
sim_parallel_3b : perform the analysis
sim_parrallel_3c : encapsulate the analysis model into a single function so thaht it can be run in parallel
fonctions : the fonctions used.