# Metapopulations
Simulation of disease systems consisting of spatially separated meta-populations.

The base code for the simulation is in metapop.cpp. This code takes the example command line

./metapop beta 8.0 maxvaccprob 0.84 minvaccprob 0.1 timesteps 480 patchpop 10000 birthrate 40.0 popstddev 0.0 images 1 logs 1 fn pop10kdev0.1 iter 10 stochastic 0

With inputs
* beta: infection parameter. Usually we have a basic reproduction rate, R0, in mind, and the infection period, gamma, is held fixed in this code at one time-step (2 weeks.) So, beta is defined
         via beta/gamma = R0, with the units in weeks.
* maxvaccprob: maximum vaccination probability
* minvaccprob: the standard deviation in vaccination probability (confusing, I know. It's named this way because of a kludge.)
* timesteps: number of timesteps to run simulation
* patchpop: average population of each patch
* birthrate: birth rate
* popstddev: standard deviation of patch population
* images: print out images?
* logs: print out logs?
* fn: filename for the stem
* iter: number of iterations to run
* stochastic: yes/no - if no, each of the iterations should be identical.

To run over a grid of possible input parameters, you can use the script gridsearch.py. To run with different mixing matrices, you need to change metapop.cpp by hand.
Make sure this is included in the beginning stuff:

#define loadbmat 1

Then, on line 607 you see that it looks to load in the mixing matrix. This is a file defined on line 33. For instance:

#include "alpha10c0.1/betamalpha10c0.1.hpp"

 These different matrices are generated using MixingMatrixGenerator.ipynb.

