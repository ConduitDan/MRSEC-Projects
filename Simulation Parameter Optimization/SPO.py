"""
A Simulation Parameter Optimizer
by Danny Goldstein 

This script will try and find an optimal set of parameters for a simulation given an experiment (target) data set.
This script accepts a configuration file that details how to run the simulation, the initial parameter values, and the specifications for the data files. 
Additionally the config file will dictate if the is running on the hpcc or not and how many runs are in a needed to create an ensemble measurement for comparison.

This script will set up folders for data and run simulation, creating scripts for the scheduler as needed. 
This script does NO post processing on data, that responsibly lies with the simulation. The output of the simulation should be data that is directily comparable with the target data set.

"""

class 