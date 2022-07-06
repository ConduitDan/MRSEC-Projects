"""
A Simulation Parameter Optimizer
by Danny Goldstein 

This script will try and find an optimal set of parameters for a simulation given an experiment (target) data set.
This script accepts a configuration file that details how to run the simulation, the initial parameter values, and the specifications for the data files. 
Additionally the config file will dictate if the is running on the hpcc or not and how many runs are in a needed to create an ensemble measurement for comparison.

This script will set up folders for data and run simulation, creating scripts for the scheduler as needed. 
This script does NO post processing on data, that responsibly lies with the simulation. The output of the simulation should be data that is directily comparable with the target data set.
"""

from enum import Enum

class SPOStatus(Enum):
    STARTING = 0
    WAITING = 1
    READY = 2


class SimulationParameterOptimizer:
    def __init__(self,configFile):
        pass
    
    # Check the log file to determine status of the runs.
    def checkStatus(self):
        logData = self.parseLog(self.logFile)
        
        if logData is None:
            return SPOStatus.STARTING
        



class SPOFileParser:
    def __init__(self,configFile):
        pass
    def parseData(self,dataFile):
        pass
    def parseLog(self,logFile):
        logData = None
        if
        pass

class SPOSimulationRunner:
    def __init__ (self):
        pass




