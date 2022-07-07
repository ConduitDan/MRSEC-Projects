"""
A Simulation Parameter Optimizer
by Danny Goldstein 

This script will try and find an optimal set of parameters for a simulation 
given an experiment (target) data set.

This script accepts a configuration file that details how to run the simulation,
the initial parameter values, and the specifications for the data files.

Additionally the config file will dictate if the is running on the hpcc or not
and how many runs are in a needed to create an ensemble measurement for comparison.

This script will set up folders for data and run simulation, creating scripts 
for the scheduler as needed. 

This script does NO post processing on data, that responsibly lies with the 
simulation. The output of the simulation should be data that is directily
comparable with the target data set.
"""

from enum import Enum
import re
class SPOStatus(Enum):
    STARTING = 0
    WAITING = 1
    READY = 2
    FINISHED = 3

class SimulationParameterOptimizer:
    def __init__(self,configFile,nameAppends = None):
        configParser = SPOFileParser(configFile)
        (self.name, self.simulationCommand, self.parameters, self.dataSpec, self.runsOn, self.method) = configParser.parseConfigFile()

        self.logfile = self.name + "Log.txt"

        #now add things like the iteration number to the name 
        if nameAppends is not None:
            self.name += nameAppends

    def run(self):
        # the main run function
        
        # first check the log file, this tells us if we're just starting the 
        # process, if we're waiting for other nodes to finish, if we're ready 
        # for the next iteration, or if we're done with the process
        
        status = self.checkStatus()
        match status:
            case SPOStatus.STARTING:
                #We're just starting up

                #create the log file

                pass
            
            case SPOStatus.WAITING: 
                # a simulation run is finished but not all the simulation in the
                # ensemble have finished. We should just exit
                return

            case SPOStatus.READY:
                # all the simulation in this ensemble have finished, 
                # run an iteration of the optimizer and start the next simulation
                pass

            case SPOStatus.FINISHED:
                # we have reached the end of the optimization loop, collate the
                # results, plotting if requested  
                pass

    
    # Check the log file to determine status of the runs.
    def checkStatus(self):
        logData = self.parseLog(self.logFile)
        
        if logData is None:
            return SPOStatus.STARTING


class SPOFileParser:
    def __init__(self,fileName):
        self.file = open(fileName,"r")

    def __del__(self):
        self.file.close()

    def parseConfigFile(self):
        #grab the name
        assert(self._nextLine()=="name:","Configuration File Format Error, Config file must start with 'Name:'")
        name = self._nextLine()
        
        #grab the simulation command
        assert(self._nextLine()=="Simulation:","Configuration File Format Error, Couldn't find expected 'Simulation:'")
        simulationCommand = self._nextLine()

        #grab the initial parameters string
        parameters = self._parseParameters()

        #advance to the data spec
        dataSpec = self._parseDataSpec()

        #grab the RunningOn option
        assert(self._nextLine()=="Runs on:","Configuration File Format Error, didn't find 'Runs on:'")
        runsOn = self._nextLine()

        assert(self._nextLine()=="Method:","Configuration File Format Error, didn't find 'Method:'")
        method = self._nextLine()
        #return all these to the SPO
        return (name, simulationCommand, parameters, dataSpec, runsOn, method)

    def parseData(self):
        pass
    def parseLog(self):
        logData = None
        pass
    def _parseParameters(self):
        initialParams = {}
        
        # we expect a parameter name could be anything followed by a :
        # we expect the value to be a number
        
        # this regex matches (any combination of non whitespace characters)
        # followed by possible white space then a colon followed by possibly more white space 
        # then it matches (0 or more digits 0 or 1 decimal points, at least one digit)
        while True:
            line = self._nextLine()
            param = re.search(line,"(.*):\s*(\d*\.?\d+")
            if param:
                initialParams[param.group(1)] = float(param.group(2))
            elif line == "Data:":
                break
            else: 
                raise Exception("Formatting error in initial parameters")
        return initialParams

    def _nextLine(self):
        #grabs the next line that isn't just a new line or comment ('#')
        line = ""
        while len(line) == 0:
            line = self._trimLine(next(self.file))
        return line

    def _trimLine(self, line):
        #look for the comment symbol
        comment = re.search(line,"(.*)#")
        # if its found
        if comment:
            line = comment.group(1) # grab the first group as the line
        
        line = line.strip() #remove any leading or trailing whitespace
        return line



class SPOSimulationRunner:
    def __init__ (self):
        pass
    def setupFolder(self):
        pass
    def writeScript(self):
        pass
    def runScript(self):
        pass

class SPOOptimizer:
    def __init__(self):
        pass
    def 




