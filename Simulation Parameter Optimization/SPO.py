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
import os
import subprocess
from scipy.optimize import minimize
class SPOStatus(Enum):
    STARTING = 0
    WAITING = 1
    READY = 2
    FINISHED = 3

class SPORunsOn(Enum):
    Desktop = 0
    HPCC = 1


class SimulationParameterOptimizer:
    def __init__(self,configFile,ensembleNo = None):
        configParser = SPOFileParser(configFile)
        (self.name, self.simulationCommand, self.parameters, self.dataSpec, self.runsOn, self.method) = configParser.parseConfigFile()

        self.logFileName = self.name + "Log.txt"
        self.optimizerLogFile = self.name + "OptimizerLog.txt"

        #now add things like the iteration number to the name 
        if ensembleNo is not None:
            self.ensembleNo = ensembleNo

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
                self.createLogFile()

                #start a run
                self.setupAndStartRun()

                #and exit
                return
            
            case SPOStatus.WAITING: 
                # a simulation run is finished but not all the simulation in the
                # ensemble have finished. We should just exit
                return

            case SPOStatus.READY:
                # all the simulation in this ensemble have finished, 
                # run an iteration of the optimizer and start the next simulation
                fitness = self.compareData()

                # if we're below tolerance then finish up
                if fitness<self.tol:
                    self.finishUp()

                else:
                    # otherwise update the parameters
                    # first read in the parameters and residuals
                    (parameters,residuals) = self.readLog()
                    self.updateParameters(parameters,residuals)

                    # and start the next run
                    self.setupAndStartRun()

    def writeLogFileHeader(self):
        logFile = open(self.logFileName,'a')
        logFile.write("# "+self.name+"\n")
        logFile.write("########################\n")
        logFile.write("Command: " + self.simulationCommand+"\n")
        logFile.write("Target Data: ")
        for spec in self.dataSpec:
            logFile.write(spec[1])
        logFile.write("\n")
        logFile.close()

        self.step = 0

    def setupAndStartRun(self):
        #makes the folder name/step
        os.mkdir(self.name)
        os.mkdir(self.name + "/" + self.step)

        myRunner = SPOSimulationRunner(self.simulationCommand,self.parameters,self.runsOn)
        myRunner.createScript()
        myRunner.callScript()

    def compareData(self):
        pass
    def updateParameters(self,parameters,residuals):
        # read the log to get the list of parameters and residuals
        myOptimizer = SPOOptimizer(self.method,self.maxSteps)
        newParamList = myOptimizer.get_next_parameters(parameters,residuals)
        
        for i in range(len(self.parameters)):
            self.parameters[i][1] = newParamList[i]

    def finishUp(self):
        pass
    
    # Check the log file to determine status of the runs.
    def checkStatus(self):
        logParser = SPOFileParser(self.logFile)
        logData = logParser.parseLog()
        
        if logData is None:
            return SPOStatus.STARTING
        pass
    


class SPOFileParser:
    def __init__(self,fileName):
        self.file = open(fileName,"r")
        #keep a list of line positions so we can rewind
        self.lastLinePos = []
        self.lastLinePos.append(self.file.tell())

    def __del__(self):
        self.file.close()

    def parseConfigFile(self):
        #grab the name
        # header =
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
        dataSpec = []
        while True:
            line = self._nextLine()
            fileNames = re.match(line,"(.*)\s+(.*)")
            if fileNames:
                entry = (fileNames.group(1),fileNames.group(2),self._nextLine)
                dataSpec.append(entry)
            else:
                # we didn't find a file name, check that this is a section header and then rewind and return
                sectionHeader = re.match(".*:")
                if sectionHeader:
                    self._rewind()
                    break
                else:
                    raise Exception("Failed to read data section.")
        return dataSpec

    
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
        self.lastLinePos.append(self.file.tell())
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
    def _rewind(self):
        self.file.seek(self.lastLinePos.pop())


class SPOSimulationRunner:
    def __init__ (self,SPO):
        self.command = SPO.command
        self.parameters = SPO.parameters
        self.runsOn = SPO.runsOn
        self.name = SPO.name
        if self.runsOn[0] == SPORunsOn.HPCC:
            self.partition = SPO.partition
        self.maxJobs = SPO.maxJobs
        self.scriptName = "scriptRunner.sh"

    def setupFolder(self):
        pass
    
    def createCommand(self):
        commandWithParams = self.command
        for param in self.parameters.keys():
            commandWithParams = re.sub("(\w)"+param+"(\w?)","$1" + self.parameters[param]+ "$2",commandWithParams)
        return commandWithParams
    
    def writeScript(self):
        file = open(self.scriptName,"w")

        file.write("#!/usr/bin/bash\n")
        match self.runsOn[0]:
            case SPORunsOn.Desktop:
                self.writeDesktopScript(file)
            case SPORunsOn.HPCC:
                self.writeHPCCScript(file,self.runsOn[1])
        file.close()

    def writeDesktopScript(self,file):
        file.write(self.createCommand()+"\n")
        file.write("python3"+str(__file__))
    
    def writeHPCCScript(self,file,ensembleSize):
        # for not using an ensemble
        file.write("#SBatch --job-name="+self.name+"\n")
        file.write("#SBatch --partition="+self.partition+"\n")

        if ensembleSize>1:
            self.writeDesktopScript(file)
        else:
            #use a job array to submit the ensemble, 
            file.write("#SBatch --array=1-"+ensembleSize+"%"+self.maxJobs+"\n")
            self.writeDesktopScript(file)
            # after the simulation finished call this script again with the 
            # ensemble ID
            
            file.wrire(" $SLURM_ARRAY_TASK_ID")

    def runScript(self):
        runString = ""
        if self.runsOn[0] == SPORunsOn.HPCC:
            runString = "sbatch "+ self.scriptName + " &"
        else:
            runString = "./"+ self.scriptName + " &"

        subprocess.run(runString,shell=True)

class SPOOptimizer:
    def __init__(self,method,maxSteps):
        self.method = method
        self.step = 0
        self.maxSteps = maxSteps
    def get_next_parameters(self,parameters,residual):
        minimize(self.pastValues,parameters[0],args=(parameters,residual),
                method = self.method, options={"maxiter":self.maxSteps})
        return self.newParam

    def pastValues(self,newParam,parameters,residual):
        # if this is a new parameter
        if self.step>len(parameters):
            # save the values
            self.newParam = newParam
            # return 0 to stop the optimization
            return 0
        
        # if this is not a new parameter make sure our parameter trajectory 
        # is valid
        if newParam != parameters[self.step]:
            raise Exception("ERROR: encountered Nondeterministic solver")
        # and return the residual
        retVal = residual[self.step]
        self.step+=1
        return retVal


