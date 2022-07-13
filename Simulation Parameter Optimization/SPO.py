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


#TODO Juggle the cwd better

from enum import Enum
from logging import exception
import re
import os
import sys
import subprocess
import fileinput
import pandas as pd
import numpy as np


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
        self.configFile = configFile
        (self.name, self.simulationCommand, self.parameters, self.dataSpec, self.runsOn, self.method) = configParser.parseConfigFile()

        self.logFileName = self.name +"/" + self.name + "Log.txt"
        self.optimizerLogFile = self.name + "OptimizerLog.txt"

        #now add things like the iteration number to the name 
        if ensembleNo is not None:
            self.ensembleNo = ensembleNo


        self.maxJobs = 100

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
                self.writeLogFileHeader()
                self.writeStartRunLog()

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
                residue = self.compareData()
                self.writeResidueLog(residue)

                # if we're below tolerance then finish up
                if residue<self.tol:
                    self.finishUp()

                else:
                    # otherwise update the parameters
                    # first read in the parameters and residuals
                    (parameters,residuals) = self.readLog()
                    self.updateParameters(parameters,residuals)
                    self.writeStartRunLog()

                    # and start the next run
                    self.setupAndStartRun()
    def writeStartRunLog(self):
        logFile = open(self.logFileName,'a')
        logFile.write("Step "+str(self.step)+"    {")
        for param in self.parameters:
            logFile.write(param[0]+":"+str(param[1]))
            if param != self.parameters[-1]:
                logFile.write(",")
        logFile.write("}    ")
        logFile.close()

    def writeResidueLog(self,residue):
        logFile = open(self.logFileName,'a')
        logFile.write("Residue:"+str(residue)+"\n")

    def writeLogFileHeader(self):
        os.mkdir(self.name)
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
        myRunner = SPOSimulationRunner(self)
        myRunner.createFolders()
        myRunner.createLogs()
        myRunner.createScript()
        myRunner.runScript()

    def compareData(self):
        residue = 0
        for (dataFile,targetFile) in self.dataSpec:
            simulationData = self.loadData(dataFile)
            targetData = self.loadData(targetFile)
            residue += self.calculateError(simulationData,targetData)
        residue /= len(self.dataSpec)

    def loadData(self,dataFile):
        data = pd.read_csv(dataFile)
        return data.DataFrame.to_numpy()
    
    def calculateError(self,simulationData,targetData):
        # ensure that data sizes are the same
        if simulationData.shape()!=targetData.shape():
            raise Exception("Data File and Target File have different formats")
        
        return np.linalg.norm(simulationData-targetData)
            

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
        # if we're not using an ensemble this we are starting if there is no log
        # or we are ready to go if the log exists
        if os.path.exists(self.logFileName):
            logParser = SPOFileParser(self.logFileName)
            (self.parameters,self.residuals,self.step) = logParser.parseLogFile()

        
        else:
            self.step = 0
            return SPOStatus.STARTING

        if self.ensembleNo is None:
            
                return SPOStatus.READY
        else:
            # if we are using an ensemble we should read the ensemble log,
            # mark that we are done and check if everyone else is done
            logParser = SPOEnsembleLogParser("ensembleLog.txt")
            finished = logParser.parseLog()
            if finished:
                return SPOStatus.READY
            else:
                return SPOStatus.WAITING
        
        

class SPOEnsembleLogParser:
    def __init__(self,fileName):
        self.file = fileinput.input(fileName,inplace=True)

    def __del__(self):
        self.file.close()
    
    def parseEnsembleLog(self,ensembleID):
        finished = True
        for line in fileinput.input("File3.txt", inplace=True):
                        #look for lines that end in "Running " 
            runningJob = re.search("Run (%d+) *. Running ",line)
            if runningJob:
                #we found a job that is still running
                if int(runningJob.group(1)) == ensembleID:
                    # we are the responsible of that job,
                    # update that is is finished
                    line = line.replace("Running","Finished")
                    sys.stdout.write(line)

                else:
                    finished = False
        return finished

class SPOFileParser:
    def __init__(self,fileName):
        self.file = open(fileName,"r+")
        self.fileName = fileName
        #keep a list of line positions so we can rewind
        self.lastLinePos = []
        self.lastLinePos.append(self.file.tell())
        self.floatMatch = "([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"



    def __del__(self):
        self.file.close()

    def parseLogFile(self):
        step = 0
        params = []
        residues = []
        for line in self:
            logData = re.match("Run\w+(\d*\.?\d+)\w+{(.*)}\w+E:"+self.floatMatch,line)
            if logData:
                step = int(logData.group(1))
                paramList = []
                paramPairs = logData.group(2).split(",")
                for param in paramPairs:
                    (name,value) = param.split(':')
                    paramList.append([name,float(value)])
                params.append(paramList)
                residues.append(float(logData.group(3)))
            else:
                raise Exception("Error reading log file")
        return(params,residues,step)

    def checkHeader(self,expected):
        header = self.__next__()
        if header !=expected:
            raise Exception("Configuration File Format Error, expected '"+expected+"' but found '" + header + "'.")


    def parseConfigFile(self):
        #grab the name
        self.checkHeader("Name:")
        name = self.__next__()
        
        #grab the simulation command
        self.checkHeader("Simulation:")
        simulationCommand = self.__next__()

        #grab the initial parameters string
        parameters = self._parseParameters()

        #advance to the data spec
        dataSpec = self._parseDataSpec()

        #grab the RunningOn option
        self.checkHeader("Runs on:")
        runsOn = self.__next__()

        self.checkHeader("Method:")
        method = self.__next__()
        #return all these to the SPO
        return (name, simulationCommand, parameters, dataSpec, runsOn, method)

    def parseData(self):
        dataSpec = []
        while True:
            line = self.__next__()
            fileNames = re.match(line,"(.*)\s+(.*)")
            if fileNames:
                entry = (fileNames.group(1),fileNames.group(2),self.__next__)
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

    def _parseDataSpec(self):
        dataSpec =[]
        for line in self:
            if line == "Runs on:":
                self._rewind()
                break
            dataSpec.append( line.split(" "))
        return dataSpec


    def _parseParameters(self):
        initialParams = []
        
        # we expect a parameter name could be anything followed by a :
        # we expect the value to be a number
        
        # this regex matches (any combination of non whitespace characters)
        # followed by possible white space then a colon followed by possibly more white space 
        # then it matches (0 or more digits 0 or 1 decimal points, at least one digit)
        for line in self:
            param = re.search("(.*):\s*"+self.floatMatch,line)
            if param:
                initialParams.append([param.group(1),float(param.group(2))])
            elif line == "Data:":
                break
            else:
                print(line)
                raise Exception("Formatting error in initial parameters")
        return initialParams

    def __next__(self):
        #grabs the next line that isn't just a new line or comment ('#')
        self.lastLinePos.append(self.file.tell())
        line = ""
        while len(line) == 0:
            line = self._trimLine(self.file.readline())
        return line
        
    def __iter__(self):
        return self
    def _trimLine(self, line):
        #look for the comment symbol
        comment = re.search("(.*)#",line)
        # if its found
        if comment:
            line = comment.group(1) # grab the first group as the line
        line = line.strip() #remove any leading or trailing whitespace
        return line
    def _rewind(self):
        self.file.seek(self.lastLinePos.pop())


class SPOSimulationRunner:
    def __init__ (self,SPO):
        self.command = SPO.simulationCommand
        self.parameters = SPO.parameters
        self.runsOn = SPO.runsOn
        self.name = SPO.name
        self.step = SPO.step
        self.configFile = SPO.configFile
        self.ensembleSize = 1
        if self.runsOn[0] == SPORunsOn.HPCC:
            self.partition = SPO.partition
            self.ensembleSize = self.runsOn[1]
        self.maxJobs = SPO.maxJobs
        self.scriptName = "scriptRunner.sh"

    def createFolders(self):
        self.path = self.name+"/parameter_step_"+str(self.step)
        os.mkdir(self.path)
        if self.ensembleSize>1:
            self.setupEnsembleFolders()

    def createLogs(self):
        self.write_parameters()
        if self.ensembleSize>1:
            self.setupEnsembleLog()

    def write_parameters(self):
        paramLog = open(self.path+"/paramlog.txt","w")
        paramLog.write(str(self.parameters))
        paramLog.close()

    def setupEnsembleLog(self):
        ensembleLog = open(self.path+"ensembleLog.txt",'w')
        for i in range(self.ensembleSize):
            ensembleLog.write("Run "+str(i) + " of "+ str(self.ensembleSize) + "Running ")

    def setupEnsembleFolders(self):
        for i in range(self.ensembleSize):
            os.mkdir(str(i))
    
    def createCommand(self):
        commandWithParams = self.command
        commandWithParams += " "
        for param in self.parameters:
            commandWithParams = re.sub(" "+param[0]+" "," "+str(param[1])+" ",commandWithParams)
        return commandWithParams
    
    def createScript(self):
        file = open(self.path+"/"+self.scriptName,"w")

        file.write("#!/usr/bin/bash\n")
        match self.runsOn:
            case "Desktop":
                self.writeDesktopScript(file)
            case SPORunsOn.HPCC:
                self.writeHPCCScript(file,self.runsOn[1])
        file.close()
        os.chmod(self.path+"/"+self.scriptName,0o755)

    def writeDesktopScript(self,file):
        file.write(self.createCommand()+"\n")
        
        file.write("python3 "+str(__file__).replace(" ",r'\ ') + " " + self.configFile)
    
    def writeHPCCScript(self,file):
        # for not using an ensemble
        file.write("#SBatch --job-name="+self.name+"\n")
        file.write("#SBatch --partition="+self.partition+"\n")

        if self.ensembleSize>1:
            self.writeDesktopScript(file)
        else:
            #use a job array to submit the ensemble, 
            file.write("#SBatch --chdir " +self.path+"/$SLURM_ARRAY_TASK_ID")
            file.write("#SBatch --array=0-"+str(self.ensembleSize-1)+"%"+self.maxJobs+"\n")
            self.writeDesktopScript(file)
            # after the simulation finished call this script again with the 
            # ensemble ID
            
            file.write(" $SLURM_ARRAY_TASK_ID")

    def runScript(self):
        runString = ""
        if self.runsOn[0] == SPORunsOn.HPCC:
            runString = "sbatch "+self.path+"/"+ self.scriptName + " &"
        else:
            runString = self.path+"/"+ self.scriptName + " &"

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

if __name__ == "__main__":
    configFile = sys.argv[1]
    ensembleNo = None
    if len(sys.argv)>2:
        ensembleNo = int(sys.argv[2])

    mySPO = SimulationParameterOptimizer(sys.argv[1],ensembleNo)
    mySPO.run()




