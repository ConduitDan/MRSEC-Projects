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


#TODO refactor the parser
#TODO what about integer parameters
#TODO change comment character to // (so we can have slurm commands)


from enum import Enum
from logging import exception
import re
import os
import sys
import subprocess
import fileinput
import pandas as pd
import numpy as np
import SPOParser


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
        config = SPOParser.FileSpecFactory().configFileSpec()
        configValues = config.readFile(configFile)
        self.configFile = configFile
        self.configuration = configValues
        # self.name = configValues["Name:"]
        # self.simulationCommand = configValues["Simulation:"]
        # self.parameters = configValues["Parameters:"]
        # self.dataSpec = configValues["Data:"]
        # self.runsOn = configValues["Runs On:"]
        # self.partition = configValues["Partition:"]
        # self.extraCommands = configValues["Extra Script Commands:"]
        # self.method = configValues["Method:"]

        self.logFileName = self.configuration["Name:"] +"/" + self.configuration["Name:"] + "Log.txt"
        self.optimizerLogFile = self.configuration["Name:"] + "OptimizerLog.txt"

        #now add things like the iteration number to the name 
        self.ensembleNo = ensembleNo


        self.maxJobs = 100
        self.tol = 10
        self.maxSteps = 100

    def run(self):
        # the main run function
        
        # first check the log file, this tells us if we're just starting the 
        # process, if we're waiting for other nodes to finish, if we're ready 
        # for the next iteration, or if we're done with the process
        
        status = self.checkStatus()
        
        self.path = self.configuration["Name:"]+"/parameter_step_"+str(self.step)

        if status == SPOStatus.STARTING:
            #We're just starting up
            print("Starting Optimization")
            #create the log file
            self.writeLogFileHeader()
            self.writeStartRunLog()

            #start a run
            self.setupAndStartRun()

            #and exit
            return
        
        elif status == SPOStatus.WAITING: 
            # a simulation run is finished but not all the simulation in the
            # ensemble have finished. We should just exit
            return

        elif status == SPOStatus.READY:
            print("Step %d finished"%self.step)
            # all the simulation in this ensemble have finished, 
            # run an iteration of the optimizer and start the next simulation
            residue = self.compareData()
            self.step += 1
            self.writeResidueLog(residue)
            self.residuals.append(residue)

            # if we're below tolerance then finish up
            if residue<self.tol:
                print("Found solution set %s, residue is %s"\
                    %(str(self.parameterHistory[-1]),residue))
                self.finishUp()

            else:
                # otherwise update the parameters
                # first read in the parameters and residuals
                
                self.updateParameters(self.parameterHistory,self.residuals)
                print("Running step %d"%self.step)
                self.writeStartRunLog()

                # and start the next run
                self.setupAndStartRun()
    def writeStartRunLog(self):
        logFile = open(self.logFileName,'a')
        logFile.write("Step "+str(self.step)+"    {")
        for param in self.configuration["Parameters:"]:
            logFile.write(param[0]+":"+str(param[1]))
            if param != self.configuration["Parameters:"][-1]:
                logFile.write(",")
        logFile.write("}    ")
        logFile.close()

    def writeResidueLog(self,residue):
        logFile = open(self.logFileName,'a')
        logFile.write("Residue:"+str(residue)+"\n")

    def writeLogFileHeader(self):
        os.mkdir(self.configuration["Name:"])
        logFile = open(self.logFileName,'a')
        logFile.write("# "+self.configuration["Name:"]+"\n")
        logFile.write("########################\n")
        logFile.write("# Command: " + self.configuration["Simulation:"]+"\n")
        logFile.write("# Target Data: ")
        for spec in self.configuration["Data:"]:
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
        for (dataFile,targetFile) in self.configuration["Data:"]:
            simulationData = self.getSimData(dataFile)
            targetData = self.loadData(targetFile)
            residue += self.calculateError(simulationData,targetData)
        residue /= len(self.configuration["Data:"])
        return residue

    def getSimData(self,dataFile):
        try:
            if self.configuration["Runs On:"][0] == SPORunsOn.HPCC and self.configuration["Runs On:"][1]>1:
                    simulationData = self.loadData(self.path+"/0/"+dataFile)
                    for i in range(1,self.configuration["Runs On:"][1]):
                        simulationData += self.loadData(self.path+'/'+str(i) + '/'+dataFile)
                    simulationData /= self.configuration["Runs On:"][1]
            else:
                simulationData = self.loadData(self.path+"/"+dataFile)
        except FileNotFoundError:
            # if we don't find the error we should fall back to running the simulation again.
            simulationData = None

        return simulationData

    def loadData(self,dataFile):

        data = np.loadtxt(dataFile)
            
        if data.shape == (0,0):
            raise Exception("coudln't read data")


        return data
    
    def calculateError(self,simulationData,targetData):
        # ensure that data sizes are the same
        if simulationData.shape!=targetData.shape:
            raise Exception("Data File and Target File have different formats")
        
        return np.linalg.norm(simulationData-targetData)
            

    def updateParameters(self,parameters,residuals):
        # read the log to get the list of parameters and residuals
        myOptimizer = SPOOptimizer(self.configuration["Method:"],self.maxSteps,self.step)
        newParamList = myOptimizer.get_next_parameters(parameters,residuals)
        
        for i in range(len(self.configuration["Parameters:"])):
            self.configuration["Parameters:"][i][1] = newParamList[i]

    def finishUp(self):
        print("Finishing up")

    
    # Check the log file to determine status of the runs.
    def checkStatus(self):
        # if we're not using an ensemble this we are starting if there is no log
        # or we are ready to go if the log exists
        if os.path.exists(self.logFileName):

            log = SPOParser.FileSpecFactory().logFileSpec()
            vals = log.readFile(self.logFileName)
            linesWithResidue = vals[0]
            linesWithOutResidue = vals[1]
            params = list(val[0] for val in linesWithResidue)
            params.append(linesWithOutResidue[0])

            residues = list(val[1] for val in linesWithResidue)


            self.step = linesWithOutResidue[1]
            self.parameterHistory = params
            self.residuals = residues

        
        else:
            self.step = 0
            return SPOStatus.STARTING

        if self.ensembleNo is None:
            
                return SPOStatus.READY
        else:
            # if we are using an ensemble we should read the ensemble log,
            # mark that we are done and check if everyone else is done
            logParser = SPOParser.SPOEnsembleLogParser("ensembleLog.txt")
            finished = logParser.parseLog()
            if finished:
                return SPOStatus.READY
            else:
                return SPOStatus.WAITING



class SPOSimulationRunner:
    def __init__ (self,SPO):
        self.configuration = SPO.configuration
        # self.command = SPO.simulationCommand
        (self.caller,self.sim,self.cmdLineArgs) = self.splitCommand()
        # self.parameters = SPO.parameters
        # self.runsOn = SPO.runsOn
        # self.name = SPO.name
        self.step = SPO.step
        self.configFile = SPO.configFile
        # self.configuration["Runs On:"][1] = 1
        # if self.runsOn[0] == SPORunsOn.HPCC:
            # self.partition = SPO.partition
            # self.configuration["Runs On:"][1] = self.runsOn[1]
            # self.extraCommands = SPO.extraCommands
        self.maxJobs = SPO.maxJobs
        self.scriptName = "scriptRunner.sh"
    def splitCommand(self):
        #first try and match ./
        bashCommand = re.match("./(.*)",self.configuration["Simulation:"])
        caller = ""
        command = ""
        if bashCommand:
            caller = "./"
            command = bashCommand.group(1)
        else:
            callerMatch = re.match("(.*?\s+(?:-\S+\s+)*)(.*)",self.configuration["Simulation:"])
            if callerMatch:
                caller = callerMatch.group(1)
                command = callerMatch.group(2)
            else:
                raise Exception("command not in a recognizable format")
        sim = command.split(' ')[0]
        cmdLineArgs = command.replace(sim+ ' ','',1)
        #sim = os.path.abspath(sim).replace(' ','\ ')
        return (caller,sim,cmdLineArgs)


    def createFolders(self):
        self.path = self.configuration["Name:"]+"/parameter_step_"+str(self.step)
        os.mkdir(self.path)
        if self.configuration["Runs On:"][1]>1:
            self.setupEnsembleFolders()

    def createLogs(self):
        self.write_parameters()
        if self.configuration["Runs On:"][1]>1:
            self.setupEnsembleLog()

    def write_parameters(self):
        paramLog = open(self.path+"/paramlog.txt","w")
        paramLog.write(str(self.configuration["Parameters:"]))
        paramLog.close()

    def setupEnsembleLog(self):

        ensembleLog = open(self.path+"/ensembleLog.txt",'w')
        for i in range(self.configuration["Runs On:"][1]):
            ensembleLog.write("Run "+str(i) + " of "+ str(self.configuration["Runs On:"][1]) + " Running\n")
        ensembleLog.close()

    def setupEnsembleFolders(self):
        for i in range(self.configuration["Runs On:"][1]):
            os.mkdir(self.path+'/' +str(i))
    
    def createCommand(self):
        commandWithParams = self.cmdLineArgs
        commandWithParams = " "+commandWithParams+" "
        for param in self.configuration["Parameters:"]:
            commandWithParams = re.sub(" "+param[0]+" "," "+str(param[1])+" ",commandWithParams)
        return self.caller + self.sim + ' ' + commandWithParams+"\n"
    
    def createScript(self):
        file = open(self.path+"/"+self.scriptName,"w")

        file.write("#!/usr/bin/bash\n")
        # don't know why these enums can't be equal
        # so we're just gonna use the values
        if self.configuration["Runs On:"][0].value == SPORunsOn.Desktop.value:
            self.writeDesktopScript(file)
        elif self.configuration["Runs On:"][0].value is SPORunsOn.HPCC.value:
            self.writeHPCCScript(file)
        else:
            raise Exception("Unrecognized Runs On option.")

        file.close()
        os.chmod(self.path+"/"+self.scriptName,0o755)

    def writeDesktopScript(self,file):
        file.write("cd " + self.path+"\n")
        file.write(self.createCommand())
        file.write("cd ../.. \n")
        file.write("python3 "+str(__file__).replace(" ",r'\ ') + " " + self.configFile+"\n")
        

    
    def writeHPCCScript(self,file):
        # for not using an ensemble
        file.write("#SBATCH --job-name="+self.configuration["Name:"]+"\n")
        file.write("#SBATCH --partition="+self.configuration["Partition:"]+"\n")

        if not self.configuration["Runs On:"][1]>1:
            for line in self.configuration["Extra Commands:"]:
                file.write(line+'\n')

            self.writeDesktopScript(file)
        else:
            #use a job array to submit the ensemble, 

            #file.write("#SBATCH --chdir " +self.path+"/$SLURM_ARRAY_TASK_ID\n")
            file.write("#SBATCH --array=0-"+str(self.ensembleSize-1)+"%"+str(self.maxJobs)+"\n")
            for line in self.extraCommands:
                file.write(line+'\n')
            file.write("cd " + self.path+"/$SLURM_ARRAY_TASK_ID\n")
            file.write(self.createCommand())
            file.write("cd ../../..\n")
            # after the simulation finished call this script again with the 
            # ensemble ID
            file.write("python3 "+str(__file__).replace(" ",r'\ ') + " " + self.configFile+" $SLURM_ARRAY_TASK_ID\n")

    def runScript(self):
        runString = ""
        if self.configuration["Runs On:"][0] == SPORunsOn.HPCC:
            runString = "sbatch "+self.path+"/"+ self.scriptName + " &"
        else:
            runString = self.path+"/"+ self.scriptName + " &"

        subprocess.run(runString,shell=True)

class SPOOptimizer:
    def __init__(self,method,maxSteps,currentStep):
        self.method = method
        self.step = 0
        self.maxSteps = maxSteps
        self.currentStep = currentStep

    def get_next_parameters(self,parameters,residual):
        paramNumbers =[]
        for param in parameters:
            paramNumbers.append([x[1] for x in param])
        try:
            outPut = minimize(self.pastValues,paramNumbers[0],args=(paramNumbers,residual),
                method = self.method, options={"maxiter":self.maxSteps,"maxfun":self.currentStep+1,"eps":0.05})
            if not hasattr(self,'newParam'):
                print(outPut)
        except StopIteration:
            print(outPut)    
        return self.newParam

    def pastValues(self,newParam,parameters,residual):
        # if this is a new parameter
        if self.step>=len(parameters):
            # save the values
            self.newParam = newParam
            # return 0 to stop the optimization
            raise StopIteration

        # if this is not a new parameter make sure our parameter trajectory 
        # is valid
        if (abs(newParam - parameters[self.step])>1e-6).any():
            print("encounted error on step %d"%self.step)
            print("List of params is %s"%str(parameters))
            print("I tried to use %s"%str(newParam))

            print(newParam-parameters[self.step])
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
