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
import operator
from hyperopt import fmin, hp
import time



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
        self.parameters = SPOParameters(self.configuration)
        self.logFileName = self.configuration["Name:"] +"/" + self.configuration["Name:"] + "Log.txt"
        self.optimizerLogFile = self.configuration["Name:"] + "OptimizerLog.txt"

        self.configuration["EPS:"] = 1e-2
        self.configuration["Tolerance:"] = 1e-6
        self.configuration["Max Steps:"] =1000


        #now add things like the iteration number to the name 
        self.ensembleNo = ensembleNo
        self.step = 0
        self.updatePath()


        self.maxJobs = 100

    def run_with_overhead(self):
        # this runs the optimizer with the overhead of a process running in the 
        # background. This is a little bit easier as we don't have to worry 
        # about where where we are in the process

        print("Starting Optimization")
        #create the log file
        self.writeLogFileHeader()
        self.writeStartRunLog()
        myRunner = SPOSimulationRunnerWithOverhead(self)
        myOptimzer = SPOOptimizerWithOverhead(self,myRunner)
        output = myOptimzer.run()
        self.finishUp(output)

        

    def run(self):
        # the main run function
        
        # first check the log file, this tells us if we're just starting the 
        # process, if we're waiting for other nodes to finish, if we're ready 
        # for the next iteration, or if we're done with the process
        
        status = self.checkStatus()
        

        if status == SPOStatus.STARTING:
            #We're just starting up
            print("Starting Optimization")
            #create the log file
            self.writeLogFileHeader()
            self.updateParameters([])
            # self.writeStartRunLog()


            myRunner = SPOSimulationRunner(self)
            #start a run
            self.setupAndStartRun(myRunner)

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
                    %(str(self.parameters.getParamsWithLabels()[-1]),residue))
                self.finishUp()

            else:
                # otherwise update the parameters
                # first read in the parameters and residuals
                
                self.updateParameters(self.residuals)
                print("Running step %d"%self.step)
                self.writeStartRunLog()
                # make the runner
                myRunner = SPOSimulationRunner(self)

                # and start the next run
                self.setupAndStartRun(myRunner)
    def writeStartRunLog(self):
        logFile = open(self.logFileName,'a')
        logFile.write("Step "+str(self.step)+"    {")
        # latestParam  =self.parameters.getParamsWithLabels()[-1]
        latestParam  =self.parameters.initalParams
        for param in latestParam:
            logFile.write(param[0]+":"+str(param[1]))
            if param != latestParam[-1]:
                logFile.write(",")
        logFile.write("}    ")
        logFile.close()

    def writeResidueLog(self,residue):
        logFile = open(self.logFileName,'a')
        logFile.write("Residue:"+str(residue)+"\n")
        logFile.close()

    def writeLogFileHeader(self):
        try:
            os.mkdir(self.configuration["Name:"])
        except:
            pass
        logFile = open(self.logFileName,'a')
        logFile.write("# "+self.configuration["Name:"]+"\n")
        logFile.write("########################\n")
        logFile.write("# Command: " + self.configuration["Simulation:"]+"\n")
        logFile.write("# Target Data: ")
        for spec in self.configuration["Data:"]:
            logFile.write(spec[1])
        logFile.write("\n")
        logFile.close()

    def setupAndStartRun(self,myRunner):
        #makes the folder name/step
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
            if self.configuration["Runs On:"][0].value == SPORunsOn.HPCC.value and self.configuration["Runs On:"][1]>1:
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
            

    def updateParameters(self,residuals):
        # read the log to get the list of parameters and residuals
        myOptimizer = SPOOptimizer(self.configuration["Method:"],self.maxSteps,self.step,space = self.configuration["Parameter Space:"])
        newParamList = myOptimizer.get_next_parameters(self.parameters,residuals)
        self.parameters.addNewParam(newParamList)
        
    def finishUp(self,output = None):
        print("Finishing up")
        if output is not None:
            print(output)

    def updatePath(self):
        self.path = self.configuration["Name:"]+"/parameter_step_"+str(self.step)
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
            self.parameters.addParamHistory(params)
            self.residuals = residues
            self.updatePath()

        else:
            self.step = 0
            self.updatePath()
            return SPOStatus.STARTING

        if self.ensembleNo is None:
            
                return SPOStatus.READY
        else:
            # if we are using an ensemble we should read the ensemble log,
            # mark that we are done and check if everyone else is done
            logParser = SPOParser.SPOEnsembleLogParser(self.path + "/ensembleLog.txt")
            finished = logParser.parseEnsembleLog(self.ensembleNo)
            if finished:
                return SPOStatus.READY
            else:
                return SPOStatus.WAITING

class SPOParameters:
    def __init__(self,config):
        # params are of the form [[[label_0,value_0_0],[label_1,value_1_0],...],[[label_0,value_0_1],...]],...]
        self.initalParams = config["Parameters:"]
        self.paramSpace = config["Parameter Space:"]
        self.paramHistory = []
    def addParamHistory(self,paramHistory):
        self.paramHistory = paramHistory
    def getInitalConditions(self):
        # return self.removeLabels(self.initalParams)
        return [x[1] for x in self.initalParams]
    def getParams(self):
        return self.removeLabels(self.paramHistory)
    def removeLabels(self,param):
        return [[x[1] for x in paramStep] for paramStep in param]
    def getParamsWithLabels(self):
        return self.paramHistory
    def addNewParam(self,newParam):
        # new params come in the form of just numbers so we have to add the labels
        self.paramHistory.append(self.initalParams)
        for i in range(len(newParam)):
            self.paramHistory[-1][i][1] = newParam[i]


         



class SPOSimulationRunner:
    def __init__ (self,SPO):
        self.configuration = SPO.configuration
        # self.command = SPO.simulationCommand
        # (self.caller,self.sim,self.cmdLineArgs) = self.splitCommand()
        # self.parameters = SPO.parameters
        # self.runsOn = SPO.runsOn
        # self.name = SPO.name
        self.parameters = SPO.parameters
        self.step = SPO.step
        self.configFile = SPO.configFile
        self.maxJobs = SPO.maxJobs
        self.scriptName = "scriptRunner.sh"

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
        paramLog.write(str(self.parameters.getParamsWithLabels()[-1]))
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
        
        delimiter = "([^a-zA-Z0-9])"
        commandWithParams = self.configuration["Simulation:"]+' '

        # commandWithParams = self.cmdLineArgs
        # commandWithParams = " "+commandWithParams+" "
        for param in self.configuration["Parameters:"]:
            commandWithParams = re.sub(delimiter+param[0]+delimiter,'\g<1>'+str(param[1])+'\g<2>',commandWithParams)
        return commandWithParams +"\n"
    
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
        self.writeRerun(file)
    def writeRerun(self,file):
        file.write("cd ../.. \n")
        file.write("python3 "+str(__file__).replace(" ",r'\ ') + " " + self.configFile+"\n")
        

    
    def writeHPCCScript(self,file):
        # for not using an ensemble
        file.write("#SBATCH --job-name="+self.configuration["Name:"]+"\n")
        
        if self.configuration["Account:"] is not None:
            file.write("#SBATCH --account="+self.configuration["Account:"]+"\n")
        if self.configuration["Partition:"] is not None:
            file.write("#SBATCH --partition="+self.configuration["Partition:"]+"\n")

        if not self.configuration["Runs On:"][1]>1:
            for line in self.configuration["Extra Commands:"]:
                file.write(line+'\n')

            self.writeDesktopScript(file)
        else:
            #use a job array to submit the ensemble, 

            #file.write("#SBATCH --chdir " +self.path+"/$SLURM_ARRAY_TASK_ID\n")
            file.write("#SBATCH --array=0-"+str(self.configuration["Runs On:"][1]-1)+"%"+str(self.maxJobs)+"\n")
            for line in self.configuration["Extra Commands:"]:
                file.write(line+'\n')
            file.write("cd " + self.path+"/$SLURM_ARRAY_TASK_ID\n")
            file.write(self.createCommand())
            self.writeHPCRerun(file)


    def writeHPCRerun(self,file):
            file.write("cd ../../..\n")
            # after the simulation finished call this script again with the 
            # ensemble ID
            file.write("python3 "+str(__file__).replace(" ",r'\ ') + " " + self.configFile+" $SLURM_ARRAY_TASK_ID\n")

    def runScript(self):
        runString = ""
        if self.configuration["Runs On:"][0].value == SPORunsOn.HPCC.value:
            runString = "sbatch "+self.path+"/"+ self.scriptName
        else:
            runString = self.path+"/"+ self.scriptName + " &"
        subprocess.run(runString,shell=True)


class SPOSimulationRunnerWithOverhead(SPOSimulationRunner):

    def runScript(self):
        runString = ""
        if self.configuration["Runs On:"][0].value == SPORunsOn.HPCC.value:
            runString = "sbatch "+self.path+"/"+ self.scriptName
        else:
            runString = self.path+"/"+ self.scriptName
        output=subprocess.run(runString,shell=True,capture_output=True)
        # if we're running on the hpc we need to wait until all jobs are done
        if self.configuration["Runs On:"][0].value == SPORunsOn.HPCC.value:
            # first parse the the job ID
            jobMatch =  re.match("Submitted batch job (%d+)",output)
            if jobMatch:
                jobID = jobMatch.group(1)
            else:
                raise Exception("Failed to submit job")
            # TODO: Run one of the jobs on this node
            
            # now enter a loop to check every 5 minutes if the job is done
            
            #check the status of the job
            counter = 0
            twoDays = 48*60*60
            fiveMinutes = 5*60
            
            while not finished:
                finished = self.checkForFinishedJob(jobID)
                counter += fiveMinutes
                if counter>twoDays:
                    raise Exception("Simulation took longer than 2 days")
                time.sleep(fiveMinutes)
        # now all simulations have finished. we can continue

    def writeRerun(self,file):
        pass
    def writeHPCRerun(self,file):
        pass






class SPOOptimizer:
    def __init__(self,method,maxSteps,currentStep,space = None):
        self.method = method
        if self.method == "Hyperopt":
            self.space = list(hp.uniform(pspace[0],pspace[1],pspace[2]) for pspace in space)
        self.step = 0
        self.maxSteps = maxSteps
        self.currentStep = currentStep

    def get_next_parameters(self,parameters,residual):
        paramHistory = parameters.getParams()

        try:
            if self.method == "Hyperopt":
                def objectivefn(newParam):
                     self.objectiveFunction(newParam,paramHistory,residual)

                print("trying Hyperopt")
                self.outPut = fmin(objectivefn, self.space)
            else:

                self.outPut = minimize(self.objectiveFunction,parameters.getInitalConditions(),args=(paramHistory,residual),
                    method = self.method, options={"maxiter":self.maxSteps,"maxfun":self.currentStep+1,"eps":0.05})
                if not hasattr(self,'newParam'):
                    self.print(self.outPut)
        except StopIteration:
            pass
        return self.newParam

    def objectiveFunction(self,newParam,parameters,residual):
        # if this is a new parameter
        if self.step>=len(residual):
            # save the values
            self.newParam = newParam
            # return 0 to stop the optimization
            raise StopIteration

        # if this is not a new parameter make sure our parameter trajectory 
        # is valid
        paramDiff = [abs(x-y)>1e-6 for x,y in zip (newParam, parameters[self.step]) ]
        if any(paramDiff):
            print("encounted error on step %d"%self.step)
            print("List of params is %s"%str(parameters))
            print("I tried to use %s"%str(newParam))

            print([x-y for x,y in zip (newParam, parameters[self.step]) ])
            raise Exception("ERROR: encountered Nondeterministic solver")
        # and return the residual
        retVal = residual[self.step]
        self.step+=1
        return retVal

#This has a whole differnt interface apperntly so we don't inhareit... should think
# carfully about the structure of this
class SPOOptimizerWithOverhead():
    def __init__(self,SPO,runner):
        self.method = SPO.configuration["Method:"]
        if self.method == "Hyperopt":
            space = SPO.parameters.paramSpace
            self.space = list(hp.uniform(pspace[0],pspace[1],pspace[2]) for pspace in space)

        self.SPO = SPO
        self.runner = runner
    def run(self):
        #run the optimizer
        myMethod = self.SPO.configuration["Method:"]
        myEps = self.SPO.configuration["EPS:"]
        myTolerance = self.SPO.configuration["Tolerance:"]
        myMaxSteps = self.SPO.configuration["Max Steps:"]
        if self.method == "Hyperopt":
            self.outPut = fmin(self.objectiveFunction, self.space,max_evals=self.SPO.configuration["Max Steps:"])
        else:
            self.outPut = minimize(self.objectiveFunction,self.SPO.parameters.getInitalConditions(),
                method = myMethod, options={"maxiter":myMaxSteps,"eps":myEps})#,"tol":myTolerance})
        return self.outPut


    def objectiveFunction(self,newParam):
        # when we run a step of the optimizer we need do 
        self.SPO.parameters.addNewParam(newParam)
        # make the folder
        self.runner.createFolders()

        # write the new parameters to the log

        self.runner.createLogs()

        self.SPO.writeStartRunLog()

        # make the script
        self.runner.createScript()

        # run the script
        self.runner.runScript()

        # calculate the residual
        res = self.SPO.compareData()

        # write the residual to the log
        self.SPO.writeResidueLog(res)
        self.SPO.step+=1
        self.runner.step+=1
        self.SPO.updatePath()
        return res




if __name__ == "__main__":
    configFile = sys.argv[1]
    ensembleNo = None
    if len(sys.argv)>2:
        ensembleNo = int(sys.argv[2])

    mySPO = SimulationParameterOptimizer(sys.argv[1],ensembleNo)
    mySPO.run_with_overhead()
    # mySPO.run()
