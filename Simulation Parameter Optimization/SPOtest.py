from asyncio import subprocess
import unittest
import os
import SPO
import SPOParser
import subprocess
import types
import shutil
import numpy as np


class ParserSpec(unittest.TestCase):
    
    def test_read_line(self):
        parser = SPOParser.SPOFileParser("TestFolder/ParserTest/TestParserBasic.txt")
        self.assertEqual(list(parser),["Test Line"])

    def test_whitespace_skip(self):
        parser = SPOParser.SPOFileParser("TestFolder/ParserTest/TestParserWhiteSpace.txt")
        self.assertEqual(list(parser),["test Line","next Line"])

    def test_comments(self):
        parser = SPOParser.SPOFileParser("TestFolder/ParserTest/TestParserComment.txt")
        self.assertEqual(list(parser),["content"])

class FileSpec(unittest.TestCase):
    def test_parseConfigFile(self):
        configFile = SPOParser.FileSpecFactory().configFileSpec()
        vals = configFile.readFile("TestFolder/TestConfigFile.txt")
        self.assertEqual(vals["Name:"],"SimpleSimTest")
        self.assertEqual(vals["Simulation:"],"python3 ~/github/MRSEC\ Projects/Simulation\ Parameter\ Optimization/TestFolder/testSim.py x y z")
        self.assertEqual(vals["Parameters:"],[["x",1],["y",2.2],["z",-2]])
        self.assertEqual(vals["Data:"],[("testSimData.txt", "testTargetData.txt")])
        self.assertEqual(vals["Runs On:"],[SPO.SPORunsOn.Desktop,1])
        self.assertEqual(vals["Method:"],"L-BFGS-B")

    # def test_parseConfigFileMatlab(self):
    #     # Look I know this test way more than just the parser
    #     # The point is this will be refactored so all this logic is just in the 
    #     # parser and the File and Field Specs
    #     mySPO = SPO.SimulationParameterOptimizer("TestFolder/TestConfigFileMatlab.txt",None)
    #     mySPO.step = 0
    #     runner = SPO.SPOSimulationRunner(mySPO)
    #     self.assertEqual(runner.caller,"matlab -nojvm -batch ")


    def test_read_log(self):
        log = SPOParser.FileSpecFactory().logFileSpec()
        vals = log.readFile("TestFolder/ParserTest/TestingLog.txt")
        linesWithResidue = vals[0]
        linesWithOutResidue = vals[1]
        params = list(val[0] for val in linesWithResidue)
        params.append(linesWithOutResidue[0])

        residues = list(val[1] for val in linesWithResidue)

        step = linesWithOutResidue[1]
        self.assertEqual(params,[[["x",1],["y",2.2],["z",-2]],[["x",2],["y",3],["z",4]]])
        self.assertEqual(residues,[47.630624000000005])
        self.assertEqual(step,1)
        
        
    def test_read_ensemble_log_ready(self):
        f = open("TestFolder/ParserTest/ensembleLogReady.txt",'w')
        f.write("Run 0 of 2 Running \n")
        f.write("Run 1 of 2 Running ")
        f.close()
        myParser = SPOParser.SPOEnsembleLogParser("TestFolder/ParserTest/ensembleLogReady.txt")
        self.assertFalse(myParser.parseEnsembleLog(1))
        self.assertTrue(myParser.parseEnsembleLog(0))


    def test_read_ensemble_log_Lvl2(self):
        # Write a fresh Log
        f = open("TestFolder/TestData/parameter_step_0/ensembleLog.txt",'w')
        f.write("Run 0 of 2 Running \n")
        f.write("Run 1 of 2 Running ")
        f.close()
        mySPO = SPO.SimulationParameterOptimizer("TestFolder/TestConfigFileHPCC.txt",0)
        mySPO.logFileName = "TestFolder/TestData/SimpleSimTestLog.txt"
        mySPO.run() #all ths should do is replace that first line with "Run 0 of 2 Finished"

        f = open("TestFolder/TestData/parameter_step_0/ensembleLog.txt",'r')
        self.assertEqual(f.readline(),"Run 0 of 2 Finished \n")
        f.close()
        



        

class data_read_tests(unittest.TestCase):
    def test_read_datum(self):
        mySPO = SPO.SimulationParameterOptimizer("TestFolder/TestConfigFile.txt",None)
        datum = mySPO.loadData("TestFolder/TestData/TestDatum.txt")
        self.assertEqual(datum,1)


    def test_read_multiple_columns(self):
        mySPO = SPO.SimulationParameterOptimizer("TestFolder/TestConfigFile.txt",None)
        datum = mySPO.loadData("TestFolder/TestData/TestDataArray.txt")
        self.assertTrue((datum==np.array([[1,2],[3,4]])).all())



    def test_read_ensemble_data(self):
        mySPO = SPO.SimulationParameterOptimizer("TestFolder/TestConfigFile.txt",None)
        mySPO.configuration["Runs On:"] = [SPO.SPORunsOn.HPCC,2]
        mySPO.path = "TestFolder/TestData/"
        ensembleData = mySPO.getSimData("TestData.txt")
        self.assertTrue((ensembleData==np.array([2,3])).all())
        


# class optimizer_tests(unittest.TestCase):
#     def test_past_value_compare(self):
#         pass

#     def test_nondetermistic_optimizer_error(self):
#         pass

class simulation_runner_tests(unittest.TestCase):
    def setUp(self):
        self.simpleSPO = types.SimpleNamespace()
        config = {"Simulation:":"matlab testcommand x1 x2 ",
                    "Name:":"Simulation_runner_tests",
                    "Runs On:":[SPO.SPORunsOn.Desktop,1],
                    "Partition:":"guest",
                    "Account:":"guest",
                    "Extra Commands:":"",
                    "Parameters:":[["x1",1],["x2",2]]}

        self.simpleSPO.configuration = config

        self.simpleSPO.step = 0
        self.simpleSPO.configFile = "TestFolder/TestConfigFile"
        self.simpleSPO.maxJobs = 10
        os.mkdir(self.simpleSPO.configuration["Name:"])

    def tearDown(self):
        shutil.rmtree(self.simpleSPO.configuration["Name:"])
        
    def test_desktop_run_folder(self):
        runner = SPO.SPOSimulationRunner(self.simpleSPO)
        runner.createFolders()
        self.assertTrue(os.path.isdir("Simulation_runner_tests/parameter_step_0"))

        #check the run script
        runner.createScript()
        self.assertTrue(os.path.exists("Simulation_runner_tests/parameter_step_0/scriptRunner.sh"))
        script = open("Simulation_runner_tests/parameter_step_0/scriptRunner.sh")

        self.assertEqual(script.readline(),"#!/usr/bin/bash\n")
        self.assertEqual(script.readline(),"cd %s/parameter_step_0\n"%self.simpleSPO.configuration["Name:"])
        file = os.getcwd() + "/testcommand"

        self.assertEqual(script.readline(),"matlab testcommand 1 2  \n")
        self.assertEqual(script.readline(),"cd ../.. \n")
        self.assertEqual(script.readline(),"python3 " + os.getcwd().replace(' ','\ ')+ "/SPO.py "+self.simpleSPO.configFile +'\n')
        script.close()

    def test_hpcc_run_folder(self):
        self.simpleSPO.configuration["Runs On:"] =  [SPO.SPORunsOn.HPCC,1]
        runner = SPO.SPOSimulationRunner(self.simpleSPO)
        runner.createFolders()
        self.assertTrue(os.path.isdir("Simulation_runner_tests/parameter_step_0"))
        runner.createScript()

                #check the run script
        runner.createScript()
        self.assertTrue(os.path.exists("Simulation_runner_tests/parameter_step_0/scriptRunner.sh"))
        script = open("Simulation_runner_tests/parameter_step_0/scriptRunner.sh")

        self.assertEqual(script.readline(),"#!/usr/bin/bash\n")
        self.assertEqual(script.readline(),"#SBATCH --job-name=Simulation_runner_tests\n")
        self.assertEqual(script.readline(),"#SBATCH --account=guest\n")
        self.assertEqual(script.readline(),"#SBATCH --partition=guest\n")

        self.assertEqual(script.readline(),"cd %s/parameter_step_0\n"%self.simpleSPO.configuration["Name:"])
        #file = os.getcwd() + "/testcommand"

        self.assertEqual(script.readline(),"matlab testcommand 1 2  \n")
        self.assertEqual(script.readline(),"cd ../.. \n")
        self.assertEqual(script.readline(),"python3 " + os.getcwd().replace(' ','\ ')+ "/SPO.py "+self.simpleSPO.configFile+'\n')
        script.close()

    def test_hpcc_ensemble_run_folder(self):
        self.simpleSPO.configuration["Runs On:"] = [SPO.SPORunsOn.HPCC,2]
        runner = SPO.SPOSimulationRunner(self.simpleSPO)
        runner.createFolders()
        self.assertTrue(os.path.isdir("Simulation_runner_tests/parameter_step_0/0"))
        self.assertTrue(os.path.isdir("Simulation_runner_tests/parameter_step_0/1"))
        runner.createScript()
        self.assertTrue(os.path.exists("Simulation_runner_tests/parameter_step_0/scriptRunner.sh"))
        script = open("Simulation_runner_tests/parameter_step_0/scriptRunner.sh")

        self.assertEqual(script.readline(),"#!/usr/bin/bash\n")
        self.assertEqual(script.readline(),"#SBATCH --job-name=Simulation_runner_tests\n")
        self.assertEqual(script.readline(),"#SBATCH --account=guest\n")
        self.assertEqual(script.readline(),"#SBATCH --partition=guest\n")
        #self.assertEqual(script.readline(),"#SBATCH --chdir Simulation_runner_tests/parameter_step_0/$SLURM_ARRAY_TASK_ID\n")
        self.assertEqual(script.readline(),"#SBATCH --array=0-1%10\n")
        
        self.assertEqual(script.readline(),"cd Simulation_runner_tests/parameter_step_0/$SLURM_ARRAY_TASK_ID\n")
        self.assertEqual(script.readline(),"matlab testcommand 1 2  \n")
        self.assertEqual(script.readline(),"cd ../../..\n")
        self.assertEqual(script.readline(),"python3 " + os.getcwd().replace(' ','\ ')+ "/SPO.py "+self.simpleSPO.configFile+' $SLURM_ARRAY_TASK_ID\n')
        script.close()



class systemTest(unittest.TestCase):
    def test_simpleSimulation(self):
        shutil.rmtree("TestFolder/SimpleSimTest/")
        os.system("(cd TestFolder; python3 ../SPO.py TestConfigFile.txt)")
    def test_hyperopt(self):
        shutil.rmtree("TestFolder/SimpleSimTestHyperopt/")
        os.system("(cd TestFolder; python3 ../SPO.py TestConfigFile2.txt)")
        


if __name__ == '__main__':
    mySystemTest=systemTest()
    mySystemTest.test_simpleSimulation()
    # mySystemTest.test_hyperopt()
    # unittest.main()