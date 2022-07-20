from asyncio import subprocess
import unittest
import os
import SPO
import subprocess
import types
import shutil
import numpy as np


class ParserSpec(unittest.TestCase):
    def test_read_line(self):
        parser = SPO.SPOFileParser("TestFolder/ParserTest/TestParserBasic.txt")
        self.assertEqual(list(parser),["Test Line"])

    def test_whitespace_skip(self):
        parser = SPO.SPOFileParser("TestFolder/ParserTest/TestParserWhiteSpace.txt")
        self.assertEqual(list(parser),["test Line","next Line"])

    def test_comments(self):
        parser = SPO.SPOFileParser("TestFolder/ParserTest/TestParserComment.txt")
        self.assertEqual(list(parser),["content"])

    def test_parseConfigFile(self):
        parser = SPO.SPOFileParser("TestFolder/TestConfigFile.txt")
        (name, simulationCommand, parameters, dataSpec, runsOn, method) = parser.parseConfigFile()
        self.assertEqual(name,"SimpleSimTest")
        self.assertEqual(simulationCommand,"python3 testSim.py x y z")
        self.assertEqual(parameters,[["x",1],["y",2.2],["z",-2]])
        self.assertEqual(dataSpec,[["testSimData.txt", "testTargetData.txt"]])
        self.assertEqual(runsOn,[SPO.SPORunsOn.Desktop])
        self.assertEqual(method,"L-BFGS-B")

    def test_read_log(self):
        parser = SPO.SPOFileParser("TestFolder/ParserTest/TestingLog.txt")
        (params,residues,step) = parser.parseLogFile()
        self.assertEqual(params,[[["x",1],["y",2.2],["z",-2]],[["x",2],["y",3],["z",4]]])
        self.assertEqual(residues,[47.630624000000005])
        self.assertEqual(step,1)
        
        
    def test_read_ensemble_log_ready(self):
        f = open("TestFolder/ParserTest/ensembleLogReady.txt",'w')
        f.write("Run 0 of 2 Running \n")
        f.write("Run 1 of 2 Running ")
        f.close()
        myParser = SPO.SPOEnsembleLogParser("TestFolder/ParserTest/ensembleLogReady.txt")
        self.assertFalse(myParser.parseEnsembleLog(1))
        self.assertTrue(myParser.parseEnsembleLog(0))
        

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
        mySPO.runsOn = [SPO.SPORunsOn.HPCC,2]
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
        self.simpleSPO.simulationCommand = "matlab testcommand x1 x2 "
        self.simpleSPO.parameters = [["x1",1],["x2",2]]
        self.simpleSPO.runsOn = [SPO.SPORunsOn.Desktop]
        self.simpleSPO.name = "Simulation_runner_tests"
        self.simpleSPO.step = 0
        self.simpleSPO.configFile = "TestFolder/TestConfigFile"
        self.simpleSPO.partition = "guest"
        self.simpleSPO.maxJobs = 10
        os.mkdir(self.simpleSPO.name)

    def tearDown(self):
        shutil.rmtree(self.simpleSPO.name)
        
    def test_desktop_run_folder(self):
        runner = SPO.SPOSimulationRunner(self.simpleSPO)
        runner.createFolders()
        self.assertTrue(os.path.isdir("Simulation_runner_tests/parameter_step_0"))

        #check the run script
        runner.createScript()
        self.assertTrue(os.path.exists("Simulation_runner_tests/parameter_step_0/scriptRunner.sh"))
        script = open("Simulation_runner_tests/parameter_step_0/scriptRunner.sh")

        self.assertEqual(script.readline(),"#!/usr/bin/bash\n")
        self.assertEqual(script.readline(),"cd %s/parameter_step_0\n"%self.simpleSPO.name)
        file = os.getcwd() + "/testcommand"

        self.assertEqual(script.readline(),"matlab " + file.replace(' ','\ ') + "  1 2  \n")
        self.assertEqual(script.readline(),"cd ../.. \n")
        self.assertEqual(script.readline(),"python3 " + os.getcwd().replace(' ','\ ')+ "/SPO.py "+self.simpleSPO.configFile +'\n')
        script.close()

    def test_hpcc_run_folder(self):
        self.simpleSPO.runsOn = [SPO.SPORunsOn.HPCC,1]
        runner = SPO.SPOSimulationRunner(self.simpleSPO)
        runner.createFolders()
        self.assertTrue(os.path.isdir("Simulation_runner_tests/parameter_step_0"))
        runner.createScript()

                #check the run script
        runner.createScript()
        self.assertTrue(os.path.exists("Simulation_runner_tests/parameter_step_0/scriptRunner.sh"))
        script = open("Simulation_runner_tests/parameter_step_0/scriptRunner.sh")

        self.assertEqual(script.readline(),"#!/usr/bin/bash\n")
        self.assertEqual(script.readline(),"#SBatch --job-name=Simulation_runner_tests\n")
        self.assertEqual(script.readline(),"#SBatch --partition=guest\n")

        self.assertEqual(script.readline(),"cd %s/parameter_step_0\n"%self.simpleSPO.name)
        file = os.getcwd() + "/testcommand"

        self.assertEqual(script.readline(),"matlab " + file.replace(' ','\ ') + "  1 2  \n")
        self.assertEqual(script.readline(),"cd ../.. \n")
        self.assertEqual(script.readline(),"python3 " + os.getcwd().replace(' ','\ ')+ "/SPO.py "+self.simpleSPO.configFile+'\n')
        script.close()

    def test_hpcc_ensemble_run_folder(self):
        self.simpleSPO.runsOn = [SPO.SPORunsOn.HPCC,2]
        runner = SPO.SPOSimulationRunner(self.simpleSPO)
        runner.createFolders()
        self.assertTrue(os.path.isdir("Simulation_runner_tests/parameter_step_0/0"))
        self.assertTrue(os.path.isdir("Simulation_runner_tests/parameter_step_0/1"))
        runner.createScript()
        self.assertTrue(os.path.exists("Simulation_runner_tests/parameter_step_0/scriptRunner.sh"))
        script = open("Simulation_runner_tests/parameter_step_0/scriptRunner.sh")

        self.assertEqual(script.readline(),"#!/usr/bin/bash\n")
        self.assertEqual(script.readline(),"#SBatch --job-name=Simulation_runner_tests\n")
        self.assertEqual(script.readline(),"#SBatch --partition=guest\n")
        self.assertEqual(script.readline(),"#SBatch --chdir Simulation_runner_tests/parameter_step_0/$SLURM_ARRAY_TASK_ID\n")
        self.assertEqual(script.readline(),"#SBatch --array=0-1%10\n")
        
        file = os.getcwd() + "/testcommand"

        self.assertEqual(script.readline(),"matlab " + file.replace(' ','\ ') + "  1 2  \n")
        self.assertEqual(script.readline(),"cd ../..\n")
        self.assertEqual(script.readline(),"python3 " + os.getcwd().replace(' ','\ ')+ "/SPO.py "+self.simpleSPO.configFile + " $SLURM_ARRAY_TASK_ID\n")
        script.close()



class systemTest(unittest.TestCase):
    def test_simpleSimulation(self):
        os.system("(cd TestFolder; python3 ../SPO.py TestConfigFile.txt)")


if __name__ == '__main__':
    unittest.main()