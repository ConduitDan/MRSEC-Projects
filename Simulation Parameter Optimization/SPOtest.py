from asyncio import subprocess
import unittest
import SPO
import subprocess

class ParserSpec(unittest.TestCase):
    # def read_line:
    #     parser = SPO.SPOFileParser("TestLine.txt")
    #     parser.

    # def skip_newline:

    # def skip_comment:
    def test_parseConfigFile(self):
        parser = SPO.SPOFileParser("TestFolder/TestConfigFile.txt")
        (name, simulationCommand, parameters, dataSpec, runsOn, method) = parser.parseConfigFile()
        self.assertEqual(name,"Testing")
        self.assertEqual(simulationCommand,"python3 testSim.py x y z")
        self.assertEqual(parameters,[["x",1],["y",2.2],["z",-2]])
        self.assertEqual(dataSpec,[["testSimData.txt", "testTargetData.txt"]])
        self.assertEqual(runsOn,"Desktop")
        self.assertEqual(method,"HyperOpt")



class systemTest(unittest.TestCase):
    def test_simpleSimulation(self):
        print(subprocess.run("python3 SPO.py TestFolder/TestConfigFile.txt"),capture_output=True,shell=True)


if __name__ == '__main__':
    unittest.main()