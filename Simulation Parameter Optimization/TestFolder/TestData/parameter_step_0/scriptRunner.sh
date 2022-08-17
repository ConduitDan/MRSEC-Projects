#!/usr/bin/bash
cd SimpleSimTest/parameter_step_0
python3 ~/github/MRSEC\ Projects/Simulation\ Parameter\ Optimization/TestFolder/testSim.py 1.0 2.2 -2.0 
cd ../.. 
python3 /home/danny/github/MRSEC\ Projects/Simulation\ Parameter\ Optimization/TestFolder/../SPO.py TestConfigFile.txt
