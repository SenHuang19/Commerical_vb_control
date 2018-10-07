import os
import subprocess
import time
from shutil import copyfile
import sys
import json,collections

weatherPath='USA_WA_Pasco-Tri.Cities.AP.727845_TMY3.epw'

modelDir='baseline\BUILDING11.idf'

cmdStr = "C:\EnergyPlusV8-4-0\energyplus -w \"%s\" -r \"%s\"" % (weatherPath, modelDir)

simulation = subprocess.Popen(cmdStr, shell=True)
simulation .wait()
copyfile('eplusout.csv', 'baseline.csv')
