import os
import subprocess
import time
from shutil import copyfile
import sys
import json,collections
import os, glob

#remove the energyplus files from the previous run
for filename in glob.glob("eplusout*"):
    os.remove(filename) 

f=open('airflow.csv','w')
f.close()
f=open('temp.csv','w')
f.close()
f=open('power.csv','w')
f.close()	

configdir='config'


modelPath=''

os.environ['BCVTB_HOME']='bcvtb'

weatherPath='USA_WA_Pasco-Tri.Cities.AP.727845_TMY3.epw'

modelDir='BUILDING11.idf'

cmdStr = "energyplus -w \"%s\" -r \"%s\"" % (weatherPath, modelDir)

simulation = subprocess.Popen(cmdStr, shell=True)
sock=subprocess.Popen('python '+'master_nn-test.py'+' '+str(modelPath)+' '+str(configdir), shell=True)


sock.wait()
#sock.terminate()
