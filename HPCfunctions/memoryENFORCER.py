#!/usr/bin/env python
# memoryENFORCER script to keep problematic jobs from using more memory than requested
# Written by Adrian Campos, QIMR Genetic Epidemiology. ("AC")

import subprocess
import re
import time
import argparse

parser = argparse.ArgumentParser(description="HPC PBS Memory enforcer input is JobID ",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input',help="Input the jobID, this scripot will kill it if it uses more memory than requested")
parser._actions[0].help="Print this help message and exit. For information on the required"
args = parser.parse_args()

jobID=args.input()
running=1
timew=60
while(running):
	submitFinished=subprocess.Popen("qstat %s"%(jobID),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	qstatSTR=''.join([i.decode("utf-8") for i in submitFinished.communicate()])
	reresult=re.search('qstat: (\d+(|\[\]))\.hpc.* Job has finished',qstatSTR)
	if reresult is None:
		submits=subprocess.Popen("qstat -f %s | grep 'resources_used.mem\|Resource_List.mem = 120gb'"%(jobID),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		out=''.join([i.decode("utf-8") for i in submits.communicate()])
		usedline,requestedline=[i for i in re.split(pattern='\n',string=out) if not i=='']
		usedmem=int(re.search('(\d+)',usedline).group(1))/1e+6
		requestedmem=int(re.search('(\d+)',requestedline).group(1))
		if requestedmem<usedmem:
			print("Job using %s but requested %s will terminate the job"%(usedmem,requestedmem))
			submits=subprocess.Popen("qdel %s"%(jobID),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		else:
			print("Job using %s but requested %s the job keeps running"%(usedmem,requestedmem))
			time.sleep(timew)
	else:
		running=0

print("The job finished without exceding its memory requirement")
