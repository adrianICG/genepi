import glob
import subprocess
import re
import time
import datetime
import sys

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

def SubmitJobArray(scriptName,valuesFile,logdir=None):
	if logdir is None:
		logdir='./'
	'''Function to submit a job array, requires the PBS scripname, the values file name and an optional logdir. Returns jobID to keep track of it
	Resource allocation should be performed in the header of the PBS script'''
	Submit=subprocess.Popen("qsub -J $(echo \"1-$(wc -l %s|cut -f1 -d ' ')\") -o %s -e %s %s"%(valuesFile,logdir,logdir,scriptName),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	outraw=Submit.communicate()
	outlines=outraw[0].decode("utf-8") + outraw[1].decode("utf-8")
	if outlines.find('qsub') >=0:
		eprint('qsub returned this message :\n%s\n'%(outlines))
		AnyErrors=True
	elif outlines=='':
		eprint('qsub returned nothing.\n')
		AnyErrors=True
	else :
		finishedJobs={}
		jobID=outlines.rstrip()
		print('.. job name is %s'%(jobID))
		finishedJobs[jobID]=0
		return(finishedJobs)

def SubmitScript(scriptName):
	''' Submits a single PBS script. Returns JobID to keep track of it'''
	import subprocess
	import re
	Submits=subprocess.Popen('qsub %s'%(scriptName),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	outlines=Submits.communicate()[0].decode("utf-8")
	if outlines.find('qsub') >=0:
		eprint('qsub returned this message :\n%s\n'%(outlines))
		AnyErrors=True
	elif outlines=='':
		eprint('qsub returned nothing.\n')
		AnyErrors=True
	else:
		finishedJobs={}
		jobID=outlines.rstrip()
		print('.. job name is %s'%(jobID))
		finishedJobs[jobID]=0
		return(finishedJobs)


def isRunning(jobID):
	submitFinished=subprocess.Popen("qstat %s"%(jobID),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	qstatSTR=''.join([i.decode("utf-8") for i in submitFinished.communicate()])
	reresult=re.search('qstat: (\d+(|\[\]))\.hpc.* Job has finished',qstatSTR)
	if reresult is None:
		return(True)
	else:
		return(False)

def EnforceJobResources(jobID,mem=True,cpus=True,timew=60):
	'''Check completiom function but also enforces a job or job array sticks to their allocated resources
	IMPORTANT! This function will KILL the job (or subjobs for an array) even if
	the memory or cpu breach is not great. CPU restrictions are a bit relaxed 
	e.g. if you request 12 cores we allow a max cpu usage of 1300.
	This function naturally keeps running as long as your job is running.'''
	if not isRunning(jobID):
		print("The job is not currently running")
		return(None)
	while(isRunning(jobID)):
		if re.search("\[",jobID):#Job array
			submitSubJobIDs=subprocess.Popen("qstat -1nt %s | cut -f 1 -d '.'"%(jobID),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			out=''.join([i.decode("utf-8") for i in submitSubJobIDs.communicate()])
			subJobIDs=[i for i in re.split('\n',out) if re.match("\d+\[\d+\]", i)]
			breaching=0
			for currjobID in subJobIDs:
				if isRunning(currjobID):
					submits=subprocess.Popen("qstat -f %s | grep 'resources_used.mem\|Resource_List.mem ='"%(currjobID),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
					submitsCPU=subprocess.Popen("qstat -f %s | grep 'resources_used.cpupercent\|Resource_List.ncpus'"%(currjobID),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
					out=''.join([i.decode("utf-8") for i in submits.communicate()])
					outCPU=''.join([i.decode("utf-8") for i in submitsCPU.communicate()])
					cpuusedline,cpurequestedline=[i for i in re.split(pattern='\n',string=outCPU) if not i=='']
					usedline,requestedline=[i for i in re.split(pattern='\n',string=out) if not i=='']
					usedmem=int(re.search('(\d+)',usedline).group(1))/1e+6
					requestedmem=int(re.search('(\d+)',requestedline).group(1))
					usedcpu=int(int(re.search('(\d+)',cpuusedline).group(1))/100)
					requestedcpu=int(re.search('(\d+)',cpurequestedline).group(1))
					if requestedmem<usedmem:
						print("Job %s using %s but requested %s will terminate the job"%(currjobID,usedmem,requestedmem))
						submits=subprocess.Popen("qdel %s"%(currjobID),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
						breaching+=1
					elif requestedcpu<usedcpu:
						print("Job %s using %s cores but requested %s will terminate the job"%(currjobID,usedcpu,requestedcpu))
						submits=subprocess.Popen("qdel %s"%(currjobID),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
						breaching+=1
			if breaching == 0:
				print("All subjobs are within memory and cpu specifications")
			else:
				print("At this round %s subjobs were killed because they use more memory than requested"%(breaching))
		else: #Single job
			submits=subprocess.Popen("qstat -f %s | grep 'resources_used.mem\|Resource_List.mem ='"%(jobID),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			submitsCPU=subprocess.Popen("qstat -f %s | grep 'resources_used.cpupercent\|Resource_List.ncpus'"%(jobID),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			out=''.join([i.decode("utf-8") for i in submits.communicate()])
			outCPU=''.join([i.decode("utf-8") for i in submitsCPU.communicate()])
			cpuusedline,cpurequestedline=[i for i in re.split(pattern='\n',string=out) if not i=='']
			usedline,requestedline=[i for i in re.split(pattern='\n',string=out) if not i=='']
			usedmem=int(re.search('(\d+)',usedline).group(1))/1e+6
			requestedmem=int(re.search('(\d+)',requestedline).group(1))
			usedcpu=int(int(re.search('(\d+)',cpuusedline).group(1))/100)
			requestedcpu=int(re.search('(\d+)',cpurequestedline).group(1))
			if requestedmem<usedmem:
				print("Job %s using %s but requested %s will terminate the job"%(jobID,usedmem,requestedmem))
				submits=subprocess.Popen("qdel %s"%(jobID),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			elif requestedcpu<usedcpu:
				print("Job %s using %s but requested %s will terminate the job"%(jobID,usedcpu,requestedcpu))
				submits=subprocess.Popen("qdel %s"%(jobID),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			else:
				print("Job currently sticking to its requirements")
		time.sleep(timew)	

# Check completion function for submitted jobs to the hpc cluster
# Input a job dictionary of job IDs and zeros (means they are running), and a predetermined wait time between cluster query
def CheckCompletion(jobDic,timew=60):
	'''Keep track of running jobs. Useful to wait for submitted jobs to finish before moving on. 
	Receives a dictionary of jobIDs (output of submission scripts above) to keep track of
	You can concatenate multiple dictionaries doing: {**dic1,**dic2,**dic3...**dicN}'''
#   print("into CheckCompletion\n")
	ids=[i for i in jobDic.keys()] #get all job ids 
	print("At %s : Watching job ids=\n%s\n"%(datetime.datetime.now(),ids))
	submstr=str("qstat %s 2>&1"%(' '.join(ids))) #submit query for all jobs
#   print("submstr=\n%s\n"%(submstr))  
	firstiter=1
	njobs=len(jobDic)
	ncomplete=0
	ncompleteprev=0
	while ncomplete < njobs:
		if firstiter != 1:
			time.sleep(timew)
		# Query status of the requested jobs.
		submits=subprocess.Popen(submstr,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		# Get query output.
		out=''.join([i.decode("utf-8") for i in submits.communicate()])
#       print('out=\n%s'%(out))
		# Split lines of output.
		OutLines=[i for i in re.split(pattern='\n',string=out) if not i=='']
		# Get ids of finished jobs. Count. Mark as complete.
		ncomplete=0
		for i in OutLines:
			reresult=re.search('qstat: (\d+(|\[\])\.hpc.*) Job has finished',i)
			if reresult != None:
				id=reresult.group(1)
				if jobDic[id] != 1:
					print("Job %s has finished.\n"%(id))
				jobDic[id] = 1
				ncomplete += 1
		# Report progress if anything has changed.
#       print("njobs=%s ncomplete=%s ncompleteprev=%s\n"%(njobs,ncomplete,ncompleteprev))
		if (firstiter==1) or (ncomplete != ncompleteprev):
			print("At %s : currently %s of %s job(s) have finished.\n"%(datetime.datetime.now(),ncomplete,njobs))
		firstiter=0
		ncompleteprev=ncomplete
	return(1) # Return something just in case you need to know.

######################################################################
# Error checker, basically greps Error in all of the files *.e* of the currwd
def CheckErrors(ext='*.e*'):
	'''Very naive way to check for errors from the script. Basically greps out the error output files from the hpc
	This is more suited for after submission of a single job'''
	import subprocess
	import re
	Submits=subprocess.Popen("grep 'Error' %s"%(ext),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	outlines=Submits.communicate()[0].decode("utf-8")
	Errors=[i for i in re.split(pattern='\n',string=outlines) if not i=='']
	return(Errors)
	
def CheckErrorsArray(logdir):
	'''Will try detect errors from the output of an array job (default ending in .ER[id]
	Input is the logfile directory'''
	import subprocess
	import re
	AnyErrors=False
	logs=glob.glob('%s/*.ER'%(logdir))
	Errors=[]
	for i in logs:
		Submits=subprocess.Popen("grep 'Error' %s"%(i),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		errlines=Submits.communicate()[0].decode("utf-8")
		Errors=[i for i in re.split(pattern='\n',string=errlines) if not i=='']
	return(Errors)