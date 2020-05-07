#!/usr/bin/env python

# Written by Adrian Campos, QIMR Genetic Epidemiology. ("AC")



##############################################################################################################################################################
#                  Use the next code to test interactively (in an ipython console, make sure you run everything but the Arg parsing section)                 #
##############################################################################################################################################################
# remember to module load python and start ipython (use the ipython command cpaste to copy and paste if the identation fails)
# Copy everything inside the ''' adjusting your parameters !!
'''
class z():
	def __init__(self):
		self.input='inputList' #gwa sumstat list
		self.jobname='CFGWAStest' #job name
		self.minRsq=0.6
		self.minMAF=0.01
		self.nosub=True
		self.skipstep1=False
		self.skipstep2=False
		self.skipstep3=False
		self.cores=12
		self.estimation='DWLS'
args=z()
'''		

#####################################################
#                  Required modules                 #
#####################################################

from sys import exit
try:
	import pandas as pd
	import re
	import argparse
	import subprocess
	import os
	import os.path
	import sys
	import glob
	import time
	import sys
	import datetime
	import socket
	import numpy as np
	import warnings
	import string
	from os import path
except ImportError as error:
	print('Python module ('+error.name+') does not seem to be loaded')
	exit()
	
	
######################################################################
#                        Pree needed functions                       #
######################################################################
######################################################################

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

def SubmitScript(scriptName):
    ''' Submits all files with an extension in the current working directory'''
    import subprocess
    import re
    Submits=subprocess.Popen('qsub %s'%(scriptName),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    outlines=Submits.communicate()[0].decode("utf-8")
    jobs={i.split('.')[0]:0 for i in re.split(pattern='\n',string=outlines) if not i==''}
    return(jobs)
# Check completion function for submitted jobs to the hpc cluster
# Input a job dictionary of job IDs and zeros (means they are running), and a predetermined wait time between cluster query
def CheckCompletion(jobDic,timew=60):
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
			reresult=re.search('qstat: (\d+(|\[\]))\.hpc.* Job has finished',i)
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
	import subprocess
	import re
	Submits=subprocess.Popen("grep 'Error' %s"%(ext),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	outlines=Submits.communicate()[0].decode("utf-8")
	Errors=[i for i in re.split(pattern='\n',string=outlines) if not i=='']
	return(Errors)

######################################################################
#                        Start of main program                       #
######################################################################

print("By Adrian Campos @ QIMR Berghofer\nplease let me know of any problems errors\nadrian.campos@qimrberghofer.edu.au\n Thanks to Jackson Thorp for example and debugging advice")
print("This script is provided as is with no warranty and under under a CC license it is a wrapper for genomicSEM to perform a common factor GWAS. Whether this analysis makes sense is something to consider BEFORE running it")

# Store working directory

normalwd=os.getcwd() 

######################################################################
#                        Arg parsing                                 #
######################################################################


#SKIP THIS WHEN TESTING INTERACTIVELY! (RUN THE SECTIO FOR INTERACTIVE TESTING INSTEAD). !!!!!!!!!!!!!!!

parser = argparse.ArgumentParser(description="INPUT: summary statistic list, must have this columns (header should be included): file(path to the file) binary(is it case control 0/1) selogit(are s.e. in logistic scale 0/1) linprob(is it a case control ran using a linear model 0/1) sampprev(NA for continuous traits) popprev(NA for continuous traits). All sumstats must contain: SNP CHR BP A1 A2 FREQ INFO BETA SE P N different order is accepted. ",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input',help="Input file")
parser.add_argument('jobname',default='PRS',help="Job Name")
parser.add_argument('-minRsq',action='store',default=0.6,type=float,help='Filter out markers with imputed Rsq below this threshold')
parser.add_argument('-minMAF',action='store',default=0.01,type=float,help='Filter out markers with MAF below this threshold')
parser.add_argument('-cores',action='store',default=1,type=int,help='Number of cores to perform the common factor GWAS (step 3)')
parser.add_argument('-estimation',action='store',default='DWLS',type=str,help='Estimation method for common factor GWAS (step 3) either ML or DWLS')
parser.add_argument('-nosub',action='store_true',help="Do not submit any jobs (only generate PBS scripts)")
parser.add_argument('-skipstep1',action='store_true',help="Do not create nor submit step1 script (use this flag if it ran correctly)")
parser.add_argument('-skipstep2',action='store_true',help="Do not create nor submit step2 script (use this flag if it ran correctly)")
parser.add_argument('-skipstep3',action='store_true',help="Do not create nor submit step3 script (Not sure why you would want this)")


parser._actions[0].help="Print this help message and exit. For information on the required columns for the input file please go to: https://github.com/MichelNivard/GenomicSEM/wiki/4.-Common-Factor-GWAS"
args = parser.parse_args()
######################################################################
#                        Initial checks                              #
######################################################################
AnyError=0
if args.estimation not in ['DWLS','ML']:
	eprint("Estimation argument not correct, please double check it is either DLWS or ML")
	AnyError=1

if args.cores > 16:
	eprint("Warning you requested too many cores for the common factot GWAS, are you sure thats necessary?")

try:
	inputDataFrame=pd.read_csv(args.input,sep="\s+",usecols=['name','file','binary','selogit','linprob','sampprev','popprev'])
except ValueError:
	eprint("Could not read list of summary statistics, are you sure all the columns are present?\n Expected columns are:name	file(path to the file) binary(is it case control 0/1) selogit(are s.e. in logistic scale 0/1) linprob(is it a case control ran using a linear model 0/1) sampprev(NA for continuous traits) popprev(NA for continuous traits) ")
	AnyError=1
except FileNotFoundError:
	eprint("Could not find the list of summary statistics, are you using the full path? Are you out of your working directory?")
	AnyError=1
	

print("Interpreting input file. Please read the following carefully:")
selogit=[]
OLS=[]
linprob=[]
prop=[]
for row in inputDataFrame.index:
	if inputDataFrame.loc[row,'binary']:
		OLS.append('F')
		if inputDataFrame.loc[row,'linprob']:
			linprob.append('T')
			selogit.append('F')
			prop.append(inputDataFrame.loc[row,'sampprev'])
			print("%s is a binary phenotype, but analyzed with a linear model, thus standard errors cannot be in the logistic scale. GenomicSEM will transform this betas to the logistic scale based on the sample prevalence\n"%(inputDataFrame.loc[row,'file']))
		else:
			linprob.append('F')
			if inputDataFrame.loc[row,'selogit']:
				selogit.append('T')
			else:
				selogit.append('F')
			prop.append('NA')
			print("%s is a binary phenotype, analyzed with a logistic model. Your input for selogit (%s) will be used to determine whether standard errors are in the logistic scale\n"%(inputDataFrame.loc[row,'file'],inputDataFrame.loc[row,'selogit']))

	else:
		OLS.append('T')
		selogit.append('F')
		linprob.append('F')
		prop.append('NA')
		print("%s is a continuous phenotype thus a linear model was used, standard errors cannot be in the logistic scale and there is no sample or population prevalence\n"%(inputDataFrame.loc[row,'file']))
		
if AnyError:
	eprint("Please double check the errors above and resubmit the job")
	exit()
######################################################################
#                     Step one munging and LDSC                      #
######################################################################

if args.skipstep1:
	print("Skipping step 1 as requested by user")
else:
	print("Generating script one")
	script1="require(GenomicSEM);setwd('%s');files<-c("%(normalwd)
	script1+=','.join(["'%s'"%(i) for i in inputDataFrame.loc[:,'file']])
	script1+=");trait.names<-c("
	script1+=','.join(["'%s'"%(i) for i in inputDataFrame.loc[:,'name']])
	script1+=");"
	script1+="munge(files = files,hm3 = '/working/lab_nickm/adrianC/CommonData/w_hm3.noMHC.snplist',trait.names = trait.names, info.filter = 0.9, maf.filter = 0.01);"
	script1+="traits<-paste(trait.names,'.sumstats.gz',sep='');"
	script1+="sample.prev<-c("
	script1+=",".join(["%f"%(i) for i in inputDataFrame.loc[:,'sampprev']])
	script1+=");"
	script1+="population.prev<-c("
	script1+=",".join(["%f"%(i) for i in inputDataFrame.loc[:,'popprev']])
	script1+=");"
	script1+="ld <- '/working/lab_nickm/adrianC/CommonData/eur_w_ld_chr/';"
	script1+="wld <- '/working/lab_nickm/adrianC/CommonData/eur_w_ld_chr/';"
	script1+="LDSCoutput_stage2 <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names);"
	script1+="save(LDSCoutput_stage2,file = 'LDSCoutput_stage2.RData')"

	print("R script one looked like this (in a one liner):\n\n%s\n\nPlease check that it looks alright"%(re.sub(";","\n",script1)))

	JobScriptHeader=r'#!/bin/bash'
	scriptname="1_GSEM_multivariateLDSC.PBS" 
	#jobarray script for submission
	eprint("Creating Step 1 job script (%s).\n"%(scriptname))
	try:
		currscript=open(scriptname, 'w')
		currscript.write("%s\n"%(JobScriptHeader))
		currscript.write("#PBS -N CFGWAS_1\n")
		currscript.write("#PBS -l mem=32GB,walltime=10:00:00\n")
		currscript.write("cd $PBS_O_WORKDIR\n")
		currscript.write("module load R/3.5.1\n")
		currscript.write("Rscript -e \"%s\"\n"%(script1))
	except IOError:
		eprint("Error : could not write file %s !\n"%(scriptname))
	finally:
		currscript.close()
		
	if args.nosub:
		print("Finished writing step 1. Flag 'nosub' detected, will not submit to the cluster")
	else:
		print("Finished writing step 1. Will now submit it to the cluster")
		finishedJobs={}
		print("Submitted step 1 at %s"%(datetime.datetime.now()))
		Step1ID=SubmitScript(scriptname) # hash table with completion info

######################################################################
#                     Step two merging sumstats                      #
######################################################################
if args.skipstep2:
	print("Skipping step 2 as requested by user")
else:
	print("Generating script two")
	script2="require(GenomicSEM);setwd('%s');files<-c("%(normalwd)
	script2+=','.join(["'%s'"%(i) for i in inputDataFrame.loc[:,'file']])
	script2+=");trait.names<-c("
	script2+=','.join(["'%s'"%(i) for i in inputDataFrame.loc[:,'name']])
	script2+=");"
	script2+="ref<-'/working/lab_nickm/adrianC/CommonData/reference.1000G.maf.0.005.txt';"
	script2+="OLS<-c("
	script2+=",".join(OLS)
	script2+=");"
	script2+="se.logit<-c("
	script2+=",".join(selogit)
	script2+=");"
	script2+="linprob<-c("
	script2+=",".join(linprob)
	script2+=");"
	script2+="prop<-c("
	script2+=",".join(["%f"%(i) for i in prop])
	script2+=");"
	script2+="maf.filter<-%s;"%(args.minMAF)
	script2+="info.filter<-%s;"%(args.minRsq)
	script2+="combined_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS=OLS,linprob=linprob,prop=prop,info.filter=info.filter,maf.filter=maf.filter,parallel=FALSE);"
	script2+="save(combined_sumstats, file='%s/%s_combinedsumstats.RData');"%(normalwd,args.jobname)
	script2+="write.table(combined_sumstats, file = '%s/%s_combinedsumstats.tsv', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\\t')"%(normalwd,args.jobname)

	print("R script two looked like this (in a one liner):\n\n%s\n\nPlease check that it looks alright"%(re.sub(";","\n",script2)))


	JobScriptHeader=r'#!/bin/bash'
	scriptname="2_GSEM_sumstats.PBS" 
	#jobarray script for submission
	eprint("Creating Step 2 job script (%s).\n"%(scriptname))
	try:
		currscript=open(scriptname, 'w')
		currscript.write("%s\n"%(JobScriptHeader))
		currscript.write("#PBS -N CFGWAS_2\n")
		currscript.write("#PBS -l mem=32GB,walltime=10:00:00\n")
		currscript.write("cd $PBS_O_WORKDIR\n")
		currscript.write("module load R/3.5.1\n")
		currscript.write("Rscript -e \"%s\"\n"%(script2))
	except IOError:
		eprint("Error : could not write file %s !\n"%(scriptname))
	finally:
		currscript.close()

	if args.nosub:
		print("Finished writing step 2. Flag 'nosub' detected, will not submit to the cluster")
	else:
		print("Finished writing step 2. Will now submit it to the cluster")
		finishedJobs={}
		print("Submitted step 2 at %s"%(datetime.datetime.now()))
		Step2ID=SubmitScript(scriptname) # hash table with completion info
		finishedJobs={**Step1ID,**Step2ID}
		CheckCompletion(finishedJobs)
		# Checking for errors using simple grep
		anyErrors=CheckErrors()
		if len(anyErrors):
			print("The following errors were reported:\n")
			for currError in anyErrors:
				print(currError)
				print(currError)
				AnyError=1
		else:
			print("No obvious errors reported")

	if AnyError:
		eprint("Please double check the errors above and resubmit the job")
		exit()
######################################################################
#                 Step three common factor GWAS                      #
######################################################################
if args.skipstep1:
	print("Skipping step 1 as requested by user")
else:
	print("Generating script three")
	script3="require(GenomicSEM);require(hms);require(data.table);setwd('%s');"%(normalwd)
	script3+="Model_D_SS <- fread('%s/%s_combinedsumstats.tsv',header=TRUE);"%(normalwd,args.jobname)
	script3+="load('LDSCoutput_stage2.RData');"
	script3+="print('Starting the CommonFactorGWAS!');"
	script3+="results<-commonfactorGWAS(SNPs=Model_D_SS,covstruc=LDSCoutput_stage2,cores=%s,estimation='%s');"%(args.cores,args.estimation)
	script3+="print('Finished the CommonFactorGWAS!');"
	script3+="write.table(results, file =%s/%s_CFGWAS.dat, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\\t');"%(normalwd,args.jobname)

	print("R script three looked like this (in a one liner):\n\n%s\n\nPlease check that it looks alright"%(re.sub(";","\n",script3)))


	JobScriptHeader=r'#!/bin/bash'
	scriptname="3_GSEM_CommonFactorGWAS.PBS" 
	#jobarray script for submission
	eprint("Creating Step 3 job script (%s).\n"%(scriptname))
	try:
		currscript=open(scriptname, 'w')
		currscript.write("%s\n"%(JobScriptHeader))
		currscript.write("#PBS -N CFGWAS_3\n")
		currscript.write("#PBS -l ncpus=%s,mem=120GB,walltime=10:00:00\n"%(args.cores))
		currscript.write("cd $PBS_O_WORKDIR\n")
		currscript.write("module load R/3.5.1\n")
		currscript.write("Rscript -e \"%s\"\n"%(script3))
	except IOError:
		eprint("Error : could not write file %s !\n"%(scriptname))
	finally:
		currscript.close()

	if args.nosub:
		print("Finished writing step 3. Flag 'nosub' detected, will not submit to the cluster")
	else:
		print("Finished writing step 3. Will now submit it to the cluster")
		finishedJobs={}
		print("Submitted step 3 at %s"%(datetime.datetime.now()))
		finishedJobs=SubmitScript(scriptname) # hash table with completion info
		CheckCompletion(finishedJobs)
		# Checking for errors using simple grep
		anyErrors=CheckErrors()
		if len(anyErrors):
			print("The following errors were reported:\n")
			for currError in anyErrors:
				print(currError)
				print(currError)
				exit()
		else:
			print("No obvious errors reported")

print("Script finished running at %s"%(datetime.datetime.now()))
