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
		self.cores=1
		self.estimation='DWLS'
		self.splitJobs=True
		self.retry=True
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
	from os import path
	from QIMRPBSutils import SubmitJobArray,SubmitScript,CheckCompletion,CheckErrors,CheckErrorsArray,eprint
except ImportError as error:
	print('Python module does not seem to be loaded')
	exit()


######################################################################
#                        Start of main program                       #
######################################################################

print("By Adrian Campos @ QIMR Berghofer\nplease let me know of any problems errors\nadrian.campos@qimrberghofer.edu.au\n Thanks to Jackson Thorp for example and debugging advice")
print("This script is provided as is with no warranty and under under a CC license it is a wrapper for genomicSEM to perform a common factor GWAS. Whether this analysis makes sense is something to consider before running it")

# Store working directory

normalwd=os.getcwd() 

######################################################################
#                        Arg parsing                                 #
######################################################################


#SKIP THIS WHEN TESTING INTERACTIVELY! (RUN THE SECTIO FOR INTERACTIVE TESTING INSTEAD). !!!!!!!!!!!!!!!

parser = argparse.ArgumentParser(description="INPUT: summary statistic list, must have these columns (header should be included):name(trait name) file(path to the file) binary(is it case control 0/1) selogit(are s.e. in logistic scale 0/1) linprob(is it a case control ran using a linear model 0/1) sampprev(NA for continuous traits) popprev(NA for continuous traits). All sumstats must contain: SNP CHR BP A1 A2 FREQ INFO BETA SE P N different order is accepted. ",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input',help="Input file") 
parser.add_argument('jobname',default='PRS',help="Job Name")
parser.add_argument('-minRsq',action='store',default=0.6,type=float,help='Filter out markers with imputed Rsq below this threshold')
parser.add_argument('-minMAF',action='store',default=0.01,type=float,help='Filter out markers with MAF below this threshold')
parser.add_argument('-cores',action='store',default=1,type=int,help='Number of cores to perform the common factor GWAS (step 3). Raise this if you do not use -splitJobs')
parser.add_argument('-estimation',action='store',default='DWLS',type=str,help='Estimation method for common factor GWAS (step 3) either ML or DWLS')
parser.add_argument('-nosub',action='store_true',help="Do not submit any jobs (only generate PBS scripts)")
parser.add_argument('-skipstep1',action='store_true',help="Do not create nor submit step1 script (use this flag if it ran correctly)")
parser.add_argument('-skipstep2',action='store_true',help="Do not create nor submit step2 script (use this flag if it ran correctly)")
parser.add_argument('-skipstep3',action='store_true',help="Do not create nor submit step3 script (Not sure why you would want this)")
parser.add_argument('-splitJobs',action='store_true',help="Split the common factor GWAS into chunks of 100k SNPs this is faster but less stable")
parser.add_argument('-retry',action='store_true',help="When the common factor GWAS jobs are split, subjobs very often fail without reason, but finish when ran again. This flag will reattempt to run non-finishing jobs up to 10 times which should be more than enough to get them to finish")



parser._actions[0].help="Print this help message and exit. HIGHLY RECOMMENDED TO INCLUDE: -splitJobs and -retry flags to speed up computation. For information on the required columns for the input file please go to: https://github.com/MichelNivard/GenomicSEM/wiki/4.-Common-Factor-GWAS"
args = parser.parse_args()
######################################################################
#                        Initial checks                              #
######################################################################
AnyError=0
if args.estimation not in ['DWLS','ML']: #If estimation is not one of the correct ones
	eprint("Estimation argument not correct, please double check it is either DLWS or ML") # eprint is defined in QIMRPBSutils
	AnyError=1

if args.cores > 16: #If user asked for too many cores (unused if -splitJobs)
	eprint("Warning you requested too many cores for the common factot GWAS, are you sure thats necessary?")

try: #try reading input fail, error if couldn't read
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
for row in inputDataFrame.index: # For each of the input sumstats
	if inputDataFrame.loc[row,'binary']: # If its binary
		OLS.append('F') #then it was not OLS (for GSEM purposes)
		if inputDataFrame.loc[row,'linprob']: # If its binary and a linear probabiliy model was used (e.g. Bolt or any other linear regression)
			linprob.append('T') # save true for linprob
			selogit.append('F') # standard errors can't be in logistic scare
			prop.append(inputDataFrame.loc[row,'sampprev']) #store population and sample prevalences
			print("%s is a binary phenotype, but analyzed with a linear model, thus standard errors cannot be in the logistic scale. GenomicSEM will transform this betas to the logistic scale based on the sample prevalence\n"%(inputDataFrame.loc[row,'file']))
		else: #if binary bot not linprob (logistic)
			linprob.append('F') 
			if inputDataFrame.loc[row,'selogit']: #We need to know hether errors are in logistic scare or not, can't assume that
				selogit.append('T')
			else:
				selogit.append('F')
			prop.append('NA')
			print("%s is a binary phenotype, analyzed with a logistic model. Your input for selogit (%s) will be used to determine whether standard errors are in the logistic scale\n"%(inputDataFrame.loc[row,'file'],inputDataFrame.loc[row,'selogit']))

	else: #If its not binary the folllowing is true:
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
	Step1ID={} #we need this blank to enable compatiblity further on
else:
	print("Generating script one") #Below we are creating a string with the R script for script 1.
	script1="require(GenomicSEM);setwd('%s');files<-c("%(normalwd)
	script1+=','.join(["'%s'"%(i) for i in inputDataFrame.loc[:,'file']])
	script1+=");trait.names<-c("
	script1+=','.join(["'%s'"%(i) for i in inputDataFrame.loc[:,'name']])
	script1+=");"
	script1+="munge(files = files,hm3 = '/working/lab_nickm/adrianC/CommonData/w_hm3.noMHC.snplist',trait.names = trait.names, info.filter = 0.9, maf.filter = 0.01);"
	script1+="traits<-paste(trait.names,'.sumstats.gz',sep='');"
	script1+="sample.prev<-c("
	script1+=",".join(["NA" if np.isnan(i) else "%f"%(i) for i in inputDataFrame.loc[:,'sampprev']])
	script1+=");"
	script1+="population.prev<-c("
	script1+=",".join(["NA" if np.isnan(i) else "%f"%(i) for i in inputDataFrame.loc[:,'popprev']])
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
	#Below we create the actual script from the string that was constructed above
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
	#Only submit if submission is requested
	if args.nosub:
		print("Finished writing step 1. Flag 'nosub' detected, will not submit to the cluster")
	else:
		print("Finished writing step 1. Will now submit it to the cluster")
		finishedJobs={}
		print("Submitted step 1 at %s"%(datetime.datetime.now()))
		Step1ID=SubmitScript(scriptname) # SubmitScript can be explored in QIMRPBSutils python package
		
#Step one and two can be done in parallel so we continue:

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
	script2+=",".join(["%s"%(i) for i in prop])
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
		Step2ID=SubmitScript(scriptname)
		finishedJobs={**Step1ID,**Step2ID} #This table contains the jobIDs of both step1 and step2 jobs, to track their completion
		CheckCompletion(finishedJobs) #This function will track all the job completion
		anyErrors=CheckErrors() # Checking for errors using simple grep
		if len(anyErrors):
			print("The following errors were reported:\n")
			for currError in anyErrors:
				print(currError)
				print(currError)
				AnyError=1
		else:
			print("No obvious errors reported")

	if AnyError:
		eprint("Please double check the errors above and resubmit the jobs")
		exit()
		
		
######################################################################
#                 Step three common factor GWAS                      #
######################################################################

if args.skipstep3:
	print("Skipping step 3 as requested by user")
else:
	if args.splitJobs: #If the GWAS calculation steps will be paralelised (highly recommended) 
		print("Will split the common factor GWAS jobs into subjobs. This has shown great speed but some jobs never finish and need to be ran again")
		print("Creating out directory and logfiles directory")
		print("Creating values file based on combined sumstats number of lines")
		
		###################################### Values file for JOB Array ####################################################
		submitWC=subprocess.Popen("wc -l %s/%s_combinedsumstats.tsv"%(normalwd,args.jobname),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE) #Creating a working directory
		wcSTR=''.join([i.decode("utf-8") for i in submitWC.communicate()])
		reresult=int(re.search('(\d+)',wcSTR).group(1))-1
		start=[]
		stop=[]
		# Create blocks of SNPs to process, current size is 10,000 I found it a sweet spot for speed.
		for i in range(1,reresult,10000):
			start.append(i)
			stop.append(i+9999)
		stop[len(stop)-1]=reresult
		input=["%s/%s_combinedsumstats.tsv"%(normalwd,args.jobname) for i in range(len(start))]
		output=["%s/CommonFactorGWASout/%s_%s_%s.dat"%(normalwd,args.jobname,start[i],stop[i]) for i in range(len(start))]
		values=pd.DataFrame({"input":input,"start":start,"stop":stop,"output":output}) #create dataframe with columns to be sued for the job array
		values.to_csv("values",sep=' ',index=False,header=False) #save that table with the name  "values"
		subprocess.call(['mkdir','-p','CommonFactorGWASout','logfiles']) #create output directories and logfile directories
		
		###################################### Base script (R) for JOB Array ####################################################
		print("Generating script three (Rbase script)") #This is similar to steps one and two above, but it is a base script for a PBS job array
		script3="args = commandArgs(trailingOnly=TRUE)\ninput1 <- args[1]\nn_start <- args[2]\nn_stop <- args[3]\noutput <- args[4]\n"
		script3+="require(GenomicSEM)\nrequire(hms)\nrequire(data.table)\nsetwd('%s')\n"%(normalwd)
		script3+="Model_D_SS <- fread('%s/%s_combinedsumstats.tsv',header=TRUE)\n"%(normalwd,args.jobname)
		script3+="load('LDSCoutput_stage2.RData')\n"
		script3+="Model_D_SS_trunc<-Model_D_SS[paste0(n_start):paste0(n_stop),]\n"
		script3+="Model_D_input<-addSNPs(LDSCoutput_stage2, Model_D_SS_trunc, parallel = FALSE)\n"
		script3+="print('Starting the CommonFactorGWAS!')\n"
		script3+="results<-commonfactorGWAS(SNPs=NULL,covstruc=NULL,cores=1,Output=Model_D_input,estimation='%s')\n"%(args.estimation)
		script3+="print('Finished the CommonFactorGWAS!')\n"
		script3+="write.table(results, file =paste0(output), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\\t')\n"
		print("R script three looked like this:\n\n%s\n\nPlease check that it looks alright"%(script3))
		with open("3_GSEM_CommonFactorGWAS.R",'w') as currscript:
			currscript.write(script3)
			
		###################################### Submission script file for JOB Array ####################################################
		JobScriptHeader=r'#!/bin/bash'
		scriptname="GSEM_CommonFactor_Array.PBS" 
		#jobarray script for submission
		eprint("Creating Step 3 job script (%s).\n"%(scriptname))
		try: #writing the job arry submission script which will call "3_GSEM_CommonFactorGWAS.R" parametrised by the "values" file
			currscript=open(scriptname, 'w')
			currscript.write("%s\n"%(JobScriptHeader))
			currscript.write("#PBS -N CFGWAS_3\n")
			currscript.write("#PBS -l ncpus=1,mem=12GB,walltime=2:00:00\n")
			currscript.write("cd $PBS_O_WORKDIR\n")
			currscript.write("module load R/3.5.1\n")
			currscript.write("VARIABLES=($(awk -v var=$PBS_ARRAY_INDEX 'NR==var' values))\n")
			currscript.write("input=${VARIABLES[0]}\n")
			currscript.write("currStart=${VARIABLES[1]}\n")
			currscript.write("currEnd=${VARIABLES[2]}\n")
			currscript.write("output=${VARIABLES[3]}\n")
			currscript.write("Rscript 3_GSEM_CommonFactorGWAS.R $input $currStart $currEnd $output\n")
		except IOError:
			eprint("Error : could not write file %s !\n"%(scriptname))
		finally:
			currscript.close()
			
		###################################### First submission round  ############################################################
		if args.nosub:
			print("Finished writing step 3. Flag 'nosub' detected, will not submit to the cluster")
		else:
			print("Finished writing step 3. Will now submit it to the cluster")
			finishedJobs=SubmitJobArray(scriptname,"values","logfiles")
			print("At %s : Finished submitting jobs; waiting for job completion.\n"%(datetime.datetime.now()))
			CheckCompletion(finishedJobs,120)
			# Checking for errors using simple grep
			anyErrors=CheckErrorsArray("logfiles")
			if len(anyErrors):
				print("The following errors were reported:\n")
				for currError in anyErrors:
					print(currError)
					exit()
			else:
				print("No obvious errors reported")
				
			print("Checking that all jobs output. When paralelising genomicSEM some jobs do not finish on first attempt")
			ExpectedFiles=values.output.values
			outfiles=glob.glob('%s/CommonFactorGWASout/%s_*.dat'%(normalwd,args.jobname))
			toRun=[i for i in ExpectedFiles if not i in outfiles]
			
###################################### Iterative submission round  ############################################################
			if args.retry: #if split jobs was used, retry is strongly recommended. Jobs fail very often probably due to lavan needing to start a localhost which can cause collision in the same node
				count=1
				while len(toRun):
					count=count+1
					print("%s subjobs did not finish. Will try rerunning them until they are done; this might be indefinite (until this job is killed) if there are errors!"%(len(toRun)))
					valuesTmp=values.loc[values.output.isin(toRun),:]
					valuesTmp.to_csv("values%s"%(count),sep=' ',index=False,header=False)
					scriptnameTmp="%s_%s"%(scriptname,count)
					valuesnameTmp="values%s"%(count)
					subprocess.Popen("sed 's/values/values%s/g' %s > %s "%(count,scriptname,scriptnameTmp),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
					time.sleep(10)
					finishedJobs=SubmitJobArray(scriptnameTmp,valuesnameTmp,logdir="logfiles")
					print("At %s : Finished submitting jobs; waiting for job completion.\n"%(datetime.datetime.now()))
					CheckCompletion(finishedJobs,120)
					# Checking for errors using simple grep
					print("Checking whether all jobs finished running in this round")
					outfiles=glob.glob('%s/CommonFactorGWASout/%s_*.dat'%(normalwd,args.jobname))
					toRun=[i for i in ExpectedFiles if not i in outfiles]
					if count>9:
						eprint("Number of job resubmissions exceeded 10, this is probably due to an error with the subjobs. Please see the different values and submission files created")
						exit()
				# Checking for errors using simple grep
				anyErrors=CheckErrorsArray("logfiles")
				if len(anyErrors):
					print("The following errors were reported:\n")
					for currError in anyErrors:
						print(currError)
						exit()
				else:
					print("No obvious errors reported")
			else:
				if len(toRun):
					print("Some subjobs did not finish in time. Try adding the -retry flag to try running them several times.")
			print ("All common factor GWAS subjobs finished running. Merging results")
			finalDF=pd.concat((pd.read_csv(f,header=0,sep='\t') for f in ExpectedFiles))
			finalDF.to_csv("%s/%s_CFGWAS.dat"%(normalwd,args.jobname),sep="\t",index=False)
			

###################################### If not slipt jobs (has never finished; good luck try 24 hrs)  ############################################################
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
			currscript.write("#PBS -l ncpus=%s,mem=120GB,walltime=24:00:00\n"%(args.cores))
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
