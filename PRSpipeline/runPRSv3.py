#!/usr/bin/env python
#
# QIMR Genetic Epidemiology standard PRS calculation pipeline.
#
# Written by Adrian Campos, QIMR Genetic Epidemiology. ("AC")
# 11 Sept 2019 onward : also partially the work of Scott Gordon ("SG")
#
#
# Update 22 March 2021 [AC]
# Updated to support Release11. Testing is required 
#
#
#
# Update 15 March 2021 [AC]
# Added the possibility of unweighted polygenic risk scoring
# Added example for interactive running
#
#
# Update 02 Feb 2020 [AC]
# Added the possibility of using SBayesR (instead or aditionally to clumping + thresholding)
# Added example for interactive running
#
# Update 23 Jan 2020 [AC]
# Removed beta changing upon flipped alleles (PLINK accounts for that already; tested empirically)
#
# Update 11 Sept 2019 [SG]
# Started work adding extra and more flexible command-line options
#
# Update 28 May 2019 [AC]
# Works with release10 default release9 and doesnt require pvalue.ranges to exist 

##############################################################################################################################################################
#                  Use the next code to test interactively (in an ipython console, make sure you run everything but the Arg parsing section)                 #
##############################################################################################################################################################
# remember to module load python and start ipython (use the ipython command cpaste to copy and paste if the identation fails)
# Copy everything inside the ''' adjusting your parameters !!
'''
class z():
	def __init__(self):
		self.input='/working/lab_nickm/adrianC/Li2016QCed/jnj_antidepressant_JNJ_efficacy2-5.0jj.dat.gz.ctg_QCed.ctg' #gwa sumstats
		self.jobname='EsCitalopramResponse' #job name
		self.log=False #keep as False
		self.output='EsCitalopramResponse' #output name
		self.standardized=True #whether to std the PRS
		self.NAmarker=['.'] #NA markers, leave as it is
		self.eafield=None # effect allele column name (if not one of the canonical ones)
		self.oafield=None # other allele column name (if not one of the canonical ones)
		self.pvaluefield=None
		self.SNPfield=None
		self.chrfield=None
		self.BPfield=None
		self.betafield=None
		self.ORfield=None
		self.FREQfield=None
		self.SEfield=None
		self.Nfield=None
		self.assocstrand='unknown'
		self.dataset='Release10_HRCr1.1'
		self.analysischromosomes=range(1,22+1)
		self.minRsq=0.6
		self.minMAF=0.01
		self.maxMAF=0.5
		self.matchby='positionandalleles'
		self.rangesfile=None
		self.runmetadata=True
		self.runclumping=True
		self.runSBayesR=True
		self.runPRS=True
		self.compileblocks=True
		self.mergeCTandSB=True
		self.SBRexcludeMHC=True
		self.SBRhsq=None
		self.SBRseed=None
		self.runUnweightedPRS=False
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
# Handy function to print to stderr which is normal practice for log messages.
# AC 10/02/20: I dont think we should send ALL output to stderr (just errors; makes it easier to find the erros)
def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

######################################################################
# Check completion function for submitted jobs to the hpc cluster
# Input a job dictionary of job IDs and zeros (means they are running), and a predetermined wait time between cluster query
def CheckCompletion(jobDic,timew=60):
#   eprint("into CheckCompletion\n")
	ids=[i for i in jobDic.keys()] #get all job ids 
	eprint("At %s : Watching job ids=\n%s\n"%(datetime.datetime.now(),ids))
	submstr=str("%s %s 2>&1"%(qstatbinary,' '.join(ids))) #submit query for all jobs
#   eprint("submstr=\n%s\n"%(submstr))  
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
#       eprint('out=\n%s'%(out))
		# Split lines of output.
		OutLines=[i for i in re.split(pattern='\n',string=out) if not i=='']
		# Get ids of finished jobs. Count. Mark as complete.
		ncomplete=0
		for i in OutLines:
			reresult=re.search('qstat: (\d+(|\[\])\.hpc.*) Job has finished',i)
			if reresult != None:
				id=reresult.group(1)
				if jobDic[id] != 1:
					eprint("Job %s has finished.\n"%(id))
				jobDic[id] = 1
				ncomplete += 1
		# Report progress if anything has changed.
#       eprint("njobs=%s ncomplete=%s ncompleteprev=%s\n"%(njobs,ncomplete,ncompleteprev))
		if (firstiter==1) or (ncomplete != ncompleteprev):
			eprint("At %s : currently %s of %s job(s) have finished.\n"%(datetime.datetime.now(),ncomplete,njobs))
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

# Return the strand-reversed/flipped version of an allele [possibly multi-basepair]

# mapping function for basepairs
SingleBasepairStrandFlipMapping={'A':'T','a':'t','T':'A','t':'a','G':'C','g':'c','C':'G','c':'g','N':'N','n':'n'}
		
# .. single-basepair-specific version [returns 'None' if not convertible]
def ReverseStrand_1bp(allele):
	x=SingleBasepairStrandFlipMapping[allele]
	if x:
		return SingleBasepairStrandFlipMapping[allele]
	return None

# .. full function [returns 'None' if not convertible]
def ReverseStrand(allele) :
		if allele==None:
				return None
		l=len(allele)
		# single-basepair shortcut
		if l==1 :
				if allele in SingleBasepairStrandFlipMapping:
						return SingleBasepairStrandFlipMapping[allele]
				return None
		rev=''
		# multi-basepair (full) version
		for i in range(l): # reverse basepair order
				bp=allele[l-1-i]
				if bp in SingleBasepairStrandFlipMapping:
						rev=rev+SingleBasepairStrandFlipMapping[bp]
				else:
						return None
		return rev
# .. array of full alleles [each entry set to 'None' if not convertible]
def ReverseStrand_array(array) :
		return [ReverseStrand(i) for i in array]

#for i in ['A','C','G','T','AT','GC','AG','CT','N','AGCTNTTT','x']:
#        eprint("ReverseStrand('%s') = '%s'"%(i,ReverseStrand(i)))

#.. enable/disable features
def str2bool(v):
	if isinstance(v, bool):
	   return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')
		
######################################################################
#                        Start of main program                       #
######################################################################

eprint("By Adrian Campos @ QIMR Berghofer\nPartly based on Colodro Conde L, and Couvy-Duchesne et al bash examples\nplease let me know of any problems errors or modifications\nadrian.campos@qimrberghofer.edu.au\nUpgrades October 2019-January 2020 by Scott Gordon.\n")
eprint("This script is provided as is with no warranty and under under a CC license ...")
hostname=socket.gethostname()

# Store working directory

normalwd=os.getcwd() 

######################################################################
#                        Arg parsing                                 #
######################################################################


#SKIP THIS WHEN TESTING INTERACTIVELY! (RUN THE SECTIO FOR INTERACTIVE TESTING INSTEAD). !!!!!!!!!!!!!!!

parser = argparse.ArgumentParser(description="INPUT: summary statistics, should contain at least the following columns: SNP(Marker),CHR,BP,REF,ALT,EAF,beta,SE,pval other header names are accepted. SUMSTATS are assumed to be QCed already, specially for imputation. This script includes ambiguous strand and indel removal. This scripts tries to intrepret the header if the header column equivalents are not provided as input (see help below)",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input',help="Input file")
parser.add_argument('jobname',default='PRS',help="Job Name")
#.. PBS stuff
parser.add_argument('-l','--log',help="Create a log instead of output to std out", action = 'store_true')
parser.add_argument('-o','--output',help="Output file, Individual_PRS_$jobname by default")

#.. summary statistics results file
parser.add_argument('-namrk','--NAmarker',dest='NAmarker',nargs='+',default=['.'],help='Use this marker value(s) as missing value (separate them with a space), default is . use this flag if you missings are not nan or na (e.g. -99 -9 .) It is best to use this flag at the end as it receives an undefinite number of arguments')
parser.add_argument('-eafield','-ea','--effect_allele',dest='eafield',nargs='+',action='store',help='Field (header) name(s) for the effect allele. Use only if it is not one of the synonyms (check help)')
parser.add_argument('-oafield','-oa','--other_allele',dest='oafield',nargs='+',action='store',help='Field name(s) for the other (non_effect or reference) allele.')
parser.add_argument('-pvaluefield',nargs='+',action='store',help='Field name(s) for p-value (P)')
parser.add_argument('-SNPfield',nargs='+',action='store',help='Field name(s) for SNP/marker name (SNP')
parser.add_argument('-chrfield',nargs='+',action='store',help='Field name(s) for chromosome name (chr)')
parser.add_argument('-BPfield',nargs='+',action='store',help='Field name(s) for basepair position (BP)')
parser.add_argument('-betafield',nargs='+',action='store',help='Field name(s) for effect size (BETA)')
parser.add_argument('-ORfield',nargs='+',action='store',help='Field name(s) for Odds Ratio (OR)')
parser.add_argument('-FREQfield',nargs='+',action='store',help='Field name(s) for effec allele (!) frequency; required only for SBayesR')
parser.add_argument('-SEfield',nargs='+',action='store',help='Field name(s) for beta SE; required only for SBayesR')
parser.add_argument('-Nfield',nargs='+',action='store',help='Field name(s) for GWAS N (sample size); required only for SBayesR')
parser.add_argument('-assocstrand',action='store',default='unknown',choices=['+','-','unknown'],help='Strand alignment of summary statistic alleles, relative to target dataset ("+"=same; "-"=opposite; "unknown"=not known)')

#.. genotype dataset
parser.add_argument('-dataset',action='store',default='Release10_HRCr1.1',help="Run the pipeline using this dataset name (some combination with '_' separators, of 'Release9' or 'Release10'; reference panel 'HRCr1.1' or '1000GPhase3' or '1000GPhase1'; and/or chip family 'GSA', 'Hap', 'Hap610K', 'Omni', or 'merged' / 'merged610K')")

# .. marker filtering
parser.add_argument('-analysischromosomes',action='store',nargs='+',default=range(1,22+1),help='Which chromosomes to include in PRS (valid range integers from 1 through 22)')
parser.add_argument('-minRsq',action='store',default=0.6,type=float,help='Filter out markers with imputed Rsq below this threshold (valid range 0 to 1; 1 is no filtering)')
parser.add_argument('-minMAF',action='store',default=0.01,type=float,help='Filter out markers with MAF below this threshold (valid range 0 to 0.5; 0 is no filtering')
parser.add_argument('-maxMAF',action='store',default=0.5,type=float,help='Filter out markers with MAF above this threshold (not normally used; valid range 0 to 0.5; 0.5 is no filtering')

#.. marker name matching (apart from simple field names)
parser.add_argument('-matchby',action='store',choices=['rsID','name','positionandalleles','position'],default='positionandalleles',help="What to match markers by; positions must be the correct HG build if using 'position'. 'rsID' = 1st recorded rsID; 'name' = imputed name; 'positionandalleles' = by chromosome, position, +-strand alleles (eg. 1:23135623:A:G); 'position' = by chromosome, position only (disregard alleles; eg. 1:23135623) - not recommended if avoidable). Note that these options do not remove the need for strand alignment of scoring file to the target dataset.")

#.. pvalue-ranges file
parser.add_argument('-rangesfile',action='store',help='Use this file as the one listing pvalue ranges; will produce one PRS for each range in the file [default=./pvalue.ranges if it exists; otherwise a standard set of ranges]')

parser.add_argument('-std','--standardized',type=str2bool, default='True', action='store',help="Standardise all the PRS calculated")
parser.add_argument('-runmetadata', type=str2bool, default='True', action='store', help="Run metadata matching step ? ('True' or 'False'; only set to 'False' if nothing has changed in the input)")
parser.add_argument('-runclumping', type=str2bool, default='True', action='store', help="Run clumping step ? ('True' or 'False'; only set to 'False' if nothing has changed in the input), or if only SBayesR approach is to be used")
parser.add_argument('-runSBayesR', type=str2bool, default='True', action='store', help="Run SBayesR ? ('True' or 'False'; default False only clumping + thresholding is to be used. Running SBayesR requires rsNumbers as markers and the sumstats to have an effect allele frequency column")
parser.add_argument('-runPRS', type=str2bool, default='True', action='store', help="Run PRS calculation step ? ('True' or 'False'; only set to 'False' if nothing has changed in the input)")
parser.add_argument('-compileblocks', type=str2bool, default='True', action='store', help="Compile PRS results across the whole genome ? ('True' or 'False'; generally only set to 'False' for testing)")
parser.add_argument('-mergeCTandSB', type=str2bool, default='False', action='store', help="Compile SbayesR and Clumping+ Thresholding results in a final single file")
parser.add_argument('-runUnweightedPRS', type=str2bool, default='False', action='store', help="Run unweighted PRS ? ('True' or 'False'; default False. If true unweighted clumping + thresholding PRS is to be performed. This is additional to normal clumping + thresholding thus the flag is incompatible if clumping is False")

#.. SBayes R optional changes
parser.add_argument('-SBRexcludeMHC',type=str2bool,default='True', action='store',help="Exclude the MHC region from SBayes R, set to false to comapre with methods not excluding MHC region")
parser.add_argument('-SBRhsq',default=None, help="Optional SNP based heritability prior value")
parser.add_argument('-SBRseed',default=None, help="Optional seed for SBayesR, for debugging")

parser._actions[0].help="Print this help message and exit"
args = parser.parse_args()
 

######################################################################
#                        Initial checks                              #
######################################################################

RunMetadataMatch=args.runmetadata
RunClumping=args.runclumping
RunSBayesR=args.runSBayesR
RunPRSCalculation=args.runPRS
CompilePRSBlocks=args.compileblocks
MergeCTandSB=args.mergeCTandSB
JobNamePrefix=args.jobname
SBRexcludeMHC=args.SBRexcludeMHC
SBRhsq=args.SBRhsq
SBRseed=args.SBRseed
Unweighted=args.runUnweightedPRS


#print("RunMetadataMatch='%s'; RunClumping='%s'; RunPRSCalculation='%s'; JobNamePrefix='%s'\n"%(RunMetadataMatch,RunClumping,RunPRSCalculation,JobNamePrefix))

# Wish list prior to 11 Sept 2019 :
# Check that args are valid and check input files to be readable and interpretable. 
#e.g. missing command line articles and their compatibility.
#Pick up as many errors as once using a flag variable errorfound=0; exit at the end of checks
#eg. for clumping just mention the errors instead on each file.
#Skipping or switching sections e.g. if the clumping part ran successfully but checking the output is really there.
#Match by rs or by CHR:BP:A1:A2 

# [SG 11 Sept 2019]
# Basic checks at commencement of script, for various command-line arguments.
# Intended to prevent obvious cases where the script could run for a long time before
# producing errors, and this could have been detected earlier; eg. missing input files

AnyErrors=False

# .. check for specification/existence of input files

inputGWAS=args.input
rangesfile=args.rangesfile
GWASStrandAlignment=args.assocstrand
MarkerNameMatching=args.matchby
PRSname=args.jobname

if not inputGWAS:
	AnyErrors=true
	eprint('\nError : No summary statistics file supplied (option "input") !\n')
else:
	if not path.isfile(inputGWAS):
		AnyErrors=True
		eprint('\nError : Summary statistics file %s does not exist/is not a file !\n'%(str(inputGWAS)))

if not rangesfile:
	if path.isfile('pvalue.ranges'):
		rangesfile="pvalue.ranges"
		eprint('\nWarning : no ranges file supplied (option "-rangesfile") : using %s !\n'%(str(rangesfile)))
	else:
		eprint('\nWarning : no ranges file supplied (option "-rangesfile") - will use a default list of ranges !\n')

minRsqfilter=args.minRsq
if not(minRsqfilter):
	minRsqfilter=0.6
elif not((minRsqfilter >= 0) and (minRsqfilter <= 1)):
	AnyErrors=True
	eprint("\nError : minimum Rsq value %f (option '-minRsq') is out of range (allowed range 0 to 1) !\n"%(minRsq))

minMAFfilter=args.minMAF
if not(minMAFfilter):
	minMAFfilter=0.01
elif (minMAFfilter<0) or (minMAFfilter>0.5):
	AnyErrors=True
	eprint("\nError : minimum MAF value %f (option '-minMAF') is out of range (allowed range 0 to 0.5) !\n"%(minMAF))
maxMAFfilter=args.maxMAF
if not(maxMAFfilter):
	minMAFfilter=0.01
elif (maxMAFfilter<0) or (maxMAFfilter>0.5):
	AnyErrors=True
	eprint("\nError : maximum MAF value %f (option '-maxMAF') is out of range (allowed range 0 to 0.5) !\n"%(maxMAF))
if maxMAFfilter<minMAFfilter:
	AnyErrors=True
	eprint("\nError : maximum MAF value %f (option '-maxMAF') is less than minimum MAF value %f (option '-minMAF') !\n"%(maxMAF,minMAF))

if not((MarkerNameMatching == 'rsID') or (MarkerNameMatching == "name") or (MarkerNameMatching == "position") or (MarkerNameMatching == "positionandalleles")) :
	AnyErrors=True
	eprint('\nError : Marker name option %s is not yet implemented !'%(MarkerNameMatching))

ReleaseName='Release9'
ReferencePanelName="HRCr1.1"
ChipFamilyName="merged"
ReleaseNamecnt=0
ReferencePanelNamecnt=0
ChipFamilyNamecnt=0
JobScriptHeader=r'#!/bin/bash' # header for bash-like PBS job scripts
FileJoinerbinary="/reference/genepi/GWAS_release/Release10/Scripts/FileJoiner"
PLINKbinary="/working/genepi/software/bin/plink_1.90"
qsubbinary="/opt/pbs/bin/qsub"
qstatbinary="/opt/pbs/bin/qstat"
DatasetParts=args.dataset.split("_")

for sub in DatasetParts:
	lsub = sub.lower()
	Valid=False
	# Release name
	if (lsub == 'release10') or (sub == "r10"):
		ReleaseNamecnt += 1
		Valid=True
		ReleaseName='Release10'
	elif (lsub == "release9") or (sub == "r9"):
		ReleaseNamecnt += 1
		Valid=True
		ReleaseName='Release9'
	elif (lsub == "release11") or (sub == "r11"):
		ReleaseNamecnt += 1
		Valid=True
		ReleaseName='Release11'
	# Reference panel name
	elif (lsub == "hrcr1.1") or (sub == "hrc"):
		ReferencePanelNamecnt += 1
		Valid=True
		RerencePanelName="HRCr1.1"
	elif (lsub == "1000gphase3") or (sub == "1000gp3") or (sub == "1p3") or (sub == "1kg"):
		ReferencePanelNamecnt += 1
		Valid=True
		ReferencePanelName="1000GPhase3"
	elif (lsub == "1000gphase1") or (sub == "1000gp1") or (sub == "1p1"):
		ReferencePanelNamecnt += 1
		Valid=True
		ReferencePanelName="1000GPhase1"
	# Chip panel/dataset-merge name
	elif (lsub == "gsa"):
		ChipFamilyNamecnt += 1
		Valid=True
		ChipFamilyName="GSA"
	elif (lsub == "hap"):
		ChipFamilyNamecnt += 1
		Valid=True
		ChipFamilyName="Hap"
	elif (lsub == "hap610k"):
		ChipFamilyNamecnt += 1
		Valid=True
		ChipFamilyName="Hap610K"
	elif (lsub == "omni"):
		ChipFamilyNamecnt += 1
		Valid=True
		ChipFamilyName="Omni"
	elif (lsub == "merged"):
		ChipFamilyNamecnt += 1
		Valid=True
		ChipFamilyName="merged"
	elif (lsub == "merged610k"):
		ChipFamilyNamecnt += 1
		Valid=True
		ChipFamilyName="merged610K"
	# Error check [this component]
	if not(Valid):
		AnyErrors=True
		eprint("Error : unrecognised component '%s' in dataset name '%s' !\n"%(sub,dataset))

# Analysis chromosomes list

# .. check legitimacy
tmp_chr=set()
for chr in args.analysischromosomes:
	chr2=str(chr)
	if chr2 in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']:
		tmp_chr.add(chr2)
	elif chr2 in ['X','XY','Y']:
		AnyErrors=True
		eprint("Error : cannot currently use sex chromosomes (you requested '%s') !\n"%(chr2))
	elif chr2 in ['MT','MITO','M']:
		AnyErrors=True
		eprint("Error : cannot currently use mitochondrial markers (you requested '%s') !\n"%(chr2))
	else:
		eprint("Error : unknown/unsupported chromosome '%s' requested !\n"%(chr2))
		AnyErrors=True

# .. Convert to a suitably-sorted list of analysis chromosomes (regardless of user-supplied order)
AnalysisChromosomes=[]
for chr in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','XY','Y','MT']:
	if chr in tmp_chr:
		AnalysisChromosomes.append(chr)
if not(AnyErrors):
	eprint("Using these chromosomes in analysis : %s\n"%(' '.join(AnalysisChromosomes)))

# Error check
if ReleaseNamecnt > 1 :
	AnyErrors=True
	eprint("Error : Multiple GWAS release names found in dataset name '"+DatasetVersion+"' (expected one of 'Release9', 'Release10') !")
if ReferencePanelNamecnt > 1 :
	AnyErrors=True
	eprint("Error : Multiple reference panel names found in dataset name '"+DatasetVersion+"' (expected one of 'HRCr1.1'; '1000GPhase3'; '1000GPhase1') !")
if ChipFamilyNamecnt > 1 :
	AnyErrors=True
	eprint("Error : Multiple chip family names found in dataset name '"+DatasetVersion+"' (expected one of 'GSA', 'Hap', 'Omni', 'Hap610K' or 'merged' dataset names 'merged' or 'merged610K') !")


# determination of formal genotype dataset name; some basic related switches

IsMergedDataset=False
ChipFamilySuffix="_%s"%(ChipFamilyName)
if ChipFamilyName == "merged":
	IsMergedDataset=True
	ChipFamilySuffix=""
elif ChipFamilyName == "merged610K":
	IsMergedDataset=True
DatasetName='%s_%s%s'%(ReleaseName,ReferencePanelName,ChipFamilySuffix)
ReleaseLocation='/reference/genepi/GWAS_release/%s'%(ReleaseName)
DatasetLocation='%s/%s'%(ReleaseLocation,DatasetName)
metadata_file='%s/info/Metadata_allchr.txt'%(DatasetLocation)
BlockMarkerCountFile='%s/info/MarkerBlockCountsPerChromosome.txt'%(DatasetLocation)


eprint('\nUsing genotype dataset %s at location %s\n'%(DatasetName,DatasetLocation))
eprint("Using metadata file %s\n"%(metadata_file))
if not path.isfile(metadata_file):
	AnyErrors=True
	eprint("\nError : metadata file %s does not exist/is not a file !"%(metadata_file))


if AnyErrors :
		exit()



######################################################################
#            End of initial command-line checking stage              #
#            Start of GWA sumstats preprocessing                     #
######################################################################


if not RunClumping and RunSBayesR:
	eprint("Only SBayesR based scoring will be performed, no clumping or thresholding.")
	
if Unweighted:
	eprint("Unweighted PRS has also been requested. This is based on clumping and won't work if clumping is not performed or has been previously performed ")
	

if not rangesfile:
	rangesfile='pvalue.ranges'
	cutoffnames=['S1','S2','S3','S4','S5','S6','S7','S8']
	col2=[0,0,0,0,0,0,0,0]
	col3=[0.00000005,0.00001,0.001,0.01,0.05,0.1,0.5,1.0]
	tmp=pd.DataFrame(np.array([cutoffnames,col2,col3]))
	tmp=tmp.transpose()
	tmp.to_csv(rangesfile,sep=' ',index=False,header=None)
	eprint('The ranges file (%s) was not found inside the working directory. The classical QIMR ranges will be assumed.'%(str(rangesfile)))
else:
	try:
		rf=pd.read_csv(rangesfile,sep=' ',header=None)
	except Exception as e:
		eprint("Error : failed to read p-value ranges file %s !\nError was :\n%s"%(str(rangesfile),e))
		exit(1)
	cutoffnames=rf.iloc[:,0]
	eprint("Read these pvalue range/cutoff names : %s\n"%(' '.join(cutoffnames)))
	

###################File preprocessing ####################
#Open summary statistics

eprint("\nWill run a PRS analysis, using the following GWA summary statistics %s"%(inputGWAS))
eprint("Interpreting column headers")
colConverter={}

####### NOTE that in this code effect allele equals the allele that has an effect on the phenotype, whereas the reference allele is considered the basis (dosage of 0)

# Now has fully configurable field names [SG 11 Sept 2019]

HeaderSynonyms={}

if not args.pvaluefield:
		args.pvaluefield=['p','pvalue','p-value','pval','pvalues','p_bolt_lmm','mtag_pval']
for x in args.pvaluefield:
		HeaderSynonyms[x.lower()]='pval'

if not args.SNPfield:
	if MarkerNameMatching == 'rsID' :
		args.SNPfield=['snp','snpid','rsnumber','rs','snprs','rsid','marker','markername']
	else :
		args.SNPfield=['snp','snpid','marker','markername']
for x in args.SNPfield:
	HeaderSynonyms[x.lower()]='Marker'

if not args.eafield:
	args.eafield=['allele1','alt','a1','ea','effect_allele'] # SG 11Sept2019 : left out 'minor_allele' as not valid.
for x in args.eafield:
	HeaderSynonyms[x.lower()]='effect_allele'

if not args.oafield:
	args.oafield=['allele0','allele2','ref','a2','reference','oa','reference_allele']
for x in args.oafield:
	HeaderSynonyms[x.lower()]='reference_allele'

if not args.chrfield:
	args.chrfield=['chr','chromosome']
for x in args.chrfield:
	HeaderSynonyms[x.lower()]='CHR'

if not args.BPfield:
	args.BPfield=['BP','pos','position']
for x in args.BPfield:
	HeaderSynonyms[x.lower()]='BP'

if not args.betafield:
	args.betafield=['beta','effsize','logor','log(or)','mtag_beta','zscore','b']
for x in args.betafield:
	HeaderSynonyms[x.lower()]='beta'

if not args.ORfield:
	args.ORfield=['or']
for x in args.ORfield:
	HeaderSynonyms[x.lower()]='or'

if not args.FREQfield:
	args.FREQfield=['freq','maf','eaf','af']
for x in args.FREQfield:
	HeaderSynonyms[x.lower()]='freq' 

if not args.SEfield:
	args.SEfield=['se','stderr']
for x in args.SEfield:
	HeaderSynonyms[x.lower()]='se'

if not args.Nfield:
	args.Nfield=['n','samplesize']
for x in args.Nfield:
	HeaderSynonyms[x.lower()]='N'

# [SG 11 Sept 2019]
# "Do we need column ..."
needMarker=True
needeffect_allele=True
needreference_allele=True
needCHR=True
needBP=False
needbeta=True
needOR=False
needFREQ=False
needSE=False
needN=False
if MarkerNameMatching != 'name':
	needMarker=False
	needCHR=True
	needBP=True
	needeffect_allele=True
	needreference_allele=True

if RunSBayesR:
	needFREQ=True
	needMarker=True
	needSE=True
	needN=True
	
if MergeCTandSB and not RunSBayesR:
		eprint("Warning! You requested to merge clumping + thresholding and SBayesR results, but not to runSBayesR")
		
if MergeCTandSB and not RunClumping:
		eprint("Warning! You requested to merge clumping + thresholding and SBayesR results, but not to run clumping!")
# Test read of file

assocpreviewlength=1000              

eprint("\nTrying a test open of summary statistics file "+inputGWAS+"; will only read the first "+str(assocpreviewlength)+" records")

dtypeDict={'pval':object,'Marker':object,'effect_allele':object,'reference_allele':object,'CHR':object,'BP':object,'beta':object,'or':object,'freq':object,'se':object,'N':object}
numericCols=['pval','beta','or','freq','se','N']
#tmpDF=pd.read_table(inputGWAS,header=0,nrows=assocpreviewlength,sep=r"\s+")
tmpDF=pd.read_csv(inputGWAS,header=0,nrows=assocpreviewlength,sep=r"\s+")
colnames=[i for i in tmpDF.columns]
colConverter={i:HeaderSynonyms.get(i.lower(),None) for i in colnames}
colConverterRev={x:y for y,x in colConverter.items()}
ColsToUse=[i for i in colConverter if not colConverter[i]==None]
dtypes= {i:dtypeDict.get(colConverter[i],None) for i in colConverter}
toDel=[i for i in dtypes if dtypes[i]==None]
for key in toDel:
	del(dtypes[key])

eprint("\nInterpreted the columns in the following way:\n")
for i in colConverter:
	eprint(i,"treated as",colConverter[i])
eprint("")

# "Did we get field ... ?" and "Did we get two copies of field ... ?"
AnyErrors=False
haveMarker=False
haveeffect_allele=False
havereference_allele=False
haveCHR=False
haveBP=False
havebeta=False
haveOR=False
haveFREQ=False
haveSE=False
haveN=False
#.. count number of columns mapping to key 'key' // AC: below function could be a one liner return(len([i for i in colConverter.values() if i==value and i!=None]))
def colcounter(value):
	cnt=0
	for i in colConverter.values():
		if i != None:
			#eprint("value \'"+value+"\' i="+i)
			if i == value:
				cnt = cnt+1
	#eprint("colcounter : value \'"+value+"\' -> count "+str(cnt))
	return cnt
		
x=colcounter('Marker')
if x>0 :
	haveMarker=True
	if x>1 :
		AnyErrors=True
		eprint('Error : in summary statistics : found multiple \'marker name\' fields [Try option -SNPfield] !')
x=colcounter('effect_allele')
if x>0 :
	haveeffect_allele=True
	if x>1 :
		AnyErrors=True
		eprint('Error : in summary statistics : found multiple \'effect allele\' fields [Try option -eafield] !')
x=colcounter('reference_allele')
if x>0 :
	havereference_allele=True
	if x>1 :
		AnyErrors=True
		eprint('Error : in summary statistics : found multiple \'reference allele\' fields [Try option -oafield] !')
x=colcounter('CHR')
if x>0 :
	haveCHR=True
	if x>1 :
		AnyErrors=True
		eprint('Error : in summary statistics : found multiple \'chromosome name\' fields  [Try option -chrfield]!')
x=colcounter('BP')
if x>0 :
	haveBP=True
	if x>1 :
		AnyErrors=True
		eprint('Error : in summary statistics : found multiple \'basepair position allele\' fields  [Try option -BPfield]!')
x=colcounter('beta')
if x>0 :
	havebeta=True
	if x>1 :
		AnyErrors=True
		eprint('Error : in summary statistics : found multiple \'beta allele\' fields  [Try option -betafield]!')
x=colcounter('or')
if x>0 :
	haveOR=True
	if x>1 :
		AnyErrors=True
		eprint('Error : in summary statistics : found multiple \'odds ratio\' fields  [Try option -ORfield]!')
x=colcounter('freq')
if x>0 :
	haveFREQ=True
	if x>1 :
		AnyErrors=True
		eprint('Error : in summary statistics : found multiple \'freq allele\' fields (effect allele frequency)  [Try option -FREQfield]!')
x=colcounter('se')
if x>0 :
	haveSE=True
	if x>1 :
		AnyErrors=True
		eprint('Error : in summary statistics : found multiple \'SE allele\' fields (beta standard error)  [Try option -SEfield]!')
x=colcounter('N')
if x>0 :
	haveN=True
	if x>1 :
		AnyErrors=True
		eprint('Error : in summary statistics : found multiple \'N allele\' fields (beta standard error)  [Try option -Nfield]!')
transformEff=0
if haveOR:
	eprint("Detected an Odds ratio column, checking first %s elements"%(assocpreviewlength))
	if sum(tmpDF[colConverterRev['or']]<0):
		eprint("Detected at least one negative value, assuming OR to be log(OR)")
		transformEff=0
	else:
		eprint("Detected only positive values, will transform OR to log(OR)")
		transformEff=1
	key=[i for i in colConverter if colConverter[i]=='or'][0]
	colConverter[key]='beta'
	havebeta=True

# Do we have necessary columns ? Complain if not

AnyErrors=False
if needMarker and not(haveMarker):
	AnyErrors=True
	eprint('Error : In summary statistics : need a \'marker name\' column, but don\'t find it [Try option -SNPfield] !\n This can happen if you are matching by rsNumber or using SBayesR')
if needeffect_allele and not(haveeffect_allele):
	AnyErrors=True
	eprint('Error : In summary statistics : need a \'effect allele\' column, but don\'t find it [Try option -eafield] !')
if needreference_allele and not(havereference_allele):
	AnyErrors=True
	eprint('Error : In summary statistics : need a \'reference allele\' column, but don\'t find it [Try option -oafield] !')
if needCHR and not(haveCHR):
	AnyErrors=True
	eprint('Error : In summary statistics : need a \'chromosome name\' column, but don\'t find it [Try option -chrfield] !')
if needBP and not(haveBP):
	AnyErrors=True
	eprint('Error : In summary statistics : need a \'basepair position\' column, but don\'t find it [Try option -BPfield] !')
if needbeta and not(havebeta):
	AnyErrors=True
	eprint('Error : In summary statistics : need a \'beta\' column, but don\'t find it [Try option -betafield] !')
if needOR and not(haveOR):
	AnyErrors=True
	eprint('Error : In summary statistics : need a \'odds ratio\' column, but don\'t find it [Try option -ORfield] !')
if needFREQ and not(haveFREQ):
	AnyErrors=True
	eprint('Error : In summary statistics : need a \'FREQ\' column (because you called SBayesR), but don\'t find it [Try option -FREQfield] !')
if needSE and not(haveSE):
	AnyErrors=True
	eprint('Error : In summary statistics : need a \'SE\' column (because you called SBayesR), but don\'t find it [Try option -SEfield] !')
if needN and not(haveN):
	AnyErrors=True
	eprint('Error : In summary statistics : need an \'N\' column (because you called SBayesR), but don\'t find it [Try option -Nfield] !')
if AnyErrors:
		eprint('\nProblems were found with necessary columns in summary statistics file '+inputGWAS+' : please check the header and \'field name\' command-line options (use with \'-help\' for a list) !')
		exit()
if RunSBayesR:
	numRS=len([i for i in tmpDF[colConverterRev['Marker']] if str(i).startswith('rs')])
	if numRS < 500:
		eprint("SBayesR option was selected but very few inital SNPs (%s out of 1000) have rsNumber as markers"%(str(numRS)))
		AnyErrors=1

		

eprint("At %s : Reading summary statistics.\n"%(datetime.datetime.now()))
naMarkers=['NA','NAN','na','nan']
naMarkers.extend(args.NAmarker)
GWAsumstats = pd.read_csv(inputGWAS,header=0,index_col=None,usecols=ColsToUse,sep="\s+",dtype=dtypes,na_values=naMarkers,compression='infer')   
GWAsumstats.columns=[colConverter[i] for i in GWAsumstats.columns]
for col in GWAsumstats.columns:
	if col in numericCols:
		GWAsumstats.loc[:,col]=pd.to_numeric(GWAsumstats.loc[:,col],errors='coerce')

eprint("Read %s snps from the summary statistics"%(len(GWAsumstats.index)))

if sum(GWAsumstats['CHR']=='x') or sum(GWAsumstats['CHR']=='X') or sum(GWAsumstats['CHR']=='23'):
	eprint("X chromosome is not used by default (not currently implemented)")
	GWAsumstats=GWAsumstats.loc[GWAsumstats['CHR'].str.lower()!='x']
	GWAsumstats=GWAsumstats.loc[GWAsumstats['CHR']!='23']
	eprint("%s variants remaining in sumstats after X chr removal"%(len(GWAsumstats.index)))

# Convert CHR and BP o int to remove 0s
GWAsumstats['CHR']=pd.to_numeric(GWAsumstats.CHR,errors='coerce').astype(pd.Int64Dtype()).astype(str)
GWAsumstats['BP']=pd.to_numeric(GWAsumstats.BP,errors='coerce').astype(pd.Int64Dtype()).astype(str)


GWAsumstats['SNPvariant']=GWAsumstats['CHR'].apply(str)+':'+GWAsumstats['BP'].apply(str)

if transformEff:
	if sum(GWAsumstats['beta']<0):
		eprint('Found negative values, assuming OR already transformed')
	else:
		eprint("Transforming OR to log(OR) using np.log ")
		GWAsumstats['beta']=np.log(GWAsumstats['beta'])

## Convert allele names to upper case
GWAsumstats['reference_allele']=GWAsumstats.loc[:,'reference_allele'].str.upper()
GWAsumstats['effect_allele']=GWAsumstats.loc[:,'effect_allele'].str.upper()

######################################################################
#                         Sumstats QC steps                          #
######################################################################

#######Indel removal
#AC consider whether this is still OK to do? I think we could score with indels

eprint("\nRemoving indels from sumstats by keeping variants where: ref allele = [A|T|C|G]")
canonBases=['A','T','C','G']
canonicalBases=[True if i in canonBases else False for i in GWAsumstats['reference_allele']]
GWAsumstats=GWAsumstats.loc[canonicalBases,:]
if MarkerNameMatching == 'rsID' :
		eprint("\nRemoving markers with names not starting with 'rs'.")
		startsWithrs=[True if re.match('rs',str(i)) else False for i in GWAsumstats['Marker']]
		GWAsumstats=GWAsumstats.loc[startsWithrs,:]
eprint("%s snps left after indels removal"%(len(GWAsumstats.index)))


####### Ambiguous strand removal; flipping to '+' if possible at this stage.
# We create a second set of effect, reference alleles which are on opposite strand; if strand is unknown.

assocstrand=args.assocstrand
if args.assocstrand == '-':
	eprint("\nReversing '-' strand markers to '+' strand [or 'None' allele values]")
	GWAsumstats.effect_allele = ReverseStrand_array(GWAsumstats.effect_allele)
	GWAsumstats.reference_allele = ReverseStrand_array(GWAsumstats.reference_allele)
	GWAsumstats['effect_allele_2'] = [None for i in range(len(GWAsumstats.Marker))]
	GWAsumstats['reference_allele_2'] = [None for i in range(len(GWAsumstats.Marker))]
	assocstrand='+'

elif args.assocstrand == '+':
	eprint("\nSummary statistics are already on '+' strand : leaving alleles unchanged.")
	GWAsumstats['effect_allele_2'] = [None for i in range(len(GWAsumstats.Marker))]
	GWAsumstats['reference_allele_2'] = [None for i in range(len(GWAsumstats.Marker))]

else:
	eprint("\nFiltering strand-ambiguous markers from summary statistics")
	GWAsumstats['effect_allele_2'] = ReverseStrand_array(GWAsumstats.effect_allele)
	GWAsumstats['reference_allele_2'] = ReverseStrand_array(GWAsumstats.reference_allele)
	nonambig=(GWAsumstats.effect_allele != GWAsumstats.reference_allele_2)
	GWAsumstats=GWAsumstats.loc[nonambig & (GWAsumstats.reference_allele != None)]
	eprint("%s SNPs remaining after strand ambiguity drops."%(len(GWAsumstats.index)))

######## Duplicates removal
tmp=GWAsumstats.SNPvariant.duplicated()
if sum(tmp) >0:
	eprint("\nDetected %s duplicates"%(sum(tmp)))
	eprint("One of the duplicates of each marker in sumstats will be removed")
	GWAsumstats=GWAsumstats.loc[tmp==False,:]
	eprint("Remaining %s SNPs after duplicates removal"%(len(GWAsumstats)))


FNULL = open(os.devnull, 'w') # to redirect HPC output

if args.log:
	log = open("RunPRS.log", "a")
	sys.stdout = log
	sys.stderr = log



######################################################################
#                        SBayes R multivariate GWAS                  #
######################################################################

if RunSBayesR:
	eprint("Will obtain estimated multivariate SNP effects using SBayesR")
	subprocess.call(['mkdir','-p','SBayesR/logfiles']) #create SBayesR working directoy
	SBRsumstatsName='SBayesR/%s_SBayesR.inputsumstats'%(PRSname) #SBR input
	SBayesRinDF=GWAsumstats.loc[:,['Marker','effect_allele','reference_allele','freq','beta','se','pval','N']] #All these columns required for SBYR
	SBayesRinDF.dropna(inplace=True) #NAs seriously affect SBR (stops reading from first encounter, silently does this
	SBayesRinDF.columns=['SNP','A1','A2','freq','b','se','p','N'] #Header for SBR gcta style
	SBayesRinDF.to_csv(SBRsumstatsName,sep=' ',index=False,float_format='%.6f')
	scriptname="SBayesR/SUBMIT_SBayesR.PBS" 
	#jobarray script for submission
	eprint("Creating SBayesR job script (%s).\n"%(scriptname))
	SBRexcludeMHCline=''
	if SBRexcludeMHC:
		SBRexcludeMHCline='--exclude-mhc'
	SBRhsqline='--hsq %s'%(SBRhsq)
	if SBRhsq is None:
		SBRhsqline=''
	SBRseedline='--seed %s'%(SBRseed)
	if SBRseed is None:
		SBRseedline=''
	try:
		currscript=open(scriptname, 'w')
		currscript.write("%s\n"%(JobScriptHeader))
		currscript.write("#PBS -N SbayesR\n")
		currscript.write("#PBS -l mem=260GB,walltime=24:00:00,chip=Intel\n")
		currscript.write("cd $PBS_O_WORKDIR\n")
		currscript.write("chr=$PBS_ARRAY_INDEX\n")
		currscript.write("module load gctb/2.03beta\n")
		currscript.write("gctb --sbayes R \\\n--mldm /reference/genepi/public_reference_panels/ldm_ukb_50k_bigset_2.8M/ukb50k_2.8M_shrunk_sparse.mldmlist \\\n--pi 0.95,0.02,0.02,0.01 \\\n--gamma 0.0,0.01,0.1,1 \\\n--gwas-summary %s \\\n--chain-length 25000 \\\n--burn-in 5000 \\\n--out-freq 10 \\\n--out SBayesR/%s_SBayesR.output %s %s %s\n"%(SBRsumstatsName,PRSname,SBRexcludeMHCline,SBRhsqline,SBRseedline))
	except IOError:
		eprint("Error : could not write file %s !\n"%(scriptname))
		AnyErrors=True
	finally:
		currscript.close()
	if AnyErrors:
		exit(1)
		
	#Submit SBayesR job
	finishedJobs={} #hash to keep track of submission and check if done
	eprint("At %s : Submitting one array job.\n"%(datetime.datetime.now()))
	subjobname="SbayesR"
	Submit=subprocess.Popen('%s -N "%s" -o SBayesR/logfiles/ -e SBayesR/logfiles/ %s'%(qsubbinary,subjobname,scriptname),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	#Maybe all the following lines (until 860) should be wrapped in a submission function or object??
	outraw=Submit.communicate()
	outlines=outraw[0].decode("utf-8") + outraw[1].decode("utf-8")
	if outlines.find('qsub') >=0:
		eprint('qsub returned this message :\n%s\n'%(outlines))
		AnyErrors=True
	elif outlines=='':
		eprint('qsub returned nothing.\n')
		AnyErrors=True
	else :
		jobID=outlines.rstrip()
		eprint('.. job name is %s'%(jobID))
		finishedJobs[jobID]=0
	
	if AnyErrors:
		eprint('Error(s) in PBS job submission. Quiting and cancelling submitted job(s).\n')
		for i in finishedJobs:
				subprocess.Popen('qdel %s'%(i),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		exit(1)
	
	eprint("At %s : Job(s) submitted; waiting for job completion prior to continuing.\n"%(datetime.datetime.now()))


	# wait for job completion
	eprint("At %s : Finished submitting jobs; waiting for job completion.\n"%(datetime.datetime.now()))
	CheckCompletion(finishedJobs,120)
	
	
	# Checking for errors using simple grep
	AnyErrors=False
	logs=glob.glob('SBayesR/logfiles/*.ER')
	for i in logs:
		Submits=subprocess.Popen("grep 'Error' %s"%(i),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		errlines=Submits.communicate()[0].decode("utf-8")
		Errors=[i for i in re.split(pattern='\n',string=errlines) if not i=='']
		if len(Errors) > 0:
			eprint("Log file %s contains these error(s) :\n%s\n"%(i,'\n'.join(Errors)))
			AnyErrors=True
	if AnyErrors:
		exit(1)
	else:
		eprint("No errors reported in SBayesR stage.\n")
		
		
SBayesRFiles=glob.glob('SBayesR/*.snpRes')
if len(SBayesRFiles)>0:
	eprint("Compiling SBayesR results this may take a while %s"%(datetime.datetime.now()))
	SBayesROutDF=pd.concat((pd.read_csv(f,header=0,sep='\s+',usecols=['Name','A1','A1Effect']) for f in SBayesRFiles))
	SBayesROutDF.columns=['Marker','sbrA1','sbrBETA']
	eprint("Finished compiling SBayesR at %s"%(datetime.datetime.now()))
	GWAsumstats=pd.merge(GWAsumstats,SBayesROutDF,on='Marker',how='outer')
else:
	eprint("SBayesR was not run in this nor in a previous run of this program. No SBayesR based scoring will occur. (If you wanted SBayesR rerun  using -runSBayesR flag)")
	MergeCTandSB=False


######################################################################
#                Match sumstats to metadata and perform QC           #
######################################################################

# Run script to do metadata work (basically a combination of 'awk' commands and Scott Gordon's FileJoiner utility)

Rsqcol="Rsq"
if IsMergedDataset:
		Rsqcol="Rsq_rederived"

# 1. Export summary statistics in the right form

procsummfilename='./summstats_filtered_raw.txt'
markercolindex=GWAsumstats.columns.get_loc('Marker')+1
betacolindex=GWAsumstats.columns.get_loc('beta')+1
eacolindex=GWAsumstats.columns.get_loc('effect_allele')+1
oacolindex=GWAsumstats.columns.get_loc('reference_allele')+1
ea2colindex=GWAsumstats.columns.get_loc('effect_allele_2')+1
oa2colindex=GWAsumstats.columns.get_loc('reference_allele_2')+1

if not RunMetadataMatch :
	eprint("Skipping metadata matching stage, by request.")
else :
	metascriptname="./makemeta.sh"
	metascr=open(metascriptname,"w")
	if not metascr:
			eprint("Error : could not create meta-data processing script !")
			exit(1)
	metascr.write("#!/bin/bash\n")
	metascr.write('echo "Deleting old tmp metadata files"\n')
	metascr.write('rm -f tmp_meta1 tmp_meta2 alignmentresults.txt\n')

	if MarkerNameMatching == 'rsID': # marker name is a dbSNP name/rsID; need to match on all known possibilities for each marker.
			eprint('Writing temporary filtered summary statistics for matching to metadata (file %s)\n'%(procsummfilename))
			GWAsumstats.to_csv(procsummfilename,header=True,sep='\t',index=False)
			# QC filter; write out just the relevant fields after appropriate filtering
			# field list is : rsID dataset_marker_name R9name bp_Build37 REF ALT MAF Rsq (tab-delimited)
			# Assumes that the various rsID fields are in a continuous run starting at field "Marker_dbSNP_1" or "SNP_dbSNP"

			metascr.write('echo "Creating metadata file, QC-filtered and indexed by rsID (file tmp_meta1)"\n')
			metascr.write('awk -v FS="\\t" -v OFS="\\t" -v minRsq=%s -v minMAF=%s -v maxMAF=%s -v Rsqcol=%s \'BEGIN{Nfilteredout=0;NnorsID=0;Ndropped=0;Nkept=0}\n(FNR==1){for(i=1;i<=NF;++i){if(($i=="Marker")||($i=="SNP")){Markerfld=i};if($i==Rsqcol){Rsqfld=i};if($i=="MAF"){MAFfld=i};if(($i=="REF")||($i=="REF(0)")){REFfld=i};if(($i=="ALT")||($i=="ALT(1)")){ALTfld=i};if($i=="CHR"){CHRfld=i};if($i=="bp_Build37"){BPfld=i};if(($i=="Marker_dbSNP_1")||($i=="SNP_dbSNP")){firstdbSNP=i};if((substr($i,1,12)=="Marker_dbSNP")||(substr($i,1,9)=="SNP_dbSNP")){lastdbSNP=i}}}\n((FNR>1)&&(Rsqfld>=1)&&(MAFfld>=1)&&(firstdbSNP>=1)&&(lastdbSNP>=1)){\nMAF=($MAFfld>0.5)?(1-$MAFfld):$MAFfld; passesfilter=(($Rsqfld>=minRsq)&&(MAF>=minMAF)&&(MAF<=maxMAF))?1:0; if(!passesfilter){++Nfilteredout}; anydbSNP=($firstdbSNP==".")?0:1; if(!anydbSNP){++NnorsID}\nif(passesfilter && anydbSNP){++Nkept;for(i=firstdbSNP;i<=lastdbSNP;++i){if($i!="."){print $i,$Markerfld,$CHRfld":"$BPfld":"(($REFfld<$ALTfld)?($REFfld":"$ALTfld):($ALTfld":"$REFfld)),$CHRfld,$BPfld,$REFfld,$ALTfld,$MAFfld,$Rsqfld}}}else{++Ndropped}\n}\nEND{print "Dropped "Nfilteredout" markers due to filtering, and "NnorsID" due to lack of rsID; "Ndropped" total; kept "Nkept > "/dev/stderr"}\' %s > tmp_meta1\n'%(minRsqfilter,minMAFfilter,maxMAFfilter,Rsqcol,metadata_file))
			# uniqueness filter
			metascr.write('echo "Extracting only rsIDs found just once in metadata (file tmp_meta2)"\n')
			metascr.write('(echo -e "rsID\tdataset_Marker\tR9name\tdataset_CHR\tdataset_bp\tREF\tALT\tMAF\tRsq"\n')
			metascr.write('sort -k1 tmp_meta1 | awk -v FS="\\t" -v OFS="\\t" \'BEGIN{Ndropped=0;prevtxt="";prevname="";N=0}(FNR==1){print $0}(FNR>1){if(prevname==$1){prevtxt=(N>0)?prevtxt"\\n"$0:$0;++N}else{if(N==1){print prevtxt}else{++Ndropped};N=1;prevtxt=$0};prevname=$1}END{if(N==1){print prevtxt}else{++Ndropped};print "Dropped "Ndropped" duplicate rsIDs" > "/dev/stderr" }\'\n')
			metascr.write(')> tmp_meta2\n')

	elif MarkerNameMatching == 'name' : # marker name is as used in the GWAS dataset, so no lookup or conversion to do.
			eprint('Writing temporary filtered summary statistics for matching to metadata (file %s)\n'%(procsummfilename))
			GWAsumstats.to_csv(procsummfilename,header=True,sep='\t',index=False)
			# QC filter; write out just the relevant fields after appropriate filtering
			# field list is : marker dataset_marker_name bp_Build37 REF ALT MAF Rsq (tab-delimited)
			# Produces directly the output which has to go via a second 'uniqueness' filter for the rsID version.

			metascr.write('echo "Creating metadata file, QC-filtered and indexed by dataset marker name (file tmp_meta2)"\n')
			metascr.write('awk -v FS="\\t" -v OFS="\\t" -v minRsq=%s -v minMAF=%s -v maxMAF=%s -v Rsqcol=%s \'BEGIN{Nfilteredout=0;NnorsID=0;Ndropped=0;Nkept=0; print "Marker","dataset_Marker","R9name","dataset_CHR","dataset_bp","REF","ALT","MAF","Rsq"}\n(FNR==1){for(i=1;i<=NF;++i){if(($i=="Marker")||($i=="SNP")){Markerfld=i};if($i==Rsqcol){Rsqfld=i};if($i=="MAF"){MAFfld=i};if(($i=="REF")||($i=="REF(0)")){REFfld=i};if(($i=="ALT")||($i=="ALT(1)")){ALTfld=i};if($i=="CHR"){CHRfld=i};if($i=="bp_Build37"){BPfld=i}}}\n((FNR>1)&&(MAFfld>=1)){\nMAF=($MAFfld>0.5)?(1-$MAFfld):$MAFfld; passesfilter=(($Rsqfld>=minRsq)&&(MAF>=minMAF)&&(MAF<=maxMAF))?1:0; if(!passesfilter){++Nfilteredout}\nif(passesfilter){++Nkept;print $Markerfld,$Markerfld,$CHRfld":"$BPfld":"(($REFfld<$ALTfld)?($REFfld":"$ALTfld):($ALTfld":"$REFfld)),$CHRfld,$BPfld,$REFfld,$ALTfld,$MAFfld,$Rsqfld}else{++Ndropped}\n}\nEND{print "Dropped "Nfilteredout" markers due to filtering; kept "Nkept > "/dev/stderr"}\' %s > tmp_meta2\n'%(minRsqfilter,minMAFfilter,maxMAFfilter,Rsqcol,metadata_file))

	elif MarkerNameMatching == 'positionandalleles' : # marker name is chr:position:A1:A2 where A1, A2 are '=' -strand REF, ALT alleles in some order (Build 37)
			eprint('Writing temporary filtered summary statistics for matching to metadata (file %s)\n'%(procsummfilename))
			if RunSBayesR:
				GWAsumstats['rsID']=GWAsumstats['Marker'].copy()
			GWAsumstats['Marker']=GWAsumstats.CHR+':'+GWAsumstats.BP+':'+GWAsumstats.effect_allele+':'+GWAsumstats.reference_allele
			GWAsumstats.to_csv(procsummfilename,header=True,sep='\t',index=False)
			# QC filter; write out just the relevant fields after appropriate filtering
			# field list is : marker dataset_marker_name R9name marker bp_Build37 REF ALT MAF Rsq (tab-delimited)
			# with a 2nd copy of each marker with the alleles swapped, to allow matching on both possibilities.
			# Produces directly the output which has to go via a second 'uniqueness' filter for the rsID version.

			metascr.write('echo "Creating metadata file, QC-filtered and indexed by chr:position:A1:A2 marker name (file tmp_meta2)"\n')
			metascr.write('awk -v FS="\\t" -v OFS="\\t" -v minRsq=%s -v minMAF=%s -v maxMAF=%s -v Rsqcol=%s \'BEGIN{Nfilteredout=0;NnorsID=0;Ndropped=0;Nkept=0; print "Marker","dataset_Marker","R9name","dataset_CHR","dataset_bp","REF","ALT","MAF","Rsq"}\n(FNR==1){for(i=1;i<=NF;++i){if(($i=="Marker")||($i=="SNP")){Markerfld=i};if($i==Rsqcol){Rsqfld=i};if($i=="MAF"){MAFfld=i};if(($i=="REF")||($i=="REF(0)")){REFfld=i};if(($i=="ALT")||($i=="ALT(1)")){ALTfld=i};if($i=="CHR"){CHRfld=i};if($i=="bp_Build37"){BPfld=i}}}\n((FNR>1)&&(MAFfld>=1)){\nMAF=($MAFfld>0.5)?(1-$MAFfld):$MAFfld; passesfilter=(($Rsqfld>=minRsq)&&(MAF>=minMAF)&&(MAF<=maxMAF))?1:0; if(!passesfilter){++Nfilteredout}\nif(passesfilter){++Nkept;R9name=$CHRfld":"$BPfld":"(($REFfld<$ALTfld)?($REFfld":"$ALTfld):($ALTfld":"$REFfld));altname=$CHRfld":"$BPfld":"(($REFfld<$ALTfld)?($ALTfld":"$REFfld):($REFfld":"$ALTfld));print R9name,$Markerfld,R9name,$CHRfld,$BPfld,$REFfld,$ALTfld,$MAFfld,$Rsqfld;print altname,$Markerfld,R9name,$CHRfld,$BPfld,$REFfld,$ALTfld,$MAFfld,$Rsqfld}else{++Ndropped}\n}\nEND{print "Dropped "Nfilteredout" markers due to filtering; kept "Nkept > "/dev/stderr"}\' %s > tmp_meta2\n'%(minRsqfilter,minMAFfilter,maxMAFfilter,Rsqcol,metadata_file))

	elif MarkerNameMatching == 'position' : # marker name is chr:position (Build 37)
			eprint('Writing temporary filtered summary statistics for matching to metadata (file %s)\n'%(procsummfilename))
			if RunSBayesR:
				GWAsumstats['rsID']=GWAsumstats['Marker'].copy()
			GWAsumstats['Marker']=GWAsumstats.CHR+':'+GWAsumstats.BP
			GWAsumstats.to_csv(procsummfilename,header=True,sep='\t',index=False)
			# QC filter; write out just the relevant fields after appropriate filtering
			# field list is : marker dataset_marker_name R9name bp_Build37 REF ALT MAF Rsq (tab-delimited)
			# with a 2nd copy of each marker with the alleles swapped, to allow matching on both possibilities.
			# Needs a second 'uniqueness' filter as per rsID version.

			metascr.write('echo "Creating metadata file, QC-filtered and indexed by dataset marker name (file tmp_meta1)"\n')
			metascr.write('awk -v FS="\\t" -v OFS="\\t" -v minRsq=%s -v minMAF=%s -v maxMAF=%s -v Rsqcol=%s \'BEGIN{Nfilteredout=0;NnorsID=0;Ndropped=0;Nkept=0}\n(FNR==1){for(i=1;i<=NF;++i){if(($i=="Marker")||($i=="SNP")){Markerfld=i};if($i==Rsqcol){Rsqfld=i};if($i=="MAF"){MAFfld=i};if(($i=="REF")||($i=="REF(0)")){REFfld=i};if(($i=="ALT")||($i=="ALT(1)")){ALTfld=i};if($i=="CHR"){CHRfld=i};if($i=="bp_Build37"){BPfld=i}}}\n((FNR>1)&&(MAFfld>=1)){\nMAF=($MAFfld>0.5)?(1-$MAFfld):$MAFfld; passesfilter=(($Rsqfld>=minRsq)&&(MAF>=minMAF)&&(MAF<=maxMAF))?1:0; if(!passesfilter){++Nfilteredout}\nif(passesfilter){++Nkept;R9name=$CHRfld":"$BPfld":"(($REFfld<$ALTfld)?($REFfld":"$ALTfld):($ALTfld":"$REFfld));print $CHRfld":"$BPfld,$Markerfld,R9name,$CHRfld,$BPfld,$REFfld,$ALTfld,$MAFfld,$Rsqfld}else{++Ndropped}\n}\nEND{print "Dropped "Nfilteredout" markers due to filtering; kept "Nkept > "/dev/stderr"}\' %s > tmp_meta1\n'%(minRsqfilter,minMAFfilter,maxMAFfilter,Rsqcol,metadata_file))

			# uniqueness filter
			metascr.write('echo "Extracting only chr:position names found just once in metadata (file tmp_meta2)"\n')
			metascr.write('(echo -e "Marker\tdataset_Marker\tR9name\tdataset_CHR\tdataset_bp\tREF\tALT\tMAF\tRsq"\n')
			metascr.write('sort -k1 tmp_meta1 | awk -v FS="\\t" -v OFS="\\t" \'BEGIN{Ndropped=0;prevtxt="";prevname="";N=0}(FNR==1){print $0}(FNR>1){if(prevname==$1){prevtxt=(N>0)?prevtxt"\\n"$0:$0;++N}else{if(N==1){print prevtxt}else{++Ndropped};N=1;prevtxt=$0};prevname=$1}END{if(N==1){print prevtxt}else{++Ndropped};print "Dropped "Ndropped" duplicate chr:position names" > "/dev/stderr" }\'\n')
			metascr.write(')> tmp_meta2\n')

	# match to summary statistics by name
	# adds strand alignment and REF:ALT alignment calls (REFALT: effect_allele=REF, reference_allele=ALT; ALTREF: reverse)
	metascr.write('summstatfile=%s\n'%(procsummfilename))
	metascr.write('FileJoiner="%s"\n'%(FileJoinerbinary))
	metascr.write('echo "Matching summary statistics (file ${summstatfile}) to metadata"\n')
#   metascr.write('${FileJoiner} -quiet ${summstatfile},sep=tabs,key=%s,headerline tmp_meta2,sep=tabs,headerline,key=1 | awk -v FS="\\t" -v OFS="\\t" -v betacol=%s -v eac=%s -v oac=%s -v ea2c=%s -v oa2c=%s  \'(FNR==1){print $0,"aligned_beta","Strand_alignment","REFALT_alignment","aligned_effect_allele","aligned_other_allele"}(FNR>1){ REF=$(NF-3);ALT=$(NF-2);matchplus1=(($eac==REF)&&($oac==ALT));matchplus2=(($eac==ALT)&&($oac==REF));matchminus1=(($ea2c==REF)&&($oa2c==ALT));matchminus2=(($ea2c==ALT)&&($oa2c==REF));matchplus=(matchplus1||matchplus1);matchminus=(matchminus1||matchminus2);strandalign=(matchplus&&matchminus)?"?":((matchplus&&!matchminus)?"+":((matchminus&&!matchplus)?"-":"none")); match1=matchplus1||matchminus1;match2=matchplus2||matchminus2;refalign=(match1&&match2)?"?":((match1&&!match2)?"REFALT":((match2&&!match1)?"ALTREF":"none")); modbeta="."; if(refalign=="ALTREF"){modbeta=$betacol}else if(refalign="REFALT"){modbeta=-$betacol};print $0,modbeta,strandalign,refalign,(refalign=="REFALT")?REF:((refalign=ALTREF)?ALT:"."),(refalign=="REFALT")?ALT:((refalign==ALREF)?REF:".")}\' > alignmentresults.txt\n'%(markercolindex,betacolindex,eacolindex,oacolindex,ea2colindex,oa2colindex))
	metascr.write('${FileJoiner} -quiet ${summstatfile},sep=tabs,key=%s,headerline tmp_meta2,sep=tabs,headerline,key=1 | awk -v FS="\\t" -v OFS="\\t" -v betacol=%s -v eac=%s -v oac=%s -v ea2c=%s -v oa2c=%s  \'(FNR==1){print $0,"aligned_beta","Strand_alignment","REFALT_alignment","aligned_effect_allele","aligned_other_allele"}(FNR>1){ REF=$(NF-3);ALT=$(NF-2);matchplus1=(($eac==REF)&&($oac==ALT));matchplus2=(($eac==ALT)&&($oac==REF));matchminus1=(($ea2c==REF)&&($oa2c==ALT));matchminus2=(($ea2c==ALT)&&($oa2c==REF));matchplus=(matchplus1||matchplus2);matchminus=(matchminus1||matchminus2);strandalign=(matchplus&&matchminus)?"?":((matchplus&&!matchminus)?"+":((matchminus&&!matchplus)?"-":"none")); match1=matchplus1||matchminus1;match2=matchplus2||matchminus2;refalign=(match1&&match2)?"?":((match1&&!match2)?"REFALT":((match2&&!match1)?"ALTREF":"none")); modbeta="."; if(refalign=="ALTREF"){modbeta=$betacol}else if(refalign=="REFALT"){modbeta=$betacol};print $0,modbeta,strandalign,refalign,(refalign=="REFALT")?REF:((refalign=="ALTREF")?ALT:"."),(refalign=="REFALT")?ALT:((refalign=="ALTREF")?REF:".")}\' > alignmentresults.txt\n'%(markercolindex,betacolindex,eacolindex,oacolindex,ea2colindex,oa2colindex))
	metascr.write('echo "Filtering out markers with no/ambiguous unique alignment to target markers (by strand or REF/ALT allele match)"\n')
	metascr.write('cat alignmentresults.txt|awk -v FS="\t" -v OFS="\t" \'(FNR==1){print $0;for(i=1;i<=NF;++i){if($i=="aligned_beta"){fld_aligned_beta=i};if($i=="REFALT_alignment"){fld_REFALT_alignment=i};if($i=="Strand_alignment"){fld_Strand_alignment=i}}}((FNR>1)&&(fld_aligned_beta>0)&&(fld_REFALT_alignment>0)&&(fld_Strand_alignment>0)&&($fld_aligned_beta != ".")&&(($fld_Strand_alignment=="+")||($fld_Strand_alignment=="-"))&&(($fld_REFALT_alignment=="REFALT")||($fld_REFALT_alignment=="ALTREF"))){print $0}\' > alignmentresults_unambig.txt')
	metascr.close()
	os.chmod(metascriptname,0o755)
	eprint("At %s : Extracting and filtering metadata.\n"%(datetime.datetime.now()))
	os.system(". ./%s"%(metascriptname))

eprint("At %s : Reading summary stats / metadata match results.\n"%(datetime.datetime.now()))
metamatchfilename='./alignmentresults_unambig.txt'
metamatch=pd.read_csv(metamatchfilename, sep='\t', header=0, na_values=['.'])
eprint("\nRead data for %s markers in summary stats & genotype data"%(len(metamatch.index)))
print(metamatch.dtypes)

#removing duplicates from GWAavail

GWAavailAll=metamatch[['dataset_CHR','R9name','dataset_Marker','dataset_bp','aligned_beta','pval']] # keep only desired columns
GWAavailAll.columns=['CHR','R9name','SNP','BP','beta','P'] # change column name

GWAavailAll=GWAavailAll.drop_duplicates() # Drop duplicates if WHOLE row is the same ()
Duplicated=GWAavailAll.loc[GWAavailAll.SNP.duplicated(keep=False)].SNP

if len(Duplicated):
	eprint("Detected %s pairs of duplicated SNP IDs with different effect_alleles and effect sizes. This should not happen.\n"%(len(Duplicated)))
	eprint(GWAavailAll.loc[GWAavailAll.SNP.duplicated()])
	eprint("Duplicated IDs even after matching with metadata.\n")
	exit(1)

GWAavail=GWAavailAll[['CHR','R9name','SNP','BP','P']] # keep only columns wanted in clumping data file

# needs a test for duplication put back in here

eprint("Detected %s overlapping between summary stats and metadata."%(len(GWAavail.index)))
if len(GWAavail.index)<100000:
	eprint("Very few overlapping SNPs, this might be due to input files based on different genome build.\nPlease ensure you are using hg19 (Build 37).")
	raise RuntimeError('Not enough overlapping SNPs.')

######################################################################
#  Get list of marker blocks (also is the basis of 'sub-job table')  #
######################################################################

blocksoutdirSBR="SBPRS_calc/PRS_out"
blockslogdirSBR="SBPRS_calc/logfiles"
blocksoutdir="PRS_calc/PRS_out"
blockslogdir="PRS_calc/logfiles"

eprint("Determining counts of numbers of blocks per chromosome.\n")
try:
	blockdata=pd.read_csv(BlockMarkerCountFile,header=0,sep=' ',dtype={'chr':str})
except Exception as e:
	eprint("Error : could not open block/marker count file %s !\n"%(BlockMarkerCountFile))
	eprint("Exception was :\n%s"%(e))
	exit(1)
eprint("blockdata.dtypes :\n%s\n"%(blockdata.dtypes))

blockdatadict={(str(blockdata.chr[i])):(blockdata.loc[i]) for i in range(len(blockdata))}

AnyErrors=False
numsubjobs=0
subjobtable = pd.DataFrame(columns=['subjobnumber','chr','block','prefix'])
subjobtableSBR= pd.DataFrame(columns=['subjobnumber','chr','block','prefix'])
blocknames=[]
for chr in AnalysisChromosomes:
	nblocks=blockdatadict[chr].Nblocks_PLINK
	if nblocks == None:
		eprint("Error : file %s has no Nblocks_PLINK for chromosome %s !\n"%(BlockMarkerCountFile,chr))
		AnyErrors=True
	else:
		for blocknum in range(1,nblocks+1):
			numsubjobs += 1
			blockname="%s.%s"%(chr,blocknum)
			subjobtable.loc[len(subjobtable)] = [numsubjobs,chr, blockname, "%s/PRS_chr%s"%(blocksoutdir,blockname)]
			subjobtableSBR.loc[len(subjobtable)] = [numsubjobs,chr, blockname, "%s/PRS_chr%s"%(blocksoutdirSBR,blockname)]
			blocknames.append(blockname)
if AnyErrors:
	exit(1)
eprint("There will be/were %s sub jobs.\n"%(numsubjobs))

######################################################################
#               Submit SBayesR scoring using plink                   #
######################################################################

if not RunPRSCalculation:
	eprint("Skipping PRS calculation for individual marker blocks - by request.")
elif 'sbrA1' in metamatch.columns: #Only run below code if SBayesR was implicitly requested (e.g. ran in a previous run etc.)
	prefix='SBPRS_calc/'
	eprint('Creating SBayesR betas directory')
	subprocess.call(['mkdir','-p','SBayesRbetas'])
	for chr in [int(i) for i in AnalysisChromosomes]:
		tmpDF=metamatch.loc[metamatch.CHR==chr,['dataset_Marker','sbrA1','sbrBETA']]
		tmpDF.dropna(inplace=True)
		tmpDF.to_csv('SBayesRbetas/chr%s_betas.txt'%(chr),sep=' ',index=False,header=False)
	# Generate PRS calculation files
	eprint('Creating two new working directories: SBPRS_calc and SBPRS_calc/PRS_out')
	subprocess.call(['mkdir','-p',blocksoutdirSBR])
	subprocess.call(['mkdir','-p',blockslogdirSBR])
	subjobtablefile="SBPRS_calc/PRS_calc_subjobtable.txt"
	try:
		subjobtableSBR.to_csv(subjobtablefile,index=False,header=False,sep=' ')
	except exception as e:
		eprint("Error creating subjob table file !\nError message %s\n"%(e))
		exit(1)

	scriptname="SBPRS_calc/PRS_calc.sh"
	eprint("Creating PRS calculation job script (%s).\n"%(scriptname)) 
	try:
		currscript = open(scriptname, 'w')
		currscript.write("%s\n"%(JobScriptHeader))
		currscript.write("cd %s\n"%(normalwd))
		currscript.write("chr=$1\n") #Not sure why we do this?
		currscript.write("BlockName=$2\n") # and this
		currscript.write("prefix=$3\n") #and this
		currscript.write("subindex=$PBS_ARRAY_INDEX \n")
		currscript.write("echo \"subindex=${subindex}\"\n")
		# bash [[ -n var ]] True if a variable is not empty [[ -z var ]] True if var is empty. && executes next step if previous command didnt fail (return true)
		currscript.write("[[ -n ${subindex} ]] && VARIABLES=($(awk -v var=${subindex} 'NR==var' %s)) && chr=${VARIABLES[1]} && BlockName=${VARIABLES[2]} && prefix=${VARIABLES[3]}\n "%(subjobtablefile))
		currscript.write("[[ ((-z ${chr}) || (-z ${BlockName})) ]] && echo \"PBS job ID '${PBS_JOBID}' : No chr or BlockName in subjobtable\" && exit 1\n")
		currscript.write("(echo \"PBS job ID '${PBS_JOBID}' : chr=${chr} BlockName=${BlockName} prefix=${prefix}\"\n")
		currscript.write(" %s --dosage %s/PLINK_dosage/BlockPLINK_chr${BlockName}.dose.gz format=1 --fam %s/PLINK_dosage/GWAS.fam --keep %s/info/GWAS_Genotyped_FullIDs.txt --score SBayesRbetas/chr${chr}_betas.txt include-cnt sum --out ${prefix} >${prefix}.log1 2>&1\n "%(PLINKbinary,DatasetLocation,DatasetLocation,DatasetLocation))
		currscript.write(") > ${prefix}.log 2>&1\n")
		os.chmod(scriptname,0o755)
	except IOError:
		eprint("Error : could not write file %s !\n"%(scriptname))
		AnyErrors=True
	finally:
		currscript.close()
	if AnyErrors:
		exit(1)
	
	# Submit an array job for parallel computing
	
	finishedJobs={} #hash to keep track of submission and check if done
	eprint("At %s : Submitting one array job.\n"%(datetime.datetime.now()))
	
	subjobname="SBPRScalc"
	if JobNamePrefix != "":
		subjobname="%s_%s"%(JobNamePrefix,subjobname)
	Submit=subprocess.Popen('%s -l select=1:ncpus=1:mem=8GB -J 1-%s -l walltime=2:0:0 -r y -N "%s" -j oe -o SBPRS_calc/logfiles/ %s'%(qsubbinary,numsubjobs,subjobname,scriptname),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	outraw=Submit.communicate()
	outlines=outraw[0].decode("utf-8") + outraw[1].decode("utf-8")
	if outlines.find('qsub') >=0:
		eprint('qsub returned this message :\n%s\n'%(outlines))
		AnyErrors=True
	elif outlines=='':
		eprint('qsub returned nothing.\n')
		AnyErrors=True
	else :
		jobID=outlines.rstrip()
		eprint('.. job name is %s'%(jobID))
		finishedJobs[jobID]=0
	
	if AnyErrors:
		eprint('Error(s) in PBS job submission. Quiting and cancelling submitted job(s).\n')
		for i in finishedJobs:
				subprocess.Popen('qdel %s'%(i),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		exit(1)
	
	eprint("At %s : Job(s) submitted; waiting for job completion prior to continuing.\n"%(datetime.datetime.now()))
	
	# wait for job completion
	eprint("At %s : Finished submitting jobs; waiting for job completion.\n"%(datetime.datetime.now()))
	CheckCompletion(finishedJobs,120)

	logs=glob.glob('SBPRS_calc/logfiles/*.OU')
	# Checking for errors using simple grep
	AnyErrors=False
	for i in logs:
		Submits=subprocess.Popen("grep 'Error' %s"%(i),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		errlines=Submits.communicate()[0].decode("utf-8")
		Errors=[i for i in re.split(pattern='\n',string=errlines) if not i=='']
		if len(Errors) > 0:
			eprint("Log file %s contains these error(s) :\n%s\n"%(i,'\n'.join(Errors)))
			AnyErrors=True
	if AnyErrors:
		exit(1)
	else:
		eprint("No errors reported by HPC")
		
	eprint("Concctenating SBayesR PRS results")
	SBresultsFiles=glob.glob('%s/*.profile'%(blocksoutdirSBR))
	ProfDF=pd.read_csv(SBresultsFiles.pop(),header=0,sep=r"\s+")
	ProfDF=ProfDF.dropna()
	for profFile in SBresultsFiles:
		tmpProfDF=pd.read_csv(profFile,header=0,sep=r"\s+")
		tmpProfDF=tmpProfDF.dropna()
		if sum(tmpProfDF['IID'].values==ProfDF['IID'].values)!=len(ProfDF):
			raise RuntimeError('profile IDs are not in the same order plink output was likely wrong') 
		ProfDF.SCORESUM+=tmpProfDF.SCORESUM
		ProfDF.CNT+=tmpProfDF.CNT  
	if not args.output==None:
		outSB='%s/SBayR_%s.profile'%(normalwd,args.output)
	else:
		outSB='%s/SBayR_%s.profile'%(normalwd,PRSname)
	eprint("Finished concatenating SBayesR PRS results, writing output file to %s"%(outSB))
	ProfDF.to_csv(outSB,sep=" ",index=False)
	
	
	if args.standardized:
		print("Generating standardized PRS")
		stdFinalDF=ProfDF.copy()
		stdFinalDF['SCORESUM']=(stdFinalDF['SCORESUM']-stdFinalDF['SCORESUM'].mean())/stdFinalDF['SCORESUM'].std()
		print("Finished standardizing PRS")
		stdFinalDF.to_csv('%s/Std_SBayR_%s.profile'%(normalwd,PRSname),index=None,sep=' ',na_rep='NA')
	eprint("Finished Standardizing output is in the same directory as before")
	SBRprofDF=ProfDF.copy()
	StdSBRprofDF=stdFinalDF.copy()
	tmpProfDF=None
	ProfDF=None
	stdFinalDF=None
	


######################################################################
#                       Submit Clumping Jobs                         #
######################################################################

if not RunClumping :
	eprint("Skipping clumping stage, by request.")
else :
	# create clumping dir
	eprint("Create clumping working directory (necessary to submit jobs and keep things ordered).")
	subprocess.call(['mkdir','-p','clumping']) #creating a clumping directory to put scripts and clumped files

	# create clumping file
	eprint("Saving the \"available\" SNPs inside the clumping working space.\n")
	GWAavail.to_csv('clumping/GWA_SNPs_clumping_data.txt',header=True,index=None,sep=' ') #save the available snps to a table

	# generating clumping plink scripts

	# Change plink to use only EUR_IDs.txt
	eprint("Submitting one clump job per chromosome, to PBS.")
	finishedJobs={}
	logs=[]
	AnyErrors=False
	for chr in AnalysisChromosomes:
		eprint(".. Chromosome %s"%(chr))
		# Prefix for names of all files produced for in clumping, including the folder path
		prefix="clumping/clumped_chr%s"%(chr)
	
		# Files for LD calculation : 1000 Genomes Phase 3 European population
		refdir='/reference/genepi/public_reference_panels/1000G_20101123_v3_GIANT'
		refpopIDs=refdir+'/20101123_v3_GIANT_EUR_FullIDs.txt'
		if chr == 'X':
			reffilenamebase=refdir+'/derived_plinkformat/chrX.no.auto.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL'
		else:
			reffilenamebase=refdir+'/derived_plinkformat/chr%s.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL'%(chr)
	
		# Create and submit PBS job script
		scriptname=prefix+'.sh'
		with open(scriptname, 'w') as currscript:
			# Script header line (usually #!/bin/bash or similar)
			currscript.write("%s\n"%(JobScriptHeader))
			# Change to correct working directory
			currscript.write("cd %s\n"%(normalwd))
			# write the commands
			# .. marker data [CHR SNP(R9) BP P]
			currscript.write("(awk -v chr=%s '(FNR==1){print \"CHR\",\"SNP\",\"BP\",\"P\"}((FNR>1)&&($1==chr)){print $1,$2,$4,$5}' clumping/GWA_SNPs_clumping_data.txt > clumping/GWA_SNPs_clumping_chr%s.available\n"%(chr,chr))
			# .. Release 9 onward marker name (which 1KG data has been relabelled to match)
			currscript.write("awk '(FNR>1){print $2}' clumping/GWA_SNPs_clumping_chr%s.available > clumping/GWA_SNPs_clumping_chr%s.markerlist\n"%(chr,chr))
			# .. run PLINK
			#currscript.write("%s --bim %s.bim_R9names --bed %s.bed --fam %s.fam --keep %s --extract clumping/GWA_SNPs_clumping_chr%s.markerlist --clump clumping/GWA_SNPs_clumping_chr%s.available --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 10000 --out %s\n"%(PLINKbinary,reffilenamebase,reffilenamebase,reffilenamebase,refpopIDs,chr,chr,prefix))
			#Test with 1000G all which is what we used to use
			currscript.write("%s --bim %s.bim_R9names --bed %s.bed --fam %s.fam --extract clumping/GWA_SNPs_clumping_chr%s.markerlist --clump clumping/GWA_SNPs_clumping_chr%s.available --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 10000 --out %s\n"%(PLINKbinary,reffilenamebase,reffilenamebase,reffilenamebase,chr,chr,prefix))
			currscript.write(") > %s.log1 2>&1\n"%(prefix))
			currscript.close()
			os.chmod(scriptname,0o755)
			subjobname="clump%s"%(chr)
			if JobNamePrefix != "":
				subjobname="%s_%s"%(JobNamePrefix,subjobname)
			Submit=subprocess.Popen('%s -l select=1:ncpus=1:mem=8GB -l walltime=4:0:0 -r y -N "%s" -j oe -o %s.log_job %s'%(qsubbinary,subjobname,prefix,scriptname),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			outraw=Submit.communicate()
			outlines=outraw[0].decode("utf-8") + outraw[1].decode("utf-8")
			if outlines.find('qsub') >=0:
				eprint('qsub returned this message :\n%s\n'%(outlines))
				AnyErrors=True
			elif outlines=='':
				eprint('qsub returned nothing.\n')
				AnyErrors=True
			else :
				jobID=outlines.rstrip()
				eprint('.... job name is %s'%(jobID))
				logs.append(prefix+'.log');logs.append(prefix+'.log1');logs.append(prefix+'.log_job')
				finishedJobs[jobID]=0
	
	if AnyErrors:
		eprint('Error(s) in PBS job submission. Quiting and cancelling any submitted jobs.\n')
		for i in finishedJobs:
				subprocess.Popen('qdel %s'%(i),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		exit(1)
	
	eprint("At %s : Jobs submitted; waiting for job completion prior to continuing.\n"%(datetime.datetime.now()))
	
	CheckCompletion(finishedJobs,30)
		
	eprint("At %s : All clump jobs have finished running.\n"%(datetime.datetime.now()))
	
	# Checking for errors using simple grep
	AnyErrors=False
	for i in logs:
		Submits=subprocess.Popen("grep 'Error' %s"%(i),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		errlines=Submits.communicate()[0].decode("utf-8")
		Errors=[i for i in re.split(pattern='\n',string=errlines) if not i=='']
		if len(Errors) > 0:
			eprint("Log file %s contains these error(s) :\n%s\n"%(i,'\n'.join(Errors)))
			AnyErrors=True
	if AnyErrors:
		exit(1)
	else:
		eprint("No errors reported in clumping stage.\n")
	

#make pvalue and betas files
eprint("Making working directories to store Betas and Pvalues per chromosome (for the clumped data).")
subprocess.call(['mkdir','-p','Betas','Pvals'])

# Generating clumped betas and pvalue files
eprint("Generating the betas and pvalue files")

ClumpedChromosomes=[] #This is to store a list of chromosomes in case one does not have clumping results
for chr in AnalysisChromosomes:
	prefix='clumped_chr%s'%(chr)
	try:
		tmpClumped=pd.read_csv('clumping/%s.clumped'%(prefix),header=0,sep=r"\s+") #read the plink clumped files
	except FileNotFoundError:
		eprint("There were no clumping results for chr%s"%(chr))
		continue
	GWAtmp=metamatch.loc[metamatch['R9name'].isin(tmpClumped['SNP'])] # match the gwas to the clumped files
	GWAtmp=GWAtmp[['dataset_Marker','aligned_effect_allele','aligned_beta','pval']]
	GWAtmp=GWAtmp.dropna()
	
	# Write data file for beta & allele values [no header line]
	GWAtmp[['dataset_Marker','aligned_effect_allele','aligned_beta']].to_csv('Betas/%s_betas.txt'%(prefix),index=None,header=None,sep=' ') # save data to file [betas]

	# Write data file for pvalues [no header line]
	GWAtmp[['dataset_Marker','pval']].to_csv('Pvals/%s_pvals.txt'%(prefix),index=None,header=True,sep=' ')
	
	if Unweighted:
		# Generate and write unweighted values:
		GWAtmp['unweighted_beta']=np.sign(GWAtmp['aligned_beta'])
		GWAtmp[['dataset_Marker','aligned_effect_allele','unweighted_beta']].to_csv('Betas/%s_betas.txt.uw'%(prefix),index=None,header=None,sep=' ') # save data to file [betas]
	
	# record that there are markers for this chromosome
	ClumpedChromosomes.append(chr)

######################################################################
#        PLINK PRS scoring based on clumping + thresholding          #
######################################################################

if not RunPRSCalculation :
	eprint("Skipping PRS calculation for individual marker blocks - by request.")
else:
	if len(ClumpedChromosomes)>0:
		# Generate PRS calculation files
		eprint('Creating two new working directories: PRS_calc and PRS_calc/PRS_out')
		blocksoutdir="PRS_calc/PRS_out"
		subprocess.call(['mkdir','-p','PRS_calc',blocksoutdir])

		eprint("Determining counts of numbers of blocks per chromosome.\n")
		try:
			blockdata=pd.read_csv(BlockMarkerCountFile,header=0,sep=' ',dtype={'chr':str})
		except Exception as e:
			eprint("Error : could not open block/marker count file %s !\n"%(BlockMarkerCountFile))
			eprint("Exception was :\n%s"%(e))
			exit(1)

		eprint("blockdata.dtypes :\n%s\n"%(blockdata.dtypes))

		blockdatadict={(str(blockdata.chr[i])):(blockdata.loc[i]) for i in range(len(blockdata))}

		#Create table for array job parametrisation
		subjobtablefile="PRS_calc/PRS_calc_subjobtable.txt"
		eprint("Creating index file for PRS calculation script (%s).\n"%(subjobtablefile))
		try:
			subjobtable.to_csv(subjobtablefile,index=False,header=False,sep=' ')
		except exception as e:
			eprint("Error creating subjob table file !\nError message %s\n"%(e))
			exit(1)

		logs=[]
		# Write the script
		scriptname="PRS_calc/PRS_calc.sh"
		eprint("Creating PRS calculation job script (%s).\n"%(scriptname)) 
		try:
			currscript = open(scriptname, 'w')
			currscript.write("%s\n"%(JobScriptHeader))
			currscript.write("cd %s\n"%(normalwd))
			currscript.write("chr=$1\n")
			currscript.write("BlockName=$2\n")
			currscript.write("prefix=$3\n")
			currscript.write("subindex=$PBS_ARRAY_INDEX \n")
			currscript.write("echo \"subindex=${subindex}\"\n")
			currscript.write("[[ -n ${subindex} ]] && VARIABLES=($(awk -v var=${subindex} 'NR==var' %s)) && chr=${VARIABLES[1]} && BlockName=${VARIABLES[2]} && prefix=${VARIABLES[3]}\n "%(subjobtablefile))
			currscript.write("[[ ((-z ${chr}) || (-z ${BlockName})) ]] && echo \"PBS job ID '${PBS_JOBID}' : No chr or BlockName in subjobtable\" && exit 1\n")
			currscript.write("(echo \"PBS job ID '${PBS_JOBID}' : chr=${chr} BlockName=${BlockName} prefix=${prefix}\"\n")
			currscript.write(" %s --dosage %s/PLINK_dosage/BlockPLINK_chr${BlockName}.dose.gz format=1 --fam %s/PLINK_dosage/GWAS.fam --keep %s/info/GWAS_Genotyped_FullIDs.txt --score Betas/clumped_chr${chr}_betas.txt include-cnt sum --q-score-range %s Pvals/clumped_chr${chr}_pvals.txt 1 2 header --out ${prefix} >${prefix}.log1 2>&1\n"%(PLINKbinary,DatasetLocation,DatasetLocation,DatasetLocation,rangesfile))
			if Unweighted:
				currscript.write(" %s --dosage %s/PLINK_dosage/BlockPLINK_chr${BlockName}.dose.gz format=1 --fam %s/PLINK_dosage/GWAS.fam --keep %s/info/GWAS_Genotyped_FullIDs.txt --score Betas/clumped_chr${chr}_betas.txt.uw include-cnt sum --q-score-range %s Pvals/clumped_chr${chr}_pvals.txt 1 2 header --out ${prefix}_uw >${prefix}_uw.log1 2>&1\n"%(PLINKbinary,DatasetLocation,DatasetLocation,DatasetLocation,rangesfile))
			currscript.write(") > ${prefix}.log 2>&1\n")
			os.chmod(scriptname,0o755)
			logs.append(prefix+".log")
			logs.append(prefix+".log1")
		except IOError:
			eprint("Error : could not write file %s !\n"%(scriptname))
			AnyErrors=True
		finally:
			currscript.close()
		if AnyErrors:
			exit(1)
		
		# Submit an array job for parallel computing
		
		finishedJobs={} #hash to keep track of submission and check if done
		eprint("At %s : Submitting one array job.\n"%(datetime.datetime.now()))
		
		subjobname="PRScalc"
		if JobNamePrefix != "":
			subjobname="%s_%s"%(JobNamePrefix,subjobname)
		Submit=subprocess.Popen('%s -l select=1:ncpus=1:mem=8GB -J 1-%s -l walltime=2:0:0 -r y -N "%s" -j oe -o PRS_calc/PRS_calc.log_job %s'%(qsubbinary,numsubjobs,subjobname,scriptname),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		outraw=Submit.communicate()
		outlines=outraw[0].decode("utf-8") + outraw[1].decode("utf-8")
		if outlines.find('qsub') >=0:
			eprint('qsub returned this message :\n%s\n'%(outlines))
			AnyErrors=True
		elif outlines=='':
			eprint('qsub returned nothing.\n')
			AnyErrors=True
		else :
			jobID=outlines.rstrip()
			eprint('.. job name is %s'%(jobID))
			logs.append(prefix+'.log');logs.append(prefix+'.log1');logs.append(prefix+'.log_job')
			finishedJobs[jobID]=0
		
		if AnyErrors:
			eprint('Error(s) in PBS job submission. Quiting and cancelling submitted job(s).\n')
			for i in finishedJobs:
					subprocess.Popen('qdel %s'%(i),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			exit(1)
		
		eprint("At %s : Job(s) submitted; waiting for job completion prior to continuing.\n"%(datetime.datetime.now()))
		
		# wait for job completion
		eprint("At %s : Finished submitting jobs; waiting for job completion.\n"%(datetime.datetime.now()))
		CheckCompletion(finishedJobs,120)
	
		# Checking for errors using simple grep
		AnyErrors=False
		for i in logs:
			Submits=subprocess.Popen("grep 'Error' %s"%(i),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			errlines=Submits.communicate()[0].decode("utf-8")
			Errors=[i for i in re.split(pattern='\n',string=errlines) if not i=='']
			if len(Errors) > 0:
				eprint("Log file %s contains these error(s) :\n%s\n"%(i,'\n'.join(Errors)))
				AnyErrors=True
		if AnyErrors:
			exit(1)
		else:
			eprint("No errors reported by HPC")
	else:
		eprint("Clumping was not performed in this nor in a previous run. PLINK scoring based on clumping + thresholding cant be performed (if you didn't want this, try -runclumping True )")
		eprint("Unweighted scoring is also dependent on clumping")
		Unweighted=False
		MergeCTandSB=False


######################################################################
#               Concatenating PLINK profiles                         #
######################################################################

eprint("Will merge (sum) all of the plink profiles, separated by p-value threshold")
profileDir=normalwd+'/PRS_calc'+'/PRS_out/'
eprint("Profile directory is %s"%(profileDir))

finalDF=pd.DataFrame(columns=['FID', 'IID', 'PHENO'])

compPRSprefix='PRS_calc/Individual_PRS'
if not args.output==None:
	compPRSprefix='PRS_calc/%s'%(args.output)
if not CompilePRSBlocks:
	eprint("Skipping PRS compilation stage, by request.")
else:
	if not os.path.exists('PRS_calc/'):
		eprint("Clumping + Thresholding Scoring was not run in this or a previous attempt of this program. Will not concatenate CT PRS results")
	else:
		finishedJobs={}
		logs=[]
		AnyErrors=False
		for cutoff in cutoffnames:
			eprint("Concatenating profiles for cutoff %s using 'CompileProfileFiles' utility (starting at %s).\n"%(cutoff,datetime.datetime.now()))
			currProfileList=["%s/PRS_chr%s.%s.profile"%(blocksoutdir,blockname,cutoff) for blockname in blocknames]
			prefix='%s_%s'%(compPRSprefix,cutoff)
			scriptname="%s.sh"%(prefix)
			OutName="%s.profile"%(prefix)
			eprint("Creating PRS compilation job script (%s).\n"%(scriptname)) 
			eprint(".. output compiled scores will be in file %s\n"%(prefix))
			try:
				for f in [scriptname,prefix+'.log',prefix+'.log1',prefix+'.profile']:
					if os.path.exists(f):
						os.remove(f)
				currscript = open(scriptname, 'w')
				currscript.write("%s\n"%(JobScriptHeader))
				currscript.write("cd %s\n"%(normalwd))
				currscript.write("%s/Scripts/CompileProfileFiles.sh -fastmode PRS_calc/PRS_out/PRS_chr{%s}.%s.profile > %s\n"%(ReleaseLocation,','.join(blocknames),cutoff,OutName))
				if Unweighted:
					currscript.write("%s/Scripts/CompileProfileFiles.sh -fastmode PRS_calc/PRS_out/PRS_chr{%s}_uw.%s.profile > %s.unweighted\n"%(ReleaseLocation,','.join(blocknames),cutoff,OutName))
				currscript.close()
			except IOError:
				eprint("Error : could not write file %s !\n"%(scriptname))
				AnyErrors=True
			finally:
				os.chmod(scriptname,0o755)
				subjobname="compile%s"%(cutoff)
				if JobNamePrefix != "":
					subjobname="%s_%s"%(JobNamePrefix,subjobname)
				Submit=subprocess.Popen('%s -l select=1:ncpus=1:mem=6GB -l walltime=1:0:0 -r y -N "%s" -j oe -o %s.log_job %s'%(qsubbinary,subjobname,prefix,scriptname),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
				outraw=Submit.communicate()
				outlines=outraw[0].decode("utf-8") + outraw[1].decode("utf-8")
				if outlines.find('qsub') >=0:
					eprint('qsub returned this message :\n%s\n'%(outlines))
					AnyErrors=True
				elif outlines == '':
					eprint('qsub returned nothing.\n')
					AnyErrors=True
				else :
					jobID=outlines.rstrip()
					eprint('.. job name is %s; filename prefix is %s'%(jobID,prefix))
					logs.append(prefix+'.log'); logs.append(prefix+'.log1')
					finishedJobs[jobID]=0
		
		if AnyErrors :
				eprint('Error(s) in PBS job submission. Quiting and cancelling any submitted jobs.\n')
				for i in finishedJobs:
						subprocess.Popen('qdel %s'%(i),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
				exit(1)
			
		eprint("At %s : Jobs submitted; waiting for job completion prior to continuing.\n"%(datetime.datetime.now()))
			
		CheckCompletion(finishedJobs,30)
				
		eprint("At %s : All compiler jobs have finished running.\n"%(datetime.datetime.now()))
			
		# Checking for errors using simple grep
		AnyErrors=False
		for i in logs:
			Submits=subprocess.Popen("grep 'Error' %s"%(i),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			errlines=Submits.communicate()[0].decode("utf-8")
			Errors=[i for i in re.split(pattern='\n',string=errlines) if not i=='']
			if len(Errors) > 0:
				eprint("Log file %s contains these error(s) :\n%s\n"%(i,'\n'.join(Errors)))
				AnyErrors=True
		if AnyErrors:
			exit(1)
		else:
			eprint("No errors reported in compiler stage.\n")
		

######################################################################
#             Merging of cutoffs / Standardisation step              #
######################################################################
if CompilePRSBlocks:
	if not os.path.exists('PRS_calc/'):
		eprint("Will not compile CT PRS results)")
	else:
		eprint("Will now merge all scores in a single table.")
		first=1
		for currcutoff in cutoffnames:
			currinfile="%s_%s.profile"%(compPRSprefix,currcutoff)
			if not(os.path.exists(currinfile)):
				eprint("Error : cannot results for cutoff %s (file %s)"%(currcutoff,currinfile))
				exit(1)
			eprint("  Reading file %s"%(currinfile))
			x=pd.read_csv(currinfile,header=0,sep='\s+')
			x1=pd.DataFrame(data={'FID':x.FID, 'IID':x.IID, 'CNT_%s'%(currcutoff):x.CNT, currcutoff:x.SCORESUM,'Std_%s'%(currcutoff):(x.SCORESUM-x.SCORESUM.mean())/x.SCORESUM.std()})
			if first==1 :
				FinalDF=x1
			else :
				FinalDF=pd.merge(FinalDF, x1, how='outer', on=['FID','IID'])
			first=0

########### Unweighted compiling ###################
		if Unweighted:
			eprint("Will now include all the unweighted scores in the table.")
			for currcutoff in cutoffnames:
				currinfile="%s_%s.profile.unweighted"%(compPRSprefix,currcutoff)
				if not(os.path.exists(currinfile)):
					eprint("Error : cannot results for cutoff %s (file %s)"%(currcutoff,currinfile))
					exit(1)
				eprint("  Reading file %s"%(currinfile))
				x=pd.read_csv(currinfile,header=0,sep='\s+')
				x1=pd.DataFrame(data={'FID':x.FID, 'IID':x.IID, 'CNT_%s_uw'%(currcutoff):x.CNT, '%s_uw'%(currcutoff):x.SCORESUM,'Std_%s_uw'%(currcutoff):(x.SCORESUM-x.SCORESUM.mean())/x.SCORESUM.std()})
				FinalDF=pd.merge(FinalDF, x1, how='outer', on=['FID','IID'])
				
		if first==1 :
			eprint("Clumping + Thresholding attempted but nothing to compile - no cutoffs !")
			exit(1)

		OutName="ClumpThresh_%s.profile"%(PRSname)
		if MergeCTandSB:
			OutName="ClumpThresh_SBayR_%s.profile"%(PRSname)
			currinfile=outSB
			currcutoff='SBR'
			eprint("  Reading file %s"%(currinfile))
			x=pd.read_csv(currinfile,header=0,sep='\s+')
			x1=pd.DataFrame(data={'FID':x.FID, 'IID':x.IID, 'CNT_%s'%(currcutoff):x.CNT, currcutoff:x.SCORESUM,'Std_%s'%(currcutoff):(x.SCORESUM-x.SCORESUM.mean())/x.SCORESUM.std()})
			FinalDF=pd.merge(FinalDF, x1, how='outer', on=['FID','IID'])
		StdCols=['FID','IID']
		StdCols.extend([i for i in  FinalDF.columns if i.startswith('Std')])
		CNTCols=['FID','IID']
		CNTCols.extend([i for i in  FinalDF.columns if i.startswith('CNT')])
		nonStdCols=['FID','IID']
		nonStdCols.extend([i for i in  FinalDF.columns if not (i.startswith('Std') or i.startswith('CNT'))])
		FinalDF_Std=FinalDF.loc[:,StdCols]
		CNT_DF=FinalDF.loc[:,CNTCols]
		FinalDF=FinalDF.loc[:,nonStdCols]
		FinalDF.to_csv(OutName,index=None,sep=' ',na_rep='NA')
		CNT_DF.to_csv('CNT_%s'%(PRSname),index=None,sep=' ',na_rep='NA')
		if args.standardized:
			FinalDF_Std.to_csv('Std_'+OutName,index=None,sep=' ',na_rep='NA')


eprint("At %s : Finished running."%(datetime.datetime.now()))

#####################End of Script#######################
