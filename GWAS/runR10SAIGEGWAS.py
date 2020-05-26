#!/usr/bin/env python
#By Adrian Campos
# 13/11/19 Version 1.1
from sys import exit
try:
    import pandas as pd
    import re
    import argparse
    import subprocess
    import os
    import glob
    import time
    import sys
    import datetime
    import socket
    import numpy as np
    import scipy.stats as stats
    import random
    import matplotlib
    matplotlib.use('Agg')
    import seaborn as sns
    import matplotlib.pyplot as plt
    matplotlib.rcParams['agg.path.chunksize'] = 100000
    from scipy.stats import norm
    from scipy.stats import chi2
except ImportError:
    print('Python module does not seem to be loaded\nrun "module load python"')
    exit()

parser = argparse.ArgumentParser(description="By Adrian Campos for the genepi and QIMR HPC. In brief: Script to run a GWAS on a phenotype. The output is very close to the Ricopilli format. Three positional inputs are required (in the following order):")
parser.add_argument('phenofile',help="Phenotypes file: FID,IID Age Sex Covar1 ... Covarx ... qCovar1 ... qCovarX ... Phenotype\n order does not matter Missing values should be NA")
parser.add_argument('covdescript',help="Describe the variables on phenofile: ColumnName [covar|qcovar] do not decribe FID, IID nor the phenotype")
parser.add_argument('jobname',help="Job Name required to name the OUTPUT and jobs")
parser.add_argument('-outAFcasectrl','--IsOutputAFinCaseCtrl',default='FALSE',help="Whether to output case control allele frequencies either TRUE or FALSE")
parser.add_argument('-loco','--LOCO',default='FALSE',help="Whether to perform a LOCO either TRUE or FALSE")
parser.add_argument('-chroms','--chromosomes',nargs='+',default=list(range(1,23)),help="Which cromosomes to use (1-22 by default) space separated e.g. -chroms 1 2 15 22 ")
parser.add_argument('-maf','--minMAF',default=0.01,help="Minimum MAF threshold to include variant")
parser.add_argument('-numPCs','--NUMBERofPCs',default=10,help="How many PCs to include as covariates use 0 to add none")
parser.add_argument('-nosub','--DoNotSubmit',help="Do not submit jobs (only generate PBS scripts)",action='store_true')
parser.add_argument('-chip',help="Which chip? All, GSA, Hap etc; default All",default='All')
parser.add_argument('-duplicates','--DealDuplicates',help="automatically deal with dups valid is nothing (by default), remove or average",default = None)


parser._actions[0].help='Print this help message and exit. This script is a standardized way of performing a GWAS using PLINK and the release 10 R10. It is written in python and it is recommended to run it interactively (so basically submit it as a job that will submit all other jobs).\nBecause this job will wait for all the child jobs it WILL need a lot of walltime (give it 24 hours)'
args = parser.parse_args()
print("AGDS_saige_GWAS script started running at %s"%(datetime.datetime.now()))

######################### TO TEST INTERACTIVELY #####################
'''
class z():
     def __init__(self):
         self.jobname='suiatt'
         self.phenofile='SuiattAll.pheno'
         self.covdescript='covdescript'
         self.DoNotSubmit='False'
         self.NUMBERofPCs=0
         self.minMAF=0.01
         self.LOCO='FALSE'
         self.IsOutputAFinCaseCtrl='FALSE'
         self.chromosomes=range(1,23)
         self.DealDuplicates=None
args=z()

'''


#################Pre needed functions for HPC submision ##############################
normalwd=os.getcwd()
FNULL = open(os.devnull, 'w') # to redirect HPC output

def CheckCompletion(jobDic,timew=300):
    '''Check completion function for submitted jobs to the hpc cluster
    input a job dictionary of job IDs and zeros (means they are running), and a predetermined wait time between cluster query'''
    ids=[i for i in jobDic.keys()] #get all job ids
    submstr=str("for job in %s; do qstat ${job}; done"%(' '.join(ids))) #submit query for all jobs
    while sum(jobDic.values())<len(jobDic):
        time.sleep(timew)
        #submit the query
        submits=subprocess.Popen(submstr,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        #get the output of the query
        out=''.join([i.decode("utf-8") for i in submits.communicate()])
        #split lines of output
        OutsStats=[i for i in re.split(pattern='\n',string=out) if not i=='']
        #get ids of finished jobs
        finished_outStats=[i for i in OutsStats if re.search('Job has finished',i)]
        #update finished hash table
        for id in [re.search('qstat: (\d+).hpc',i).group(1) for i in finished_outStats]:
            jobDic[id]=1
    return(1) #return something just in case you need to know, this wont exit until they are finished OMG I've used a while for the first time!
    
    
def SubmitAllJobs(ext='.PBS'):
    ''' Submits all files with an extension in the current working directory'''
    import subprocess
    import re
    Submits=subprocess.Popen('for script in *%s; do qsub ${script}; done'%(ext),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    outlines=Submits.communicate()[0].decode("utf-8")
    jobs={i.split('.')[0]:0 for i in re.split(pattern='\n',string=outlines) if not i==''}
    return(jobs)
    
    
def CheckErrors(ext='*.e*'):
    '''
        Error checker, basically greps Error in all of the files *.e* of the currwd
    '''
    import subprocess
    import re
    Submits=subprocess.Popen("grep 'Error' %s"%(ext),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    outlines=Submits.communicate()[0].decode("utf-8")
    Errors=[i for i in re.split(pattern='\n',string=outlines) if not i=='']
    return(Errors)
    
    
### Print and check for errors.
print("By Adrian Campos @ QIMR Berghofer\nplease let me know of any problems errors or modifications needed\nadrian.campos@qimrberghofer.edu.au")
print("This script is provided as is with no warranty and under under a CC license ...")
if socket.gethostname()=='hpcpbs01.adqimr.ad.lan':
    print("This script opens a lot of files to memory ~80Gb, run it interactively or with qsub!")
    raise RuntimeError('You should run this script interactively or via qsub, not in login nodes')
    


################ Read pheno, args include covars etc
################### Read Arguments ############################################
Phenofile=args.phenofile  
print('Phenofile is %s'%(Phenofile))
jobName=args.jobname    
print('Job name is %s'%(jobName))
covdescriptfile=args.covdescript
print('Cov descript name is %s'%(covdescriptfile))
chrs=args.chromosomes
print('Chromosomes to be studied are  %s'%(' '.join([str(i) for i in chrs])))

numPCs=int(args.NUMBERofPCs)

###################File preprocessing ####################
phenodf=pd.read_table(Phenofile,header=0,sep=r"\s+")
phenodf.index=phenodf.IID.astype(str)


#duplicates
if args.DealDuplicates == 'remove':
    print("Removing duplicated elements from the phenotype file")
    phenodf=phenodf.loc[~phenodf.index.duplicated(),:]
elif args.DealDuplicates == 'average':
    print("Averaging duplicated elements from the phenotype file this will fail with non numerics")
    dups=list(set(phenodf.index[phenodf.index.duplicated(keep=False)]))
    for ix in dups:
        currFID=phenodf.loc[ix,:].FID.values[0]
        tmp=phenodf.loc[ix,[i for i in phenodf if i not in ['FID']]].mean()
        tmp['FID']=currFID
        tmp.name=ix
        phenodf.drop(ix,inplace=True)
        phenodf=phenodf.append(tmp)

#Open covar names
if args.covdescript!="None":
    covarnames=pd.read_table(covdescriptfile,header=None,sep=r"\s+",index_col=0)
    covarDict= covarnames.to_dict()[1]
    print ('Opened all of the input files')
#match them

print("Creating and changing to GWAS_calc working directory (to submit from there PLINK scripts)")




# Generate GWAS calculation files
print('Creating two new working directories')
subprocess.call(['mkdir','GWAS_calc'])
subprocess.call(['mkdir','GWAS_calc/GWAS_out'])

os.chdir(normalwd+'/GWAS_calc') # move there
## VARIABLES TO CREATE INPUT FILES:
phenodf['IID']=phenodf.index
if args.covdescript!="None":
    covariates=[i for i in covarDict if covarDict[i]=='covar']
    qcovariates=[i for i in covarDict if covarDict[i]=='qcovar']
else:
    covariates=[]
    qcovariates=[]
orderedCols=['FID','IID']
phenotype=[i for i in phenodf.columns if i not in covariates and i not in orderedCols and i not in qcovariates]
if len(phenotype)>1:
    print("Detected more than one phenotype, will use only:%s"%(phenotype[0]))
phenotype=phenotype[0]


tmpdf=phenodf[['FID','IID',phenotype]+covariates+qcovariates] # create phenotype file
tmpdf=tmpdf.dropna()
levels=set(tmpdf.loc[:,phenotype])
if (len(levels) == 2):
    binpheno=True
    traitTypeLine="traitType = 'binary'"
    print ("the phenotype is binary")
else:
    binpheno=False
    traitTypeLine="traitType = 'quantitative'"    
    print ("the phenotype is taken to be quantitative")
        
    
################ add  PCs ###############
print("Adding %s number of PCs max is 20"%(numPCs))
if numPCs:
    PCdata=pd.read_csv('/reference/genepi/GWAS_release/Release10/Release10_Observed/AncestryChecks/GWAS_PrincipalComponentsScores.txt',
                        sep = '\s+',index_col='INDID')
    for i in range(numPCs):
        currCol='PC'+str(i+1)
        qcovariates.append(currCol)
        tmpdf[currCol]=PCdata.loc[tmpdf.IID.astype(str).values,currCol]
    
    
    
tmpdf['IID']=tmpdf['FID'].astype(str)+'_'+tmpdf['IID'].astype(str)
tmpdf.drop('FID', inplace=True,axis=1)
tmpdf.to_csv(jobName+'.pheno',index=False,sep="\t")

################Save file ######################
if args.covdescript!="None" or numPCs:
    covarColLine='covarColList = c(\''+'\',\''.join(covariates+qcovariates)+'\')'
else:
    covarColLine='covarColList =NULL'


shebang=r'#!/bin/bash' #header to add to each file consider env bash

if args.chip=="GSA":
    plinkFileLine="/reference/genepi/GWAS_release/Release10/Scripts/SAIGE/required_files/GWAS_Release10_draft4_02102019"
else:
    plinkFileLine="/reference/genepi/GWAS_release/Release10/Scripts/SAIGE/required_files/GWAS_Release10_consensusset_SAIGE"
print("Creating the STEP1 submission script")
with open('Step1%s.PBS1'%(jobName), 'w') as currscript:
    currscript.write(shebang+'\n') #write shebang
    #write PBS headers
    currscript.write("#PBS -N SAIGEstep1\n#PBS -r n\n#PBS -l mem=40GB,walltime=30:00:00,ncpus=5\nmodule load R/3.5.1\ncd $PBS_O_WORKDIR\n")
    currscript.write("R -e \"library(SAIGE);fitNULLGLMM(plinkFile='%s',phenoFile='%s',phenoCol = '%s',%s,%s,sampleIDColinphenoFile='IID',outputPrefix = 'SAIGE_Step1_%s', nThreads = 5)\""%(plinkFileLine,jobName+'.pheno',phenotype,traitTypeLine,covarColLine,jobName))
# Submit step1

if args.DoNotSubmit==True:
    print("Do not submit flag present, will generate step2 scripts and exit")
else:
    finishedJobs={} #hash to keep track of submission and check if done
    print("Submitting the first step")
    print("Started submitting jobs %s"%(datetime.datetime.now()))
    finishedJobs=SubmitAllJobs(ext='.PBS1') # hash table with completion info
    print("Finished submitting jobs %s\n Waiting for job completion"%(datetime.datetime.now()))

    # wait for job completion
    CheckCompletion(finishedJobs,120)
    print("All GWAS jobs finished running %s +- 2 min"%(datetime.datetime.now()))


    #check for job errors
    anyErrors=CheckErrors()
    if len(anyErrors):
        print("The following errors were reported by HPC:\n")
        for currError in anyErrors:
            print(currError)
    else:
        print("No errors reported by hpc :)")

# Creating Step2

print("Creating the STEP2 submission script")
for chrom in chrs: #for each chromosome that we have clumping results for
    #obtain all blocknames
    blockNames=[i for i in glob.glob('/reference/genepi/GWAS_release/Release10/Release10_HRCr1.1/VCF/BlockVCF_chr%s.*.gz'%(chrom)) if not re.search('poly',i)]
    for block in blockNames: #for each block
        numblock=str(re.search('BlockVCF_chr\d+.(\d+).dose.vcf.gz',block).group(1)) #obtain current block number
        #open the file, the with closes it even if there are problems
        vcfFile='/reference/genepi/GWAS_release/Release10/Release10_HRCr1.1/VCF/BlockVCF_chr%s.%s.dose.vcf.gz'%(chrom,numblock)
        vcfTbix=vcfFile+'.tbi'
        sampleFile='/reference/genepi/GWAS_release/Release10/Release10_HRCr1.1/VCF/VCF_chr%s_VCFIDsequence.txt'%(chrom)
        GMMATmodelFile='SAIGE_Step1_%s.rda'%(jobName)
        varianceRatioFile='SAIGE_Step1_%s.varianceRatio.txt'%(jobName)
        with open('chrom%s_block%s'%(chrom,numblock)+'.PBS', 'w') as currscript:
            currscript.write(shebang+'\n')
            currscript.write("#PBS -N SAIGE\n#PBS -r n\n#PBS -l ncpus=1,mem=10GB,walltime=20:00:00\nmodule load R/3.5.1\ncd $PBS_O_WORKDIR\n")
            currscript.write("R -e \"library(SAIGE);SPAGMMATtest(vcfFile='%s',vcfFileIndex='%s',vcfField='DS',chrom='%s',sampleFile='%s',GMMATmodelFile='%s',varianceRatioFile='%s',SAIGEOutputFile='GWAS_out/%s_sumstats_chr%s_block%s.saige',minMAF = %s,IsOutputAFinCaseCtrl = %s,LOCO = %s)\""%(vcfFile,vcfTbix,chrom,sampleFile,GMMATmodelFile,varianceRatioFile,jobName,chrom,numblock,args.minMAF,args.IsOutputAFinCaseCtrl,args.LOCO))
            
if args.DoNotSubmit==True:
    print("Do not submit flag present, exiting now")
    exit()
finishedJobs={} #hash to keep track of submission and check if done
print("Started submitting jobs %s"%(datetime.datetime.now()))
finishedJobs=SubmitAllJobs(ext='.PBS') # hash table with completion info
print("Finished submitting jobs %s\n Waiting for job completion"%(datetime.datetime.now()))

# wait for job completion
CheckCompletion(finishedJobs,120)
print("All GWAS jobs finished running %s +- 2 min"%(datetime.datetime.now()))

#check for job errors
anyErrors=CheckErrors()
if len(anyErrors):
    print("The following errors were reported by HPC:\n")
    for currError in anyErrors:
        print(currError)
else:
    print("No errors reported by hpc :)")


print("Change to main working directory")
os.chdir(normalwd)

os.chdir('GWAS_calc/GWAS_out')
dosageFiles=glob.glob('*.saige')


print("Compiling dosage results this may take a while %s"%(datetime.datetime.now()))
finalDF=pd.concat((pd.read_csv(f,header=0,index_col='SNPID',sep='\s+') for f in dosageFiles))
print("Finished compiling dosage at %s"%(datetime.datetime.now()))


print("Matching meta data (rsNumber and imputation status) %s"%(datetime.datetime.now()))
meta_data=pd.read_table('/reference/genepi/GWAS_release/Release10/Release10_HRCr1.1/info/Metadata_allchr.txt',usecols=['Marker','Marker_dbSNP_1','Genotyped_Hap','Genotyped_GSA','Genotyped_Omni','Rsq_rederived'],sep="\s+",na_values='.')
meta_data['Genotyped']=(meta_data.loc[:,['Genotyped_Hap','Genotyped_GSA','Genotyped_Omni']]==1).any()*1
meta_data.index=meta_data.Marker
meta_data['Genotyped']=meta_data.Genotyped=='1'
meta_data.Genotyped=meta_data.Genotyped*1
finalDF['SNPrs']=meta_data.loc[finalDF.index,'Marker_dbSNP_1']
finalDF['ngt']=meta_data.loc[finalDF.index,'Genotyped']
finalDF['INFO']=meta_data.loc[finalDF.index,'Rsq_rederived']



print("Saving the compiled results to %s" %(normalwd))
os.chdir(normalwd)
finalDF.to_csv('%s_GWASsumstats.dat'%(jobName),sep='\t',na_rep='NA')

    
print("Finished running %s"%(datetime.datetime.now()))

    
