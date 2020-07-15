#!/usr/bin/env python
#By Adrian Campos
# 13/11/19 Version 1.1
#Fixed bug with reading maf and info score cutoffs as strings
#New PCs obtained from draft 3
#Updated to new Release10 GSA version. added qq and manhattan plots v2
# Matches metadata to almost final ricopilli format future version will allow different output formats
# PLINK GWAS wrapper for the Release10 (particularly for the AGDS and GWAS of non-related individuals using plink)
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

####################Arg parsing #########################
parser = argparse.ArgumentParser(description="By Adrian Campos for the genepi and QIMR HPC. In brief: Script to run a GWAS on a phenotype. The output is very close to the Ricopilli format. Three positional inputs are required (in the following order):")
parser.add_argument('phenofile',help="Phenotypes file: FID,IID Age Sex Covar1 ... Covarx ... qCovar1 ... qCovarX ... Phenotype\n order does not matter Missing values should be NA")
parser.add_argument('covdescript',help="Describe the variables on phenofile: ColumnName [covar|qcovar] do not decribe FID, IID nor the phenotype")
parser.add_argument('jobname',help="Job Name required to name the OUTPUT and jobs")
parser.add_argument('-remove','--Remove',help="Specify a list of individuals to exclude from the analysis (i.e. related and PC outliers)")
parser.add_argument('-keep','--Keep',help="Specify a list of individuals to exclude from the analysis (i.e. related and PC outliers)")
parser.add_argument('-exclude','--Exclude',help="Specify a list of variants to exclude from the analysis (e.g. failed QC)")
parser.add_argument('-duplicates','--DealDuplicates',help="automatically deal with dups valid is nothing (by default), remove or average",default = None)
parser.add_argument('-chroms','--chromosomes',nargs='+',default=list(range(1,23)),help="Which cromosomes to use (1-22 by default) space separated e.g. -chroms 1 2 15 22 ")
parser.add_argument('-numPCs','--NUMBERofPCs',default=10,help="How many PCs to include as covariates use 0 to add none")
parser.add_argument('-noimput','--DoNotAddImput',help="Do not include imputation covariates automatically (use if they are already in your phenofile)",action = 'store_true')
parser.add_argument('-nosub','--DoNotSubmit',help="Do not submit jobs (only generate PBS scripts)",action='store_true')
parser.add_argument('-qcMAF','--PerformQCMAF',help="Whether to perform QC removing alleles with a MAF < cutoff (specify cutoff), remove nans",type=float,default=0.01)
parser.add_argument('-qcINFO','--infoS',help="Whether to remove alleles with low info score",type=float,default=0.8)
parser.add_argument('-plot','--Plot',help="Do a QQ and Manhattan plot of QC results (check QC defaults above)",action='store_true')



parser._actions[0].help='Print this help message and exit. This script is a standardized way of performing a GWAS using PLINK and the release 10 R10. It is written in python and it is recommended to run it interactively (so basically submit it as a job that will submit all other jobs).\nBecause this job will wait for all the child jobs it WILL need a lot of walltime (give it 24 hours)'
args = parser.parse_args()
print("AGDSGWAS script started running at %s"%(datetime.datetime.now()))

######################### TO TEST INTERACTIVELY #####################
'''
class z():
     def __init__(self):
         self.chromosomes=list(range(1,23))
         self.Keep=None
         self.jobname='test'
         self.phenofile='QIMR_teenacnewithmild.dat'
         self.covdescript='None'
         self.log=None
         self.DoNotSubmit=True
         self.DealDuplicates=None
         self.Remove=None
         self.DoNotAddImput=False
         self.Exclude=None
         self.PerformQCMAF=0.05
         self.infoS=0.6
         self.Plot=True
         self.NUMBERofPCs=10
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
    
binpheno=False




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
phenodf.set_index('IID', inplace=True)

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
phenotypes=[i for i in phenodf.columns if i not in covariates and i not in orderedCols and i not in qcovariates]
if len(phenotypes)>1:
    print("Detected more than one phenotype, will use All of them")

tmpdf=phenodf[['FID','IID']+phenotypes] # create phenotype file
tmpdf=tmpdf.dropna()
for phenotype in phenotypes:
    levels=set(tmpdf.loc[:,phenotype])
    if (len(levels) == 2):
        binpheno=True
        if 0 in levels and 1 in levels:
            print ("the phenotype seems binary, it will be recoded to plink format (1 case 2 control)")
            tmpdf.loc[:,phenotype]=tmpdf.loc[:,phenotype].map({1:2,0:1})
    if len(levels) < 2:
        print ("there was only one level for the phenotype %s, this wont work will be dropped"%(phenotype))
        tmpdf.drop(phenotype,axis=1)
    if (len(levels) > 2) and (0 in levels):
        print("Note that 0 will be considered as missing as per plink specifications")
tmpdf.to_csv(jobName+'.pheno',index=False,sep="\t")

if args.covdescript!="None":
    print("Creating GWAS input files: covar")
    ######### Create covar file (only one needed) ##################
    cols=['FID','IID']
    cols.extend(covariates)
    cols.extend(qcovariates)
    tmpdf=phenodf[cols]

    for cov in covariates:
        levels=set(tmpdf.loc[tmpdf[cov].notna(),cov])
        print ("The levels for covariate %s are %s "%(cov,levels))
        if (len(levels) == 2) and 0 in levels and 1 in levels:
            print ("this covariate seems binary, it will be recoded to plink format (1 case 2 control)")
            tmpdf.loc[:,cov]=tmpdf.loc[:,cov].map({1:2,0:1})
        if len(levels) < 2:
            print ("there was only one level for this covariate, it will be dropped off the model")
            tmpdf.drop(cov,axis=1,inplace=True)
        if (len(levels) > 2) and (0 in levels):
            print("Note that 0 will be considered as missing as per plink specifications")
elif numPCs:
    cols=['FID','IID']
    tmpdf=phenodf[cols]
################ Define plink modifiers ##########################

removeLine=''
if args.Remove:
    removeLine='--remove %s'%(args.Remove)
        
excludeLine=''
if args.Exclude:
    excludeLine='--exclude %s'%(args.Exclude)
    
keepLine=''
if args.Keep:
    keepLine='--keep %s'%(args.Keep)
    
################ add  PCs ###############
print("Adding %s number of PCs max is 20"%(numPCs))
if numPCs:
    PCdata=pd.read_csv('/reference/genepi/GWAS_release/Release10/Release10_Observed/AncestryChecks/GWAS_PrincipalComponentsScores.txt',
                        sep = '\s+',index_col='INDID')
    for i in range(numPCs):
        currCol='PC'+str(i+1)
        qcovariates.append(currCol)
        tmpdf[currCol]=PCdata.loc[tmpdf.IID,currCol]
        
############# add imputation covariate ########
if not args.DoNotAddImput:
    print("Adding imputation covariates")
    imputcovfile='/reference/genepi/GWAS_release/Release10/Release10_HRCr1.1/info/GWAS_ImputationRunCovariates.ped'
    imputcovdf=pd.read_table(imputcovfile,header=None,sep=r"\s+",usecols=[1,5,6])
    imputcovdf.index=imputcovdf[1]
    imputcovdf.drop(1,axis=1,inplace=True)
    imputcovdf=imputcovdf.loc[tmpdf.IID,:]
    imputcovdf.columns=['Imput1','Imput2']
    imputcovdf['Imput1']=imputcovdf.Imput1.map({1:2,0:1})
    imputcovdf['Imput2']=imputcovdf.Imput2.map({1:2,0:1})
    for column in imputcovdf:
        if len(imputcovdf[column].value_counts())>1:
            tmpdf[column]=imputcovdf.loc[tmpdf.IID,column]
        else:
            print("Excluded %s due to lack of variance"%(column))

################Save file ######################
if args.covdescript!="None" or numPCs:
    tmpdf=tmpdf.dropna()
    tmpdf.to_csv(jobName+'.covariates',index=False,sep=" ")
    covarFileName=jobName+'.covariates' #--covarFile
    covarFileLine='--covar %s'%(covarFileName)
else:
    covarFileLine=''

shebang=r'#!/bin/bash' #header to add to each file consider env bash

print("Creating the GWAS PBS submission scripts, using all blocks (except poly) in:\n \"/reference/genepi/GWAS_release/Release10/Release10_HRCr1.1_GSA/PLINK_dosage/\"")
for chrom in chrs: #for each chromosome that we have clumping results for
    #obtain all blocknames
    blockNames=[i for i in glob.glob('/reference/genepi/GWAS_release/Release10/Release10_HRCr1.1_GSA/PLINK_dosage/BlockPLINK_chr%s.*.gz'%(chrom)) if not re.search('poly',i)]
    for block in blockNames: #for each block
        numblock=str(re.search('BlockPLINK_chr\d+.(\d+).dose.gz',block).group(1)) #obtain current block number
        #open the file, the with closes it even if there are problems
        with open('chrom%s_block%s'%(chrom,numblock)+'.PBS', 'w') as currscript:
            currscript.write(shebang+'\n') #write shebang
            #write PBS headers
            currscript.write("#PBS -N GWASchr%sB%s\n#PBS -r n\n#PBS -l mem=12GB,walltime=2:00:00\nmodule load plink/1.90b6.8\ncd $PBS_O_WORKDIR\n"%(chrom,numblock))
            #write plink dosage inputs
            currscript.write("plink --dosage %s format=1 case-control-freqs --fam /reference/genepi/GWAS_release/Release10/Release10_HRCr1.1_GSA/PLINK_dosage/GWAS.fam --pheno %s.pheno --all-pheno %s %s %s %s --out GWAS_out/%s_sumstats_chr%s_block%s --memory 8000 --threads 1"%(block,jobName,covarFileLine,removeLine,excludeLine,keepLine,jobName,chrom,numblock))
        
# Submit them all for parallel computing

# submit jobs
if args.DoNotSubmit==True:
    print("Do not submit flag present, exiting now")
    exit()
finishedJobs={} #hash to keep track of submission and check if done
print("Submitting one job per block this may take a while")
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


for phenotype in phenotypes:
    os.chdir('GWAS_calc/GWAS_out')
    dosageFiles=glob.glob('*.dosage')


    print("Compiling dosage results this may take a while %s"%(datetime.datetime.now()))
    finalDF=pd.concat((pd.read_csv(f,header=0,index_col='SNP',sep='\s+') for f in dosageFiles))
    print("Finished compiling dosage at %s"%(datetime.datetime.now()))


    print("Matching meta data (rsNumber and imputation status) %s"%(datetime.datetime.now()))
    meta_data=pd.read_table('/reference/genepi/GWAS_release/Release10/Release10_HRCr1.1_GSA/info/Metadata_allchr.txt',usecols=['Marker','Marker_dbSNP_1','Genotyped'],sep="\s+")
    meta_data.index=meta_data.Marker
    meta_data['Genotyped']=meta_data.Genotyped=='1'
    meta_data.Genotyped=meta_data.Genotyped*1
    finalDF['SNPrs']=meta_data.loc[finalDF.index,'Marker_dbSNP_1']
    finalDF['ngt']=meta_data.loc[finalDF.index,'Genotyped']
    finalDF['CHR']=[int(re.split(":",i)[0]) for i in finalDF.index]
    finalDF['BP']=[int(re.split(":",i)[1]) for i in finalDF.index]

    if binpheno:
        print("Obtaining number of cases and controls")
        logFile=dosageFiles[0].replace('.assoc.dosage','.log')
        with open(logFile,'r') as logf:
            for line in logf:
                if re.match('Among remaining phenotypes,',line):
                    numbers=re.findall(r'\d+', line)
                    cases=numbers[0]
                    controls=numbers[1]
                    cases=int(cases)
                    controls=int(controls)
                    break
        newNames={'FRQ_A':'FRQ_A_'+str(cases),'FRQ_U':'FRQ_U_'+str(controls)}
    else:
        newNames={}
        cases=1
        controls=1

    finalDF.columns=[newNames.get(i,i) for i in finalDF.columns]
    print("Saving the compiled results to %s" %(normalwd))
    os.chdir(normalwd)
    finalDF.to_csv('%s_%s_GWASsumstats.dat'%(jobName,phenotype),sep='\t',na_rep='NA')

    ## QC steps 

    #AF 
    if binpheno:
        AF=(finalDF['FRQ_A_'+str(cases)]*cases+finalDF['FRQ_U_'+str(controls)]*controls)/(cases+controls)
    else:
        AF=(finalDF['FRQ_A']*cases+finalDF['FRQ_U']*controls)/(cases+controls)
        
    print('QC step removing sumstats with MAF < %s'%(args.PerformQCMAF))
    finalDF=finalDF.loc[(AF>args.PerformQCMAF) & (AF<(1-args.PerformQCMAF)),:]

    #Info
    print('QC step removing sumstats with info < %s'%(args.infoS))
    finalDF=finalDF.loc[finalDF.INFO>=args.infoS,:]
    finalDF.to_csv('%s_%s_GWASsumstats_QCed.dat'%(jobName,phenotype),sep='\t',na_rep='NA')

    if args.Plot:
    ## qqplot based on script provided by Enda Byrne
        finalDF.dropna(inplace=True)
        print('Now plotting qqplot and manhattan plot')
        z=norm.ppf(finalDF['P']/2)
        lambGC = round(np.median(z**2)/chi2.ppf(0.5,df=1),3)
        p = 2*norm.cdf(-abs(z))
        p.sort()
        expected = list(range(len(p)))
        lobs = np.array(-np.log10(p))
        lexp = np.array(-np.log10(np.array(expected) / (len(expected)+1)))
        # plots all points with p < 0.05
        fig,ax = plt.subplots(1)
        p_sig = [i for i in p if i<0.05]
        ax.plot(lexp, lobs, 'd', color='steelblue',alpha=0.75,label='Observed')
        xmin,xmax=ax.get_xlim()
        ax.plot([0,xmax+0.2], [0,xmax+0.2],linestyle='dashed', color='black',label='expected')
        ax.set_xlim(xmin,xmax)
        ax.set_xlabel('-log10(expected P)')
        ax.set_ylabel('-log10(observed P)')
        ax.text(0.2,6,'lambda = %s'%(str(lambGC)))
        fig.savefig('%s_%s_QCed.qqplot.png'%(jobName,phenotype))
        
    ## Manhattan plot
        GWAsumstats1=finalDF
        GWAsumstats1.sort_values(by=['CHR','BP'],inplace=True)
        # calculate new positions for plot:
        dftmp=GWAsumstats1
        dftmp=dftmp.loc[:,['SNPrs','CHR','BP']]
        dftmp['SNP']=dftmp.index
        dftmp.drop_duplicates(inplace=True)

        if sum(dftmp.SNPrs.duplicated()):
            print("There are duplicated variants in sumstats with different positions. One will randomly be removed")
            dftmp=dftmp.loc[~dftmp.SNPrs.duplicated(),:]

        pastlen=len(dftmp)
        dftmp.dropna(inplace=True)
        GWAsumstats1.dropna(inplace=True)
        if len(dftmp)!=pastlen:
            print("There were NAs in the sumstats, they were removed automatically but you might want to double check")

        chrs=set(args.chromosomes) #make sure its sorted
        pastChr=chrs.pop()
        maxBP=max(dftmp.loc[dftmp.CHR==pastChr,'BP'])
        for chrom in chrs:
            tmpBP=dftmp.loc[dftmp.CHR==chrom,'BP']
            tmpBP=tmpBP+maxBP
            dftmp.loc[tmpBP.index,'BP']=tmpBP
            maxBP=max(dftmp.loc[dftmp.CHR==chrom,'BP'])

        dftmp.set_index(dftmp.SNPrs,inplace=True)


        ##Plot
        fig = plt.figure() 
        topAx=fig.add_subplot(111)
        fig.set_size_inches(25,10)
        #Manhattan 1
        colors1=['#a92e4a' if i%2 else '#782135' for i in GWAsumstats1.CHR.values]
        topAx.scatter(dftmp.loc[GWAsumstats1.SNPrs,'BP'],-np.log10(GWAsumstats1.P),c=colors1,s=30,rasterized=True,alpha=0.75)
        sns.despine(ax=topAx)# removes top and right lines
        chrs=set(dftmp.CHR) #To make xticks
        medians=[] #postions
        labels=[] #labels
        for chrom in chrs:
            medians.append(np.median(dftmp.loc[dftmp.CHR==chrom,'BP']))
            labels.append(str(chrom))
        #Manhattan 2
        topAx.set_xticks(medians) #add xticklabels
        topAx.set_xticklabels(labels,fontsize=15)
        topAx.set_ylabel('-log10(pvalue)',fontsize=30)
        topAx.set_title(jobName,fontsize=15)
        topAx.axhline(y=-np.log10(5e-8),color='red')
        topAx.axhline(y=-np.log10(1e-5),linestyle='--',color='blue')

        fig.tight_layout()
        fig.savefig('%s_%s_QCed.manhattan.png'%(jobName,phenotype))
        
        
    
print("Finished running %s"%(datetime.datetime.now()))

    
    