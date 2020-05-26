#!/usr/bin/env python
#Update 25/03/20 
#Allows release 10 adding of imputation covariates
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
except ImportError:
    print('Python module does not seem to be loaded')
    exit()

####################Arg parsing #########################
parser = argparse.ArgumentParser(description="input positional parameters: PRS file, phenotypes file (includes covars and qcovars) and variable description file. Check --help for required IDs")
parser.add_argument('prsfile',help="PRS file: FID IID S1 S2 ....")
parser.add_argument('phenofile',help="Phenotypes file: FID IID Age Sex Covar1 ... Covarx ... Pheno1 Pheno2 Pheno3 ... PhenoK\n Missing values shopuld be NA or NAN not -99 nor -9 nor .")
parser.add_argument('covdescript',help="Describe the variables on phenofile: ColumnName [covar|qcovar]")
parser.add_argument('jobname',help="Job Name")
parser.add_argument('-r9','--UseRelease9',help="Use release9 grm (default is release10)",action = 'store_true')
parser.add_argument('-std','--Standardize',help="Scale the phenotypes to mean=0 std =1",action = 'store_true')
parser.add_argument('-grm','--UseThisGRM',help="Use a specific GRM (useful for diagplus in R9) must be full path!")
parser.add_argument('-grmext','--GRMextension',help="Sepcify the GRM extension, by default gz",default='gz')
parser.add_argument('-noimput','--DoNotAddImput',help="Do not include imputation covariates automatically (use if they are already in your phenofile)",action = 'store_true')
parser.add_argument('-noPC','--DoNotAddPC',help="Do not include PC covariates automatically (use if they are already in your phenofile)",action = 'store_true')
parser.add_argument('-noStdID','--NotStdIDs',help="Do not standardize the twin ids IIDs to length 7",action = 'store_true')
parser.add_argument('-sexage','--SexAgeInt',help="automatically add sex * age interaction as qcovariate",action = 'store_true')
parser.add_argument('-agesqrd','--AgeSquared',help="automatically add age squared as qcovariate",action = 'store_true')
parser.add_argument('-agesqrdsex','--AgeSquaredSex',help="automatically add age squared time sex as qcovariate",action = 'store_true')
parser.add_argument('-duplicates','--DealDuplicates',help="automatically deal with dups valid is nothing (by default), remove or average",default = None)
parser.add_argument('-namrk','--NAmarker',nargs='+',default=['.'],help='Use this marker(s) as missing value (separate them with a space), default is . use this flag if you missings are not nan or na (e.g. -99 -9 .) It is best to use this flag at the end as it receives an undefinite number of arguments')


parser._actions[0].help='Print this help message and exit'
args = parser.parse_args()
######################### TO TEST INTERACTIVELY #####################
#  1.  Fill in the fields below with your input data
#  2. Start an interactive hpc session (load python module and start ipython) and run the first 20 lines of this file (skip from 21 to 36)
#  3. Run the par below without including the ''' lines
#  4. Keep running everything below  and observing the output
'''
class z():
    def __init__(self):
        self.prsfile='StdIndividual_PRS_MDDr10.txt'
        self.phenofile='vars.txt'
        self.covdescript='COVS'
        self.jobname='DepPRS_nickvars'
        self.UseRelease9=False
        self.Standardize=False
        self.UseThisGRM=None
        self.GRMextension='gz'
        self.DoNotSubmit=False
        self.DoNotAddImput=True
        self.DoNotAddPC=False
        self.NotStdIDs=True
        self.SexAgeInt=True
        self.AgeSquared=False
        self.AgeSquaredSex=False
        self.DealDuplicates='remove'
        self.NAmarker=['.']
args=z()
'''
#########################Pre needed functions #######################

normalwd=os.getcwd()

def CheckCompletion(jobDic,timew=300):
    '''Check completion function for submitted jobs to the hpc cluster
    input a job dictionary of job IDs and zeros (means they are running), and a predetermined wait time between cluster query'''
    ids=[i for i in jobDic.keys()] #get all job ids
    submstr=str("for job in %s; do qstat ${job}; done"%(' '.join(ids))) #submit query for all jobs
    while sum(jobDic.values())<len(jobDic): #now understand why whiles exist
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
    return(1) #return something just in case you need to know
    
    
def SubmitAllJobs(ext='.PBS'):
    ''' Submits all files with an extension in the current working directory'''
    import subprocess
    import re
    # Submit the instruction to run all .PBS scripts on cluster
    Submits=subprocess.Popen('for script in *%s; do qsub ${script}; done'%(ext),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    #Obtain the output from the submissions (to parse job ids)
    outlines=Submits.communicate()[0].decode("utf-8")
    # find the job ids and store them in a has table 
    jobs={i.split('.')[0]:0 for i in re.split(pattern='\n',string=outlines) if not i==''}
    return(jobs) #return the completion hash table (to be used to check completion if needed)
    
    
def CheckErrors(ext='*.e*'):
    '''
        Error checker, basically greps Error in all of the files *.e* of the currwd
    '''
    import subprocess
    import re
    # Submit the grep command looking for errors reported by HPC
    Submits=subprocess.Popen("grep 'Error' %s"%(ext),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    # Process the output
    outlines=Submits.communicate()[0].decode("utf-8")
    #Store all errors and return them
    Errors=[i for i in re.split(pattern='\n',string=outlines) if not i=='']
    return(Errors)

    
###################File preprocessing ####################
PRSfile=args.prsfile  #'EA3PRSScores02.tab' 
Phenofile=args.phenofile  #'GCTA_data.tab'
jobName=args.jobname    #'SelfHarmBeh'
covdescriptfile=args.covdescript
if args.UseRelease9:
    grmfile='/reference/genepi/GWAS_release/Release9/Release9_Observed/GRM/GRM_GCTA_autosomes'
else:
    grmfile='/reference/genepi/GWAS_release/Release10/Release10_Observed/GRM/GRM_GCTA_autosomes_GSA_European'
    
if args.UseThisGRM:
    grmfile=args.UseThisGRM
        
grmext=args.GRMextension
#open both files
NaMarkers=[]
for element in args.NAmarker:
    try:
        NaMarkers.append(float(element))
    except ValueError:
        NaMarkers.append(element)
print("The following values will be taken as NA")
for i in NaMarkers:
    print(i)

prsdf=pd.read_table(PRSfile,header=0,sep=r"\s+",dtype={'IID':str,'FID':str},na_values=NaMarkers) #check if IID or iid before, or use try, except
prsdf.set_index('IID', inplace=True)
phenodf=pd.read_table(Phenofile,header=0,sep=r"\s+",dtype={'IID':str,'FID':str},na_values=args.NAmarker)
phenodf.set_index('IID', inplace=True)

#ensure 7 digit IDs
if not args.NotStdIDs:
    phenodf.index=['0000000'[0:7-len(str(i))]+str(i) for i in phenodf.index] #adding 0's to the left of IDs
    prsdf.index=['0000000'[0:7-len(str(i))]+str(i) for i in prsdf.index]

#Open covar names
covarnames=pd.read_table(covdescriptfile,header=None,sep=r"\s+",index_col=0)
covarDict= covarnames.to_dict()[1]
print ('Opened all of the input files, attempting to match them:')

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

#match them
availablePhenoID=[i for i in phenodf.index if i in prsdf.index]
print("Of the %s IDs in phenotype %s matched to prs table"%(len(phenodf.index),len(availablePhenoID)))
prsdf=prsdf.loc[availablePhenoID]
phenodf=phenodf.loc[availablePhenoID]



print("Creating and changing to lmm working directory (to submit from there GCTA scripts)")

#create lmm working subdirectory
subprocess.call(['mkdir','lmm_calc'])
os.chdir(normalwd+'/lmm_calc') # move there


## Add imputation covariates make (optional)
if not args.DoNotAddImput:
    print("Adding imputation covariates")
    if args.UseRelease9:
        imputcovfile='/reference/genepi/GWAS_release/Release9/Release9_HRCr1.1/info/GWAS_ImputationRunCovariates.ped'
    else:
        imputcovfile='/reference/genepi/GWAS_release/Release10/Release10_HRCr1.1/info/GWAS_ImputationRunCovariates.ped'
    imputcovdf=pd.read_table(imputcovfile,header=None,index_col=1,sep=r"\s+")

## Add PC qcovars (optional)
if not args.DoNotAddPC:
    PCcols=['PC'+str(i) for i in [1,2,3,4,5,6,7,8,9,10]]
    PCcols.append('INDID')
    print("Adding PC qcovariates")
    if args.UseRelease9:
        pcdf=pd.read_table('/reference/genepi/GWAS_release/Release9/Release9_Observed/AncestryChecks/GWAS_filtered_PrincipalComponentsScores.txt',header=0,usecols=PCcols,sep=r"\s+",index_col='INDID')
    else:
        pcdf=pd.read_table('/reference/genepi/GWAS_release/Release10/Release10_Observed/AncestryChecks/GWAS_PrincipalComponentsScores.txt',header=0,usecols=PCcols,sep=r"\s+",index_col='INDID')
    for i in pcdf.columns:
        phenodf[i]=pcdf.loc[phenodf.index,i]
        print("%s Subjects will be removed because they had no PC information"%(str(len([i for i in phenodf.index if not i in pcdf.index]))))
        covarDict[i]='qcovar'
        
## VARIABLES TO CREATE INPUT FILES:
phenodf['IID']=phenodf.index
covariates=[i for i in covarDict if covarDict[i]=='covar']
qcovariates=[i for i in covarDict if covarDict[i]=='qcovar']
orderedCols=['FID','IID']

### Obtain phenotypes (not covars nor qcovars)
orderedCols.extend([i for i in phenodf.columns if i not in covariates and i not in orderedCols and i not in qcovariates])

print("Creating GCTA input files: pheno")
######### Create GCTA pheno inputs ##############################
phenos=[]
for i in range(2,len(orderedCols)):
    tmpdf=phenodf[['FID','IID',orderedCols[i]]]
    tmpdf=tmpdf.dropna()
    if args.Standardize==True:
        phenodf.loc[:,orderedCols[i]]=(phenodf.loc[:,orderedCols[i]]-np.mean(phenodf.loc[:,orderedCols[i]]))/np.std(phenodf.loc[:,orderedCols[i]])
        tmpdf.loc[:,orderedCols[i]]=(tmpdf.loc[:,orderedCols[i]]-np.mean(tmpdf.loc[:,orderedCols[i]]))/np.std(tmpdf.loc[:,orderedCols[i]])
    tmpdf.dropna(inplace=True)
    tmpdf.to_csv(jobName+orderedCols[i]+'.pheno',index=False,sep=" ",na_rep='NA')
    phenos.append(orderedCols[i])

print("Creating GCTA input files: covar")
########## Create GCTA covar file (only one needed) ##################
cols=['FID','IID']
cols.extend(covariates)
tmpdf=phenodf[cols]
if not args.DoNotAddImput:
    for i in imputcovdf.loc[:,5:].columns:
        tmp=imputcovdf.loc[availablePhenoID,i]
        tmp.dropna(inplace=True)
        if len(set(tmp.values))<2:
            print("One of the imputation covariates was a constant, it wont be added as GCTA would fail\nThis probably means your sample only contains individuals genotyped in two chips.")
            continue
        tmpdf['ImputRun'+str(i)]=imputcovdf.loc[availablePhenoID,i] #matching by IID
        noImputRunInfo=sum(tmpdf['ImputRun'+str(i)].isna())
        if noImputRunInfo:
            print("No imputation info detected for %s individuals check the GWAS_ImputationRunCovariate.ped"%(noImputRunInfo))
        
tmpdf=tmpdf.dropna()
tmpdf.to_csv(jobName+'.covariates',index=False,sep=" ",na_rep='NA')


print("Creating GCTA input files: qcovars")
#########Crate GCTA qcovars (adding the PRS to each) ##################
nonPrs=['FID','IID','PHENO','CNT']
PRS=[i for i in prsdf.columns if i not in nonPrs]
cols=['FID','IID']
cols.extend(qcovariates)
tmpdf=phenodf[cols]
cutoffs=[]
for score in PRS:
    #reorder the columns because we want PRS to be first qcovar to find it systematically later
    tmpcols=['FID','IID']
    tmpdf2=tmpdf.copy()
    tmpdf2[score]=prsdf[[score]].loc[tmpdf2.index]
    tmpcols.append(score)
    tmpcols.extend([i for i in tmpdf2.columns if i not in tmpcols])
    tmpdf2=tmpdf2[tmpcols]
    # finish the reordering
    if args.SexAgeInt:
        print('Adding Sex * Age interactions as covariate for %s'%(score))
        if 'sex' not in [i.lower() for i in phenodf.columns]:
            print('No sex column detected, thus not adding the interaction')
        elif 'age' not in [i.lower() for i in phenodf.columns]:
            print('No age column detected, thus not adding the interaction')
        else:
            colsex=[i for i in phenodf.columns if i.lower()=='sex'][0]
            colage=[i for i in phenodf.columns if i.lower()=='age'][0]
            tmpdf2['SexAge']=phenodf.loc[:,colsex]*phenodf.loc[:,colage]
    if args.AgeSquared:
        print('Adding age squared as a covariate for %s'%(score))
        if 'age' not in [i.lower() for i in phenodf.columns]:
            print('No age column detected, thus not adding the interaction')
        else:
            colage=[i for i in phenodf.columns if i.lower()=='age'][0]
            tmpdf2['AgeSquared']=phenodf.loc[:,colage]*phenodf.loc[:,colage]
    if args.AgeSquaredSex:
        print('Adding Sex * Age^2 interactions as covariate for %s'%(score))
        if 'sex' not in [i.lower() for i in phenodf.columns]:
            print('No sex column detected, thus not adding the interaction')
        elif 'age' not in [i.lower() for i in phenodf.columns]:
            print('No age column detected, thus not adding the interaction')
        else:
            colage=[i for i in phenodf.columns if i.lower()=='age'][0]
            phenodf.loc[:,colage]*phenodf.loc[:,colage]
            tmpdf2['AgeSquaredSex']=tmpdf2['AgeSquared']*phenodf.loc[:,colsex]
    tmpdf2=tmpdf2.dropna()
    tmpdf2.to_csv(jobName+score+'.qcovariates',index=False,sep=" ",na_rep='NA')
    cutoffs.append(score)

############### Create GCTA scripts:###########################
print("Creating GCTA output directory lmm_calc_out")
subprocess.call(['mkdir','lmm_calc_out'])

print("Creating GCTA PBS scripts")
shebang=r'#!/bin/bash' #header to add to each file consider env bash
covarfile=jobName+'.covariates'


mem=64
for phenotype in phenos:
    for score in cutoffs:
        currname='_'.join([phenotype,score])
        currphenofile=jobName+phenotype+'.pheno'
        currqcovar=jobName+score+'.qcovariates'
        with open(currname+'.PBS','w') as currscript:
            currscript.write(shebang+'\n') #write header
            currscript.write("#PBS -N lmm%s\n#PBS -r n\n#PBS -l ncpus=4,mem=%sGB,walltime=10:00:00\nmodule load GCTA/1.91.7beta\ncd $PBS_O_WORKDIR\n"%(currname,mem)) #write pbs options
            currscript.write("gcta64 --reml --grm-%s %s --pheno %s --qcovar %s --covar %s --threads 4 --out lmm_calc_out/%s --reml-est-fix"%(grmext,grmfile,currphenofile,currqcovar,covarfile,currname))

###################### submit all jobs #######################
print("Submitting jobs to the cluster GCTA PBS scripts")
finishedJobs={}
print("Started submitting %s"%(datetime.datetime.now()))
finishedJobs=SubmitAllJobs(ext='.PBS') # hash table with completion info
print("Jobs submitted %s waiting for job completion prior to continue"%(datetime.datetime.now()))
#####################Check completion and errors ##############
CheckCompletion(finishedJobs)
print("jobs finished running approx %s"%(datetime.datetime.now()))
# Checking for errors using simple grep
anyErrors=CheckErrors()
if len(anyErrors):
    print("The following errors were reported by HPC: (will try to continue)\n")
    for currError in anyErrors:
        print(currError)
else:
    print("No errors reported by hpc GCTA MAY have reported errors in its logs")
            
################# Compiling results #############################
            
print("Compiling GCTA results")            
os.chdir('lmm_calc_out')
hsqFiles=glob.glob('*.hsq')

########### Creating hash tables to store values #######
FinaldictFixedEff={}
FinaldictSE={}
FinaldictR2={}
FinaldictR2Low={}
FinaldictR2Upp={}
FinaldictPval={}

#Parse each output file to extract the PRS Feff, SE, N and calculate R2 and Pval
for currhsqname in hsqFiles:
    with open(currhsqname,'r') as currhsq: #open file
        name=re.split('_',currhsqname) #obtain name
        currpheno='_'.join(name[0:len(name)-1])
        currpheno=re.sub("_Std", "", currpheno)
        currcutoff=name[len(name)-1].replace('.hsq','')
        if re.search("_Std",currhsqname):
            currcutoff="Std_"currcutoff
        if currpheno not in FinaldictFixedEff:
        ## Initialize values on the hash tables
            FinaldictFixedEff[currpheno]={currcutoff:None}
            FinaldictSE[currpheno]={currcutoff:None}
            FinaldictR2[currpheno]={currcutoff:None}
            FinaldictR2Low[currpheno]={currcutoff:None}
            FinaldictR2Upp[currpheno]={currcutoff:None}
            FinaldictPval[currpheno]={currcutoff:None}
        else:
        ## Or store them
            FinaldictFixedEff[currpheno][currcutoff]=None
            FinaldictSE[currpheno][currcutoff]=None
            FinaldictR2[currpheno][currcutoff]=None
            FinaldictR2Low[currpheno][currcutoff]=None
            FinaldictR2Upp[currpheno][currcutoff]=None
            FinaldictPval[currpheno][currcutoff]=None
        currFixedEfflines=[]
        foundHeader=0
        N=None
        tstat=None
        for line in currhsq: # for line in current GCTA output fin N and the Fix_eff lines and store them:
            if re.match('n',line):
                N=np.float(re.split('\t',line)[1])
            if re.match('Fix_eff',line):
                foundHeader=1
            if not foundHeader:
                continue
            currFixedEfflines.append(line) # saving all lines (we just care for position 2 thats were our PRS is as the first is intercept)
        FinaldictFixedEff[currpheno][currcutoff]=float(re.split('\t',currFixedEfflines[2])[0]) # Get the Feff we care for
        FinaldictSE[currpheno][currcutoff]=float(re.split('\t',currFixedEfflines[2])[1]) # Get the SE we care for
        tstat=FinaldictFixedEff[currpheno][currcutoff]/FinaldictSE[currpheno][currcutoff] # calculate tstat
        FinaldictPval[currpheno][currcutoff]=stats.t.sf(np.abs(tstat), N-1)*2 #pvalue two sided
        FinaldictR2[currpheno][currcutoff]=(FinaldictFixedEff[currpheno][currcutoff]/np.std(phenodf.loc[:,currpheno].astype(float))*np.std(prsdf.loc[:,currcutoff]))**2 #R2
        lowFE=FinaldictFixedEff[currpheno][currcutoff]-1.96*FinaldictSE[currpheno][currcutoff]
        FinaldictR2Low[currpheno][currcutoff]=(lowFE/np.std(phenodf.loc[:,currpheno].astype(float))*np.std(prsdf.loc[:,currcutoff]))**2 #R2Low
        uppFE=FinaldictFixedEff[currpheno][currcutoff]+1.96*FinaldictSE[currpheno][currcutoff]
        FinaldictR2Upp[currpheno][currcutoff]=(uppFE/np.std(phenodf.loc[:,currpheno].astype(float))*np.std(prsdf.loc[:,currcutoff]))**2 #R2
#### Convert to tables
        
PRS_Fixed_eff=pd.DataFrame.from_dict(FinaldictFixedEff)
PRS_SE=pd.DataFrame.from_dict(FinaldictSE)
PRS_R2=pd.DataFrame.from_dict(FinaldictR2)
PRS_R2Low=pd.DataFrame.from_dict(FinaldictR2Low)
PRS_R2Upp=pd.DataFrame.from_dict(FinaldictR2Upp)
PRS_Pval=pd.DataFrame.from_dict(FinaldictPval)

#############Save output ##################
os.chdir(normalwd)
print("Saving output on the working directory:")
print(os.getcwd())
PRS_Fixed_eff.to_csv(jobName+'Feff.tab',sep='\t',na_rep='NA')
PRS_SE.to_csv(jobName+'FeffSE.tab',sep='\t',na_rep='NA')
PRS_R2.to_csv(jobName+'R2.tab',sep='\t',na_rep='NA')
PRS_R2Low.to_csv(jobName+'R2Low.tab',sep='\t',na_rep='NA')
PRS_R2Upp.to_csv(jobName+'R2Upp.tab',sep='\t',na_rep='NA')
PRS_Pval.to_csv(jobName+'Pval.tab',sep='\t',na_rep='NA')

print("finished running at %s"%(datetime.datetime.now()))

#END of script
    
    
    
    
    
    
    
    
    
    
    