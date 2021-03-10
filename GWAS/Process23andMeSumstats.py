#!/usr/bin/env python
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

    
parser = argparse.ArgumentParser(description="Input comprises four files: Sumstats and Annot files as per 23andMe standards")
parser.add_argument('-sumstats',help="Path to the 23andMe sumstats")
parser.add_argument('-gtfile',help="genotyped snps file annot")
parser.add_argument('-allfile',help="genotyped snps file annot")
parser.add_argument('-imfile',help="imputed SNPs annot")
parser.add_argument('-output',help=" file output name",default='GWASsumstats')

parser._actions[0].help='Print this help message and exit'
args = parser.parse_args()

if args.sumstats is None:
    print("Need sumstats file ")
    sys.exit()
    
if args.gtfile is None:
    print("Not all files were provided")
    sys.exit()

if args.imfile is None:
    print("Not all files were provided")
    sys.exit()

if args.allfile is None:
    print("Not all files were provided")
    sys.exit()


print("Reading im, gt and all for obtaining allele data %s"%(datetime.datetime.now()))
gtFile=pd.read_csv(args.gtfile,header=0,sep='\s+',index_col='assay.name',compression='infer',usecols=['assay.name','freq.b'])
imFile=pd.read_csv(args.imfile,header=0,sep='\s+',index_col='assay.name',compression='infer',usecols=['assay.name','freq.b'])
annotFile=pd.read_csv(args.allfile,header=0,sep='\s+',index_col='all.data.id',compression='infer',usecols=['assay.name','all.data.id','scaffold','position','alleles'])
print("Finished reading annot info at %s"%(datetime.datetime.now()))


print("Reading summary statistics %s"%(datetime.datetime.now()))
sumstats=pd.read_csv(args.sumstats,header=0,sep='\s+',compression='infer',usecols=['all.data.id','src','pvalue','effect','stderr','pass','im.num.0','im.num.1','AA.0','AB.0','BB.0','AA.1','AB.1','BB.1'])
print("Finished reading sumstats info at %s"%(datetime.datetime.now()))

sumstats=sumstats.loc[sumstats['pass']=='Y',:]

sumstats['Ncases']=sumstats['im.num.1']
sumstats.Ncases.fillna(sumstats.loc[:,['AA.1','AB.1','BB.1']].sum(axis=1),inplace=True)
sumstats['Nctrls']=sumstats['im.num.0']
sumstats.Nctrls.fillna(sumstats.loc[:,['AA.0','AB.0','BB.0']].sum(axis=1),inplace=True)

sumstats.index=annotFile.loc[sumstats.index,'assay.name']
sumstats.index.name='SNP'
annotFile.index=annotFile['assay.name']
sumstats['CHR']=annotFile.loc[sumstats.index,'scaffold'].str.replace('chr','')
sumstats[['A2','A1']]=annotFile.loc[sumstats.index,'alleles'].str.split('\/',expand=True)
sumstats['BP']=annotFile.loc[sumstats.index,'position']
sumstats['BETA']=sumstats['effect']
sumstats['SE']=sumstats['stderr']
sumstats['P']=sumstats['pvalue']
sumstats['FREQ']=imFile.loc[[i for i in sumstats.index if i in imFile.index],'freq.b']
sumstats.FREQ.fillna(gtFile['freq.b'],inplace=True)





sumstats=sumstats.loc[:,['CHR','BP','A1','A2','FREQ','BETA','SE','P','Ncases','Nctrls']]




sumstats.to_csv('%s_QCed.ctg'%(args.output),sep='\t',na_rep='NA',index=True,)


