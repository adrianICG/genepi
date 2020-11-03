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
	
parser = argparse.ArgumentParser(description="Input is the extension of the files to be merged (rowbind). Files to be merged should be in the working directory. Accepts compressions such as gzip and txt etc. NOTE this program will require enough memory to load all of your files")
parser.add_argument('-ext',help="Extension for compiling results, a wild card will be used before the extension",default='sscore')
parser.add_argument('-output',help="merged file output ",default='./prs')
parser.add_argument('-PRSice',help="PRSice stylke output",action = 'store_true')


parser._actions[0].help='Print this help message and exit'
args = parser.parse_args()

if args.ext is None:
	print("A common extension for the files to be merged is required")
	sys.exit()

scoreFiles=glob.glob('*%s'%(args.ext))




print("Compiling results this may take a while %s"%(datetime.datetime.now()))
firstFile=scoreFiles.pop()
if args.PRSice:
	cols2use=None
else:
	cols2use=['IID','ALLELE_CT','SCORE1_SUM']

PRSDF=pd.read_csv(firstFile,header=0,sep='\s+',compression='infer',usecols=cols2use)
PRSDF.index=PRSDF.IID
PRSDF.drop('IID',axis=1,inplace=True)
PRSDF.drop('FID',axis=1,inplace=True)

for f in scoreFiles:
	tmpDF=pd.read_csv(f,header=0,sep='\s+',compression='infer',usecols=cols2use)
	tmpDF.index=tmpDF.IID
	tmpDF.drop('IID',axis=1,inplace=True)
	tmpDF.drop('FID',axis=1,inplace=True)
	PRSDF=PRSDF+tmpDF
print("Finished compiling at %s"%(datetime.datetime.now()))

PRSDF.to_csv('%s.sumscores'%(args.output),sep='\t',na_rep='NA',index=True)

