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
	
parser = argparse.ArgumentParser(description="Input is the extension of the files to be merged (rowbind). Files to be merged should be in the working directory. Accepts compressions such as gzip and txt etc.")
parser.add_argument('-ext',help="Extension for compiling results, a wild card will be used before the extension",default=None)
parser.add_argument('-wildcardAfter',help="Use also a wildcard after e.g. *.txt*",action = 'store_true')
parser.add_argument('-sep',help="file delimiter ",default='\t')
parser.add_argument('-output',help="merged file output ",default='./')

parser._actions[0].help='Print this help message and exit'
args = parser.parse_args()

if args.ext is None:
	print("A common extension for the files to be merged is required")
	sys.exit()

if not args.wildcardAfter:
	dosageFiles=glob.glob('*%s'%(args.ext))
else:
	dosageFiles=glob.glob('*%s*'%(args.ext))

print("Compiling results this may take a while %s"%(datetime.datetime.now()))
finalDF=pd.concat((pd.read_csv(f,header=0,sep=args.sep,compression='infer') for f in dosageFiles))
print("Finished compiling at %s"%(datetime.datetime.now()))

finalDF.to_csv('%s.merged'%(args.output),sep=args.sep,na_rep='NA')


