#!/usr/bin/env python
#Just calculates lambda based on autosomes
## Code for doing a manhattan plot with python
# Added MarkerName as possible SNPrs synonim

from sys import exit
try:
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['agg.path.chunksize'] = 100000
    import pandas as pd
    from scipy.stats import probplot
    from scipy.stats.morestats import _calc_uniform_order_statistic_medians
    from scipy.stats import uniform
    import re
    import argparse
    import subprocess
    import seaborn
    import os
    import glob
    import time
    import sys
    import datetime
    import socket
    import numpy as np
    import warnings
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy import stats
    from scipy.stats import norm
    from scipy.stats import chi2
except ImportError:
    print('Python module does not seem to be loaded')
    exit()
    
####################Arg parsing #########################
parser = argparse.ArgumentParser(description="input summary statistics, should contain at least the following columns: pval other header names are accepted. SUMSTATS are assumed to be QCed already, specially for imputation. This script includes ambiguous strand and indel removal. Accepted synonims are: pval -> p pvalue pvalues	, snp -> snpid SNPrs rsnumber rs snprs rsid ,	alt -> effect_allele a1  minor_allele effect_allele ea	, reference-> reference_allele ref a2 oa , chr-> CHR chrom chromosome	,	pos-> position bp	,	beta->beta effsize logor log(or) or. This is all capslock insensitive thus BP==bp")
parser.add_argument('input1',help="Input file gwa sumstats 1")
parser.add_argument('-title1','--T1',help="title for plot of gwa sumstats",default="")
parser.add_argument('-o','--output',help="output name without extension",default='QQplot')
parser.add_argument('-namrk','--NAmarker',nargs='+',default=['.'],help='Use this marker(s) as missing value (separate them with a space), default is . use this flag if you missings are not nan or na (e.g. -99 -9 .) It is best to use this flag at the end as it receives an undefinite number of arguments')


parser._actions[0].help="Print this help message and exit"
args = parser.parse_args()

######################### TO TEST INTERACTIVELY #####################
'''
class z():
    def __init__(self):
        self.input1='TASI_noBSB_GWASsumstats_QCed.dat'
        self.NAmarker='.'
        self.output='TASI_qqplot'
args=z()
'''

def qqplot(P,ax=None):
     import seaborn as sns
     import matplotlib.pyplot as plt
     z=norm.ppf(P/2)
     lambGC = round(np.median(z**2)/chi2.ppf(0.5,df=1),3)
     p = 2*norm.cdf(-abs(z))
     p.sort()
     expected = list(range(len(p)))
     lobs = np.array(-np.log10(p))
     lexp = np.array(-np.log10(np.array(expected) / (len(expected)+1)))
     # plots all points with p < 0.05
     if ax is None:
         fig,ax = plt.subplots(1)
     p_sig = [i for i in p if i<0.05]
     ax.plot(lexp, lobs, 'd', color='steelblue',alpha=0.75,label='Observed')
     xmin,xmax=ax.get_xlim()
     ax.plot([0,xmax+0.2], [0,xmax+0.2],linestyle='dashed', color='black',label='expected')
     ax.set_xlim(xmin,xmax)
     ax.set_xlabel('-log10(expected P)')
     ax.set_ylabel('-log10(observed P)')
     ax.text(0.2,6,'$\lambda = %s$'%(str(lambGC)),fontsize=16)
     sns.despine()
     return(ax)


### Print and check for errors.
print("By Adrian Campos @ QIMR Berghofer\n please let me know of any problems errors or modifications\nadrian.campos@qimrberghofer.edu.au")
print("This script is provided as is with no warranty and under under a CC license ...")
if socket.gethostname()=='hpcpbs01.adqimr.ad.lan':
	print("This script opens sumstats to memory ~50Gb, run it interactively or with qsub!")
	raise RuntimeError('You should run this script interactively or via qsub, not in login nodes')

###################File preprocessing ####################
#Open summary statistics

inputGWAS1=args.input1
print("Will do a QQplot, using the following GWA summary statistics:\n%s"%(inputGWAS1))
print("Interpreting column headers for %s"%(inputGWAS1))
colConverter={}

####### NOTE that in this code effect allele equals the allele that has an effect on the phenotype, whereas the reference allele is considered the basis (dosage of 0)
HeaderSynonims={'pval':'pval','p':'pval','p.value':'pval','pvalue':'pval','pvalues':'pval','pvalues':'pval','p_bolt_lmm':'pval','mtag_pval':'pval','p-value':'pval',
				'chr':'CHR','chrom':'CHR','chromosome':'CHR'}
				
dtypeDict={'pval':object,'CHR':object}
numericCols=['pval','CHR']				
tmpDF=pd.read_table(inputGWAS1,header=0,nrows=500,sep=r"\s+")
colnames=[i for i in tmpDF.columns]
colConverter={i:HeaderSynonims.get(i.lower(),None) for i in colnames}
colConverterRev={x:y for y,x in colConverter.items()}
ColsToUse=[i for i in colConverter if not colConverter[i]==None]
dtypes= {i:dtypeDict.get(colConverter[i],None) for i in colConverter}
toDel=[i for i in dtypes if dtypes[i]==None]
for key in toDel:
	del(dtypes[key])

    
print("Interpreted the columns in the following way:")
if len(ColsToUse)<1:
	for i in colConverter:
		print(i,colConverter[i])
	raise RuntimeError('The GWA summary statisitcs header could not be interpreted, please refer to the documentation')
elif len(ColsToUse)>1:
	for i in colConverter:
		print(i,colConverter[i])
	warnings.warn('Two or more columns mapped to the same interpretation, one will be arbitrarily chosen.')
	colConverter={x:y for y,x in colConverterRev.items() if not y==None}
	ColsToUse=list(colConverter.keys())
for i in colConverter:
	if not colConverter[i]==None:
		print(i,colConverter[i])


print("Opening summary statistics")
naMarkers=['NA','NAN','na','nan']
naMarkers.extend(args.NAmarker)
GWAsumstats1 = pd.read_table(inputGWAS1,header=0,index_col=None,usecols=ColsToUse,sep="\s+",dtype=dtypes,na_values=naMarkers)    
GWAsumstats1.columns=[colConverter[i] for i in GWAsumstats1.columns]
for col in GWAsumstats1.columns:
    if col in numericCols:
        GWAsumstats1.loc[:,col]=pd.to_numeric(GWAsumstats1.loc[:,col],errors='coerce')

print("Read %s snps from the summary statistics"%(len(GWAsumstats1.index)))		
GWAsumstats1.sort_values(by='pval',inplace=True)
GWAsumstats1=GWAsumstats1.loc[GWAsumstats1.CHR<23,:]

##Plot

fig = plt.figure() 
ax=fig.add_subplot(111)
qqplot(GWAsumstats1.pval.dropna(),ax)
fig.tight_layout()
fig.savefig(args.output+'.png')


