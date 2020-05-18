#!/usr/bin/env python
## Code for doing a manhattan plot with python
# Added MarkerName as possible SNPrs synonim

from sys import exit
try:
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['agg.path.chunksize'] = 100000
    import pandas as pd
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
except ImportError:
    print('Python module does not seem to be loaded')
    exit()
    
####################Arg parsing #########################
parser = argparse.ArgumentParser(description="input summary statistics, should contain at least the following columns: SNPrs,CHR,BP,REF,ALT,beta,pval other header names are accepted. SUMSTATS are assumed to be QCed already, specially for imputation. This script includes ambiguous strand and indel removal. Accepted synonims are: pval -> p pvalue pvalues	, snp -> snpid SNPrs rsnumber rs snprs rsid ,	alt -> effect_allele a1  minor_allele effect_allele ea	, reference-> reference_allele ref a2 oa , chr-> CHR chrom chromosome	,	pos-> position bp	,	beta->beta effsize logor log(or) or. This is all capslock insensitive thus BP==bp")
parser.add_argument('input1',help="Input file gwa sumstats 1")
parser.add_argument('-title1','--T1',help="title for plot of gwa sumstats",default="")
parser.add_argument('-color1','--C1',help="Color1 ",default='#999795')
parser.add_argument('-color2','--C2',help="Color2 ",default='#a92e4a')
parser.add_argument('-sep',help="file delimiter ",default='\s+')
parser.add_argument('-circular',help="Experimental plot ciruclar manhattan plot",action='store_true')
parser.add_argument('-o','--output',help="output name without extension",default='Manhattan')
parser.add_argument('-namrk','--NAmarker',nargs='+',default=['.'],help='Use this marker(s) as missing value (separate them with a space), default is . use this flag if you missings are not nan or na (e.g. -99 -9 .) It is best to use this flag at the end as it receives an undefinite number of arguments')


parser._actions[0].help="Print this help message and exit"
args = parser.parse_args()

######################### TO TEST INTERACTIVELY #####################
'''
class z():
    def __init__(self):
        self.input1='daner_PGC_BIP32b_mds7a.loo.no.bmau'
        self.T1=''
        self.NAmarker='.'
        self.output='BIP32b'
args=z()
'''


### Print and check for errors.
print("By Adrian Campos @ QIMR Berghofer\n please let me know of any problems errors or modifications\nadrian.campos@qimrberghofer.edu.au")
print("This script is provided as is with no warranty and under under a CC license ...")
if socket.gethostname()=='hpcpbs01.adqimr.ad.lan':
	print("This script opens sumstats to memory ~50Gb, run it interactively or with qsub!")
	raise RuntimeError('You should run this script interactively or via qsub, not in login nodes')

###################File preprocessing ####################
#Open summary statistics

inputGWAS1=args.input1
print("Will do a Manhattan plot, using the following GWA summary statistics:\n%s"%(inputGWAS1))
print("Interpreting column headers for %s"%(inputGWAS1))
colConverter={}

####### NOTE that in this code effect allele equals the allele that has an effect on the phenotype, whereas the reference allele is considered the basis (dosage of 0)
HeaderSynonims={'pval':'pval','p':'pval','p.value':'pval','pvalue':'pval','pvalues':'pval','pvalues':'pval','p_bolt_lmm':'pval','mtag_pval':'pval','p-value':'pval','pval_estimate':'pval',
                'markername':'SNPrs','snp':'SNPrs','snpid':'SNPrs','rsnumber':'SNPrs','rs':'SNPrs', 'snprs':'SNPrs','rsid':'SNPrs',
				'chr':'CHR','chrom':'CHR','chromosome':'CHR',
				'pos':'BP','position':'BP','bp':'BP','bpos':'BP'}
				
dtypeDict={'pval':object,'SNPrs':object,'effect_allele':object,'reference_allele':object,'CHR':object,'BP':object,'beta':object,'or':object}
numericCols=['pval','CHR','BP']				
tmpDF=pd.read_table(inputGWAS1,header=0,nrows=500,sep=args.sep)
colnames=[i for i in tmpDF.columns]
colConverter={i:HeaderSynonims.get(i.lower(),None) for i in colnames}
colConverterRev={x:y for y,x in colConverter.items()}
ColsToUse=[i for i in colConverter if not colConverter[i]==None]
dtypes= {i:dtypeDict.get(colConverter[i],None) for i in colConverter}
toDel=[i for i in dtypes if dtypes[i]==None]
for key in toDel:
	del(dtypes[key])
	
transformEff=0
if 'or' in colConverter.values():
	print("Detected an Odds ratio column, checking first 500 elements")
	if sum(tmpDF[colConverterRev['or']]<0):
		print("Detected at least one negative value, assuming OR to be log(OR)")
		transformEff=0
	else:
		print("Detected only positive values, will transform OR to log(OR)")
		transformEff=1
	key=[i for i in colConverter if colConverter[i]=='or'][0]
	colConverter[key]='beta'

print("Interpreted the columns in the following way:")
if len(ColsToUse)<3:
	for i in colConverter:
		print(i,colConverter[i])
	raise RuntimeError('The GWA summary statisitcs header could not be interpreted, please refer to the documentation')
elif len(ColsToUse)>3:
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
GWAsumstats1 = pd.read_table(inputGWAS1,header=0,index_col=None,usecols=ColsToUse,sep=args.sep,dtype=dtypes,na_values=naMarkers)    
GWAsumstats1.columns=[colConverter[i] for i in GWAsumstats1.columns]
GWAsumstats1['SNPvariant']=GWAsumstats1['CHR'].apply(str)+':'+GWAsumstats1['BP'].apply(str)
for col in GWAsumstats1.columns:
    if col in numericCols:
        GWAsumstats1.loc[:,col]=pd.to_numeric(GWAsumstats1.loc[:,col],errors='coerce')

print("Read %s snps from the summary statistics"%(len(GWAsumstats1.index)))		


#Sort

GWAsumstats1.sort_values(by=['CHR','BP'],inplace=True)


# calculate new positions for plot:
dftmp=GWAsumstats1
dftmp=dftmp.loc[:,['SNPrs','CHR','BP','SNPvariant']]
dftmp.drop_duplicates(inplace=True)

if sum(dftmp.SNPrs.duplicated()):
    print("There are duplcated variants in sumstats with different positions. One will randomly be removed")
    dftmp=dftmp.loc[~dftmp.SNPrs.duplicated(),:]

pastlen=len(dftmp)
dftmp.dropna(inplace=True)
GWAsumstats1.dropna(inplace=True)
if len(dftmp)!=pastlen:
    print("There were NAs in the sumstats, they were removed automatically but you might want to double check")

chrs=set(dftmp.CHR) #make sure its sorted
pastChr=chrs.pop()
maxBP=max(dftmp.loc[dftmp.CHR==pastChr,'BP'])
for chr in chrs:
    tmpBP=dftmp.loc[dftmp.CHR==chr,'BP']
    tmpBP=tmpBP+maxBP
    dftmp.loc[tmpBP.index,'BP']=tmpBP
    maxBP=max(dftmp.loc[dftmp.CHR==chr,'BP'])

dftmp.set_index(dftmp.SNPrs,inplace=True)


##Plot
if not args.circular:
	fig = plt.figure() 
	topAx=fig.add_subplot(111)
	fig.set_size_inches(25,10)
	#Manhattan 1
	#colors1=['#424242' if i%2 else '#a8a7a5' for i in GWAsumstats1.CHR.values] #SCZ colors
	colors1=[args.C1 if i%2 else args.C2 for i in GWAsumstats1.CHR.values] #BIP colors
	topAx.scatter(dftmp.loc[GWAsumstats1.SNPrs,'BP'],-np.log10(GWAsumstats1.pval),c=colors1,s=30,rasterized=True)
	sns.despine(ax=topAx)# removes top and right lines
	chrs=set(dftmp.CHR) #To make xticks
	medians=[] #postions
	labels=[] #labels
	for chr in chrs:
		medians.append(np.median(dftmp.loc[dftmp.CHR==chr,'BP']))
		labels.append(str(chr))
	#Manhattan 2
	topAx.set_xticks(medians) #add xticklabels
	topAx.set_xticklabels(labels,fontsize=15)
	topAx.set_ylabel('-log10(pvalue)',fontsize=30)
	topAx.set_title(args.T1,fontsize=15)
	topAx.axhline(y=-np.log10(5e-8),color='red')
	topAx.axhline(y=-np.log10(1e-5),linestyle='--',color='blue')

	fig.tight_layout()
	fig.savefig(args.output+'.svg')
else:
	fig = plt.figure()
	max_theta = 2.0 * np.pi
	number_points = len(GWAsumstats1.pval)
	topAx=fig.add_subplot(111,projection='polar')
	fig.set_size_inches(15,15)
	#Manhattan 1
	#colors1=['#424242' if i%2 else '#a8a7a5' for i in GWAsumstats1.CHR.values] #SCZ colors
	colors1=[args.C1 if i%2 else args.C2 for i in GWAsumstats1.CHR.values] #BIP colors
	theta = np.linspace(0.0, max_theta, number_points)
	r =-np.log10(GWAsumstats1.pval)
	topAx.scatter(theta, r,c=colors1,s=30,rasterized=True)
	topAx.plot(theta,[-np.log10(5e-8) for i in theta],color='red',linestyle='solid')
	topAx.set_rorigin(-10)
	for key in topAx.spines:
		topAx.spines[key].set_visible("False")
		topAx.set_xticks([])
		topAx.set_yticks([])
	#topAx.set_ylim(0,10)
	fig.tight_layout()
	fig.savefig(args.output+'.svg')

print("Finished running at %s"%(datetime.datetime.now()))




