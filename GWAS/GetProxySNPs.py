#!/usr/bin/env python

import argparse
import re

parser = argparse.ArgumentParser(description="input a list of SNPs")
parser.add_argument('SNPlist',help="Input meta analysis file default columns")
parser.add_argument('-Nproxies',help="number of proxies per SNP",default=3,type=int)
parser.add_argument('-o','--output',help="output name ",default='proxyList')

parser._actions[0].help="Print this help message and exit"
args = parser.parse_args()

print("Opening SNPlist file: %s"%(args.SNPlist))
SNPlist=[]
with open(args.SNPlist,'r') as infile:
	for line in infile:
		SNPlist.append(line.strip())

print("Read %s SNPs from meta analysis file"%(len(SNPlist)))

with open(args.output,'w') as outfile:
	for SNP in SNPlist:
		r= requests.get('https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=%s&pop=CEU&r2_d=r2&token=22141f405e7d'%(SNP))
		lines=r.text.split('\n')
		if re.search("error",lines[0]):
			print("Could not find SNP %s"%(SNP))
			continue
		for i in range(args.Nproxies):
			try:
				outfile.write(lines[i+1])
			except IndexError:
				print("SNP %s had less than %s proxy SNPs"%(SNP,args.Nproxies))
				
print("Finished requesting proxies for %s SNPs"%(len(SNPlist)))
	
	