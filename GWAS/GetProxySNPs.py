#!/usr/bin/env python

import argparse
import re
import requests
import sys

'''
class z():
	def __init__(self):
		self.SNPlist='IndSigHits'
		self.token=None
		self.Nproxies=3
		self.output='listFor23andMe.txt'
args=z()
'''

parser = argparse.ArgumentParser(description="input a list of SNPs and the LDlink token to access the API. See help for more information")
parser.add_argument('SNPlist',help="Input meta analysis file default columns")
parser.add_argument('-token',help="LD link token, pelase request one at https://ldlink.nci.nih.gov/?tab=apiaccess",default=None)
parser.add_argument('-Nproxies',help="number of proxies per SNP",default=3,type=int)
parser.add_argument('-o','--output',help="output name ",default='proxyList')

parser._actions[0].help="Print this help message and exit"
args = parser.parse_args()

if args.token == None:
	print("Please provide an LD link token via the -token flag")
	sys.exit()
	
print("Opening SNPlist file: %s"%(args.SNPlist))
SNPlist=[]
with open(args.SNPlist,'r') as infile:
	for line in infile:
		SNPlist.append(line.strip())
SNPlist=list(set(SNPlist))
print("Read %s SNPs from input file"%(len(SNPlist)))

with open(args.output,'w') as outfile:
	for SNP in SNPlist:
		print("Requesting proxies for %s"%(SNP))
		r= requests.get('https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=%s&pop=CEU&r2_d=r2&token=%s'%(SNP,args.token))
		lines=r.text.split('\n')
		currProxy=SNP
		outfile.write("%s\n"%(currProxy))
		if re.search("error",lines[0]):
			print("Could not find SNP %s"%(SNP))
			continue
		i=0
		j=0
		while i < (args.Nproxies):
			try:
				currProxy=lines[j+2].split('\t')[0]
				j=j+1
				if currProxy!=".":
					outfile.write("%s\n"%(currProxy))
					i=i+1
			except IndexError:
				print("SNP %s had less than %s proxy SNPs"%(SNP,args.Nproxies))
				break
				
print("Finished requesting proxies for %s SNPs"%(len(SNPlist)))
	
	