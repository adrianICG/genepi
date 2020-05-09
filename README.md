# genepi
Bancroft level 7 collaborative repository.

This repository is for people in QIMR particularly Bancroft level 7 to make available standard or example scripts for performing statistical genetics analyses.

It currently contains:
* PRS pipeline for genepi Genotyped datasets
* Some GWAS scripts for genepi genotyped datasets
* HPC utilities for people to develop their own scripts using python
* A common factor GWAS wrapper for GenomicSEM

For each folder please check available documentation. You can clone this repository into the HPC and have all of this functions available in your bin.

Intructions for how to set this up working for you (Follow this steps to ensure everything runs smoothly!):

1. Navigate to your home
```bash
cd $HOME
```

2. If you dont already have a bin directory, create one (this is basically where Unix will look for the scripts when you call them)
```bash
mkdir bin
```

3. Append this new directory to your PATH variable by editing the file .bashrc:
```bash
nano ~/.bashrc
```
Add the following lines at the end of the file:
```
export PATH="$PATH:/$HOME/bin"
export PYTHONPATH="$HOME/bin/genepi/HPCfunctions/"
```
Exit (ctrl+x) and save changes. Setting Python path will be of great importance for python to find the QIMRBHPC utility functions that are common to many scripts. (They will all eventually import those functions)

You will need to restart your session for these changes to take effect.

4. Clone the github repository somewhere in your HOME (in this case ill show you inside bin):
```bash
cd $HOME/bin/
git clone https://github.com/adrianICG/genepi.git
```
You will get prompted for username and password (github credentials)

5. Create simbolic links for each of the scripts to add them to your bin (they have to be visible right in the bin, and with executable properties; double check that).
```bash
cd $HOME/bin/
ln -s genepi/HPCfunctions/memoryENFORCER.py ./
ln -s genepi/PRSpipeline/plotPRS.py ./
ln -s genepi/PRSpipeline/PRSlmm.py ./
ln -s genepi/PRSpipeline/runPRSv2.py ./
ln -s genepi/CommonFactorGWAS/RunCommonFactorGWAS.py ./
ln -s genepi/GWAS/runR10SAIGEGWAS.py ./
ln -s genepi/GWAS/runAGDSGWAS.py ./
```

Whenever there is an issue fixed or an updated version you can easily update your packages by running:
```
cd $HOME/bin/genepi
git pull
```
Do not modify things unless you are working on fixing a bug or issue. In that case please make sure you read git and githubs documentation so you are aware of how to contribute