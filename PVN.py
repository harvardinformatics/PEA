#Created by Alexandria D'Souza

import sys
import csv
import os
import re
import numpy as np
import pandas as pd
import argparse

deltaFile=sys.argv[1]
proteinMatrix=sys.argv[2]


def readDeltamz(infile):
	DeltaPandas=pd.read_csv(infile)
	DeltaPandas['Deltam/z [Da]'] = DeltaPandas['Deltam/z [Da]'].astype(float)
	DeltaPandas=DeltaPandas.groupby(['Master Protein Accessions'])['Deltam/z [Da]'].mean()
	return DeltaPandas


def readProteinMatrix(infile):
	fn = infile
	with open(fn, 'r') as myfile:
	    file = csv.reader(myfile, dialect='excel')
	    allData=[]
	    for row in file:
	        allData.append(row)
	myfile.close()
	return allData



def writeFile(infile, DataMatrix):
	fn = infile + '_PVNtest.csv'
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    i=0
	    while i<len(DataMatrix):
	        outputFile.writerows([DataMatrix[i]])
	        i+=1


def main():
	pMatrix=readProteinMatrix(proteinMatrix)
	dFile=readDeltamz(deltaFile)

	pMatrix=(list(map(list, zip(*pMatrix))))
	groups=pMatrix[len(pMatrix)-1]
	for protein in pMatrix:
		controlGroups=[]
		treatmentGroups=[]
		i=1
		while i<len(groups):
			if group=='0':
				controlGroups.append(float(protein[i]))
			elif group=='1':
				treatmentGroups.append(float(protein[i]))
			i+=1
		dFile[protein]
		#Implement Matrix Normalization, Row and Column
		#Best method: Dict? List of Lists? NumPy Array?


	

if __name__=="__main__":
	main()

