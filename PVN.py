#Created by Alexandria D'Souza

import sys
import csv
import os
import re
import numpy as np
import pandas as pd
import argparse
import math

deltaFile=sys.argv[1]
proteinMatrix=sys.argv[2]


def readDeltamz(infile):
	DeltaPandas=pd.read_csv(infile)
	DeltaPandas['Deltam/z [Da]'] = DeltaPandas['Deltam/z [Da]'].astype(float)
	DeltaPandas['Master Protein Accessions','Metric']=DeltaPandas['Master Protein Accessions'].str.split('.')
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

def applyVariance(peptideArray, deltaNorm):
	return [y - deltaNorm for y in [float(x) for x in peptideArray]]

def applyNormalization(Abundance, replArray, varianceArray):
	return math.log(Abundance/abs(np.var(replArray)-np.var(varianceArray[1])), 2)

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

	protMatrix=(list(map(list, zip(*pMatrix))))
	groups=protMatrix[len(protMatrix)-1]

	dictProtein={}

	for protein in protMatrix:
		controlGroups=[]
		treatmentGroups=[]
		i=1
		while i<len(groups):
			if groups[i]=='0':
				controlGroups.append([float(protein[i]),pMatrix[i][:-1]])
			elif groups[i]=='1':
				treatmentGroups.append([float(protein[i]),pMatrix[i][:-1]])
			i+=1
		try:
			controlReps=[]
			controlPeps=[]
			for controlAbundance in controlGroups:
				controlReps.append(controlAbundance[0])
				controlPeps.append(controlAbundance[1])
			controlNormalized=[]
			i=0
			while i<len(controlPeps):
				controlNormalized.append(applyNormalization(controlReps[i], controlReps, applyVariance(controlPeps[i],dFile[protein[0]])))
				i+=1
			controlGroups=[]

			treatmentReps=[]
			treatmentPeps=[]
			for treatmentAbundance in treatmentGroups:
				treatmentReps.append(treatmentAbundance[0])
				treatmentPeps.append(treatmentAbundance[1])
			treatmentNormalized=[]
			j=0
			while j<len(treatmentPeps):
				treatmentNormalized.append(applyNormalization(treatmentReps[j], treatmentReps, applyVariance(treatmentPeps[j],dFile[protein[0]])))
				j+=1
			treatmentGroups=[]
			dictProtein[protein[0]]=[controlNormalized,treatmentNormalized]
		except KeyError as entry:
			pass

	print(dictProtein)
if __name__=="__main__":
	main()

