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
replChannels=sys.argv[3]

def readDeltamz(infile):
	DeltaPandas=pd.read_csv(infile)
	DeltaPandas['Deltam/z [Da]'] = DeltaPandas['Deltam/z [Da]'].astype(float)
	DeltaPandas['Master Protein Accessions','Metric']=DeltaPandas['Master Protein Accessions'].str.split('-')
	DeltaPandas=DeltaPandas.groupby(['Master Protein Accessions'])['Deltam/z [Da]'].mean()
	return DeltaPandas

def readReplicatesChannels(infile):
	fn = infile
	with open(fn, 'r') as myfile:
	    file = csv.reader(myfile)
	    channel=[]
	    group=[]
	    for row in file:
	    	channel.append(row[0])
	    	group.append(row[1])
	myfile.close()
	return channel,group

def readProteinMatrix(infile):
	fn = infile
	with open(fn, 'r') as myfile:
	    file = csv.reader(myfile, dialect='excel')
	    allData=[]
	    for row in file:
	    	allData.append(row[:-1])
	myfile.close()
	return allData

def applyVariance(peptideArray, deltaNorm):
	return [y - deltaNorm for y in [float(x) for x in peptideArray]]

def applyNormalization(Abundance, replArray, varianceArray):
	return abs(math.log(Abundance/abs(np.var(replArray)-np.var(varianceArray[1])), 2))

def writeFile(DataMatrix, group, Channel):
	fn = 'PVN_ProteinNormalization_matrix_formatting.csv'
	#Need to add Channel Description
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    for key, value in DataMatrix.items():
	    	outputFile.writerows([[key]+value[0]+value[1]])
	    outputFile.writerows([['TreatmentGroups']+group])
	    outputFile.writerows([Channel])

def writeFilePMatrixTranspose(DataMatrix, group, Channel):
	fn = 'ProteinNormalization_matrix_transpose.csv'
	#Need to add Channel Description
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    i=0
	    c=1
	    t=1
	    ctTranspose=[]
	    while i<len(group):
	    	if group[i]=='0':
	    		ctTranspose.append('Control'+'_'+str(c))
	    		c+=1
	    	elif group[i]=='1':
	    		ctTranspose.append('Treatment'+'_'+str(t))
	    		t+=1
	    	i+=1
	    outputFile.writerows([['Protein']+ctTranspose])
	    for key, value in DataMatrix.items():
	    	outputFile.writerows([[key]+value[0]+value[1]])


def writeFilePNormMatrix(DataMatrix, group, Channel):
	fn = 'ProteinNormalization_matrix.csv'
	#Need to add Channel Description
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)

	    #outputFile.writerows([['TreatmentGroup']+['Channel']+list(DataMatrix.keys())])
	    transposeDict=[]
	    for key, value in DataMatrix.items():
	    	transposeDict.append([key]+value[0]+value[1])
	    transposeDict=(list(map(list, zip(*transposeDict))))

	    i=0
	    while i<len(transposeDict):
	    	outputFile.writerows([[group[i]]+[Channel[i]]+transposeDict[i]])
	    	i+=1

def writeFilePNormMatrixBox(DataMatrix, group, Channel):
	fn = 'ProteinNormalization_matrix_box.csv'
	#Need to add Channel Description
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    transposeDict=[]
	    for key, value in DataMatrix.items():
	    	transposeDict.append([key]+value[0]+value[1])
	    transposeDict=(list(map(list, zip(*transposeDict))))

	    i=0
	    while i<len(transposeDict):
	    	if group[i]=='0':
	    		outputFile.writerows([['Control']+[Channel[i]]+transposeDict[i]])
	    	elif group[i]=='1':
	    		outputFile.writerows([['Treatment']+[Channel[i]]+transposeDict[i]])
	    	i+=1


def main():
	pMatrix=readProteinMatrix(proteinMatrix)
	dFile=readDeltamz(deltaFile)
	pChannel, groups=readReplicatesChannels(replChannels)

	protMatrix=(list(map(list, zip(*pMatrix))))
	protMatrix=[groups]+[pChannel]+protMatrix


	dictProtein={}

	for protein in protMatrix:
		if all([protein[0] != 'TreatmentGroup', protein[0] != 'Channel']):
			#print(protein)
			controlGroups=[]
			treatmentGroups=[]
			formatKey=protein[0].replace('.','-')
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
					controlNormalized.append(applyNormalization(controlReps[i], controlReps, applyVariance(controlPeps[i],dFile[formatKey])))
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
					treatmentNormalized.append(applyNormalization(treatmentReps[j], treatmentReps, applyVariance(treatmentPeps[j],dFile[formatKey])))
					j+=1
				treatmentGroups=[]
				dictProtein[formatKey]=[controlNormalized,treatmentNormalized]
			except KeyError as entry:
				pass
		writeFile(dictProtein, groups[1:], pChannel)
		writeFilePMatrixTranspose(dictProtein, groups[1:], pChannel)
		writeFilePNormMatrix(dictProtein, groups, pChannel)
		writeFilePNormMatrixBox(dictProtein, groups, pChannel)

if __name__=="__main__":
	main()

