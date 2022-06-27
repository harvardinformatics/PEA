from scipy.optimize import curve_fit
import sys
import csv
import os
import re
import numpy as np
import argparse

begCol=int(sys.argv[3])
midCol=int(int(sys.argv[3])+(int(sys.argv[2])))
endCol=midCol+int(sys.argv[4])
reps=int(sys.argv[2])
reps2=int(sys.argv[4])
name=sys.argv[5].split('.')[0]

fileName=sys.argv[1]

def readFile(infile):
	fn = infile
	with open(fn, 'r') as myfile:
	    file = csv.reader(myfile, dialect='excel')

	    Trt = []
	    Ctl =[]
	    allData=[]
	    for row in file:
	        Trt.append(row[begCol:midCol])
	        Ctl.append(row[midCol:endCol])
	        allData.append(row)
	myfile.close()
	return Trt, Ctl, allData


def Technical(x, a, b, c, D):
	return a*(x**3)+b*(x**2)+c*x+D

def Experimental(x, a, b, c, D):
	return a*(x**3)+b*(x**2)+c*x+D 

def Detection(x, a, b, c, d, E):
	return a*(x**4)+b*(x**3)+c*(x**2)+d*x+E 

def writeFile(infile, DataMatrix):
	fn = infile + '_DET_Corrected.csv'
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    i=0
	    while i<len(DataMatrix):
	        outputFile.writerows([DataMatrix[i]])
	        i+=1


def main():
	dataTrt, dataCtl, dataMatrix =readFile(fileName)
	yT_Trt=[]
	TechnicalScore_Trt=[]
	yE_Trt=[]
	ExperimentalScore_Trt=[]
	yD_Trt=[]
	DetectionScore_Trt=[]
	yT_Ctl=[]
	TechnicalScore_Ctl=[]
	yE_Ctl=[]
	ExperimentalScore_Ctl=[]
	yD_Ctl=[]
	DetectionScore_Ctl=[]
	i=1
	while i<len(dataTrt):
		dataTrt[i]=[float(a) for a in dataTrt[i]]
		dataCtl[i]=[float(a) for a in dataCtl[i]]
		yTech_Trt=Technical(np.mean(dataTrt[i]),5.67,4.26,3.82,2)/10**9
		yExp_Trt=Experimental(np.mean(dataTrt[i]),3.5,4.2,5.67,3.4)/10**9
		yDet_Trt=Detection(np.mean(dataTrt[i]),4.7,2.3,5.4,4.6,2.3)/10**13
		yT_Trt.append(yTech_Trt)
		TechnicalScore_Trt=np.mean(dataTrt[i])/yTech_Trt
		yE_Trt.append(yExp_Trt)
		ExperimentalScore_Trt=np.mean(dataTrt[i])/yExp_Trt
		yD_Trt.append(yDet_Trt)
		DetectionScore_Trt=np.mean(dataTrt[i])/yDet_Trt
		yTech_Ctl=Technical(np.mean(dataCtl[i]),5.67,4.26,3.82,2)/10**9
		yExp_Ctl=Experimental(np.mean(dataCtl[i]),3.5,4.2,5.67,3.4)/10**9
		yDet_Ctl=Detection(np.mean(dataCtl[i]),4.7,2.3,5.4,4.6,2.3)/10**13
		yT_Ctl.append(yTech_Ctl)
		TechnicalScore_Ctl=np.mean(dataCtl[i])/yTech_Ctl
		yE_Ctl.append(yExp_Ctl)
		ExperimentalScore_Ctl=np.mean(dataCtl[i])/yExp_Ctl
		yD_Ctl.append(yDet_Ctl)
		DetectionScore_Ctl=np.mean(dataCtl[i])/yDet_Ctl
		dataTrt[i]=[(x/np.mean([TechnicalScore_Trt,ExperimentalScore_Trt,DetectionScore_Trt])) for x in dataTrt[i]]
		dataCtl[i]=[(x/np.mean([TechnicalScore_Ctl,ExperimentalScore_Ctl,DetectionScore_Ctl])) for x in dataCtl[i]]
		dataMatrix[i][begCol:midCol]=dataTrt[i]
		dataMatrix[i][midCol:endCol]=dataCtl[i]
		i+=1
	writeFile(name, dataMatrix)

	# arrayOpt, arrayCov = curve_fit(Technical, dataJm, dataJs)
	# print(arrayOpt)

if __name__=="__main__":
	main()

