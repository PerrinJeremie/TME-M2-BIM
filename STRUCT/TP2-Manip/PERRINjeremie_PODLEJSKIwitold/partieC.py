# -*- coding:utf-8 -*-

import re
import sys
import math
import numpy as np
import matplotlib.pyplot as pplot
import ForceField
from LecturePDB import *


def norme(l):
	return math.sqrt(l[0]**2+l[1]**2+l[2]**2)


def calculCV(lesAtomes,rc):
	n = len(lesAtomes)
	l = [[0,0,0,0] for _ in range(n)]
	res = {}
	for i in range(n):
		ai = lesAtomes[i]
		for j in range(i+1,n):
			if i != j :
				aj = lesAtomes[j]
				x = aj.x - ai.x
				y = aj.y - ai.y
				z = aj.z - ai.z
				s = math.sqrt(x**2+y**2+z**2)
				if s <= rc :
					l[j][0] += x/s
					l[j][1] += y/s
					l[j][2] += z/s
					l[j][3] += 1
					l[i][0] -= x/s
					l[i][1] -= y/s
					l[i][2] -= z/s
					l[i][3] += 1
	for i in range(n):
		if l[i][3] == 0 :
			l[i][3] = 1
		lesAtomes[i].bvalue = 1 - norme(l[i])/l[i][3]
		r = lesAtomes[i].numRes
		if res.has_key(r):
			res[r] = [res[r][0] + lesAtomes[i].bvalue,res[r][1]+1]
		else:
			res[r] = [lesAtomes[i].bvalue,1]

	for k in res:
		res[k] = res[k][0]/res[k][1]
	return res


def protuberants_enfouis(lesAtomes,x,rc):
	resCV = calculCV(lesAtomes,rc)
	resList = sorted(resCV, key=resCV.__getitem__)
	n = len(resList)
	return resList[0:int(x*n/100)],resList[int(n*(100-x)/100):]

def writeAtom(ligne,outFile,atom):
	#Autres champs
	#inscode = re_whitespaces.sub('',ligne[26])
	#atomnumber = int(ligne[6:11])
	#altloc = re_whitespaces.sub('',ligne[16])
	#occupancy = float(ligne[54:60])
	#tempfactor = float(ligne[60:66])
	#element = re_whitespaces.sub('',ligne[76:78])
	#charge = re_whitespaces.sub('',ligne[78:80])
	outFile.write(ligne[:60] + "{:6.2f}".format(atom.bvalue) + ligne[66:])

def writePDB(nomFiIn,nomFiOut,lesAtomes,chaineVoulue=["A"],modelVoulu='1',allAtoms=True):
	try:
		f = open(nomFiIn, "r")
	except IOError:
		print "writePDB:: Fichier <%s> introuvable, arret du programme"%(nomFiIn)
		sys.exit(1)

	try:
		g = open(nomFiOut, "w")
	except IOError:
		print "writePDB:: Fichier <%s> introuvable, arret du programme"%(nomFiOut)
		sys.exit(1)
	#EXPDTA peut être sur plusieurs lignes
	#Pour avoir la première, j'ajoute des espaces pour les champs continuation (colonnes 9 - 10)
	re_expdta = re.compile("^EXPDTA {4}")
	#Format: REMARK puis 3 espaces puis 2 puis un espace puis RESOLUTION et ANGSTROM à la fin
	re_resolution = re.compile("^REMARK {3}2 RESOLUTION.*ANGSTROMS")
	re_atom = re.compile("^ATOM")
	re_dbref = re.compile("^DBREF")
	re_nummdl = re.compile("^NUMMDL")
	re_helix = re.compile("^HELIX")
	re_sheet = re.compile("^SHEET")
	re_model = re.compile("^MODEL")
	exp=None
	resol=None
	nbmdl=None
	numModel=None
	lesResSheet=[]
	lesResHelix=[]
	lesChaines=[]
	i = 0
	for ligne in f:
		if re_atom.search(ligne):
			c,a=getAtom(ligne,allAtoms)
			#print modelVoulu, numModel, c, chaineVoulue, a
			if (numModel==None):
				numModel=0
			if (a!=None and (numModel==0 or str(numModel)==modelVoulu) and c in chaineVoulue) :
				writeAtom(ligne,g,lesAtomes[i])
				i+=1
		else:
			g.write(ligne)
			if re_expdta.search(ligne):
				exp=parse_expdta(ligne)
			elif re_resolution.search(ligne):
				resol=parse_resol(ligne)
			elif re_nummdl.search(ligne):
				nbmdl=parse_nummdl(ligne)
			elif re_dbref.search(ligne):
				lesChaines.append(getChainDBREF(ligne))
			elif re_model.search(ligne):
				numModel=getNumModel(ligne)
			elif re_helix.search(ligne):
				lesResHelix+=getISeqNumHelix(ligne)
			elif re_sheet.search(ligne):
				lesResHelix+=getiSeqNumSheet(ligne)


if __name__ == '__main__':
	if (len(sys.argv)<2):
		print "USAGE: ", sys.argv[0], "<pdb infile> <pdb outfile> <chains>"
		sys.exit(1)

	filename1 = sys.argv[1]
	filename2 = sys.argv[2]
	chains = sys.argv[3]
	_, _, _, _, lesAtomes1=readPDB(filename1, list(chains), allAtoms = True)
	print "Nb Atomes lus Molec.1:", len(lesAtomes1)

	_ = calculCV(lesAtomes1,20)
	print "Done for chains :",chains
	writePDB(filename1,filename2,lesAtomes1,list(chains))

