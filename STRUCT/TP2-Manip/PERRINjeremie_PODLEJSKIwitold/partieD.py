# -*- coding:utf-8 -*-

import re
import sys
import math
import numpy as np
import matplotlib.pyplot as pplot
import ForceField
from LecturePDB import *
from partieA import *

def kirchhoff_matrix(CA_list, force, cut_distance):

	matrix = np.zeros((len(CA_list),len(CA_list)))
	for index_ref, alpha_carbon_ref in enumerate(CA_list) :
		 for index, alpha_carbon in enumerate(CA_list) :
		  	if index != index_ref :
				if distanceAtoms(alpha_carbon_ref,alpha_carbon)< cut_distance :
					matrix[index_ref, index] = -force
	for i in range(len(CA_list)):
		matrix[i,i] = - np.sum(matrix[i,:])
	return matrix

def writePDB_elastic(nomFiIn,nomFiOut,lesAtomes, kirchhoff,L,chaineVoulue=["A"],modelVoulu='1',allAtoms=False,ampl = 7):
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
	re_end = re.compile("^END")
	re_conect = re.compile("^CONECT")
	re_het = re.compile("^HET")
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
			continue
		else:
			if re_end.search(ligne) or re_conect.search(ligne) or re_het.search(ligne):
				continue
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
	numModel = 1
	for i in np.arange(-np.pi/2,np.pi/2+0.01,np.pi/10):
		g.write("MODEL " + "    " + str(numModel).rjust(4) + "\n")
		l = mod_atoms(lesAtomes,L,i,ampl)
		write_atomlist(l,g)
		g.write("ENDMDL\n")
		numModel += 1
	write_CONECT(kirchhoff, lesAtomes,g)

	
def mod_atoms(lesAtomes,L,t,ampl):
	l = []
	for i in range(len(lesAtomes)):
		a = lesAtomes[i]
		f = L[i]*ampl*np.sin(t)
		p = Atom( a.x + f, a.y + f, a.z + f, a.numRes, a.resType, a.sse,a.bvalue,a.chain,a.atomname,a.atomnum)
		l.append(p)
	return l

def write_oneatom(a,file):
	file.write("ATOM  ")
	file.write(str(a.atomnum).rjust(5))
	file.write("  CA  ")
	file.write(a.resType)
	file.write(" A")
	file.write(str(a.numRes).rjust(4))
	file.write("".rjust(4))
	file.write("{:8.3f}".format(a.x))
	file.write("{:8.3f}".format(a.y))
	file.write("{:8.3f}".format(a.z))
	file.write("\n")
	
	

def write_atomlist(atomes,file):
	for a in atomes:
		write_oneatom(a,file)
	

def write_CONECT(kirchhoff, lesAtomes, file):

	for index_ref in range(len(lesAtomes)):
		write_count = 0
		for index in range(index_ref + 1,len(lesAtomes)):
			if kirchhoff[index_ref,index] != 0 :
				if write_count % 4 == 0:
					file.write("\nCONECT")
					file.write(str(lesAtomes[index_ref].atomnum).rjust(5))
				file.write(str(lesAtomes[index].atomnum).rjust(5))
				write_count += 1

	file.write('\n')


if __name__ == '__main__':
	if (len(sys.argv)<2):
		print "USAGE: ", sys.argv[0], "<pdb infile> <pdb outfile> <cut_off distance> <Amplitude>"
		sys.exit(1)

	filename1 = sys.argv[1]
	outfilename = sys.argv[2]
	coff = float(sys.argv[3])
	ampl = float(sys.argv[4])
	_, _, _, _, lesAtomes=readPDB(filename1,chaineVoulue=["A"])

	kirchhoff = kirchhoff_matrix(lesAtomes, 1, coff)

	diag, L = np.linalg.eigh(kirchhoff)
	for i in range(4):
		writePDB_elastic(filename1,"mode"+str(i)+"_"+outfilename,lesAtomes,kirchhoff,L[:,i],chaineVoulue=["A"],ampl=ampl)


