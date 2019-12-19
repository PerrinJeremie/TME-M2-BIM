# -*- coding:utf-8 -*-

import re
import sys
import math
import numpy as np
import matplotlib.pyplot as pplot
import ForceField
from LecturePDB import *

#str->list(int)
def readsel(nomFi):
	try:
		f = open(nomFi, "r")
	except IOError:
		print "readPDB:: Fichier <%s> introuvable, arret du programme"%(nomFi)
		sys.exit(1)
	l = []
	for ligne in f:
		i = int(ligne)
		l.append(i)
	l2 = []
	for i in range((len(l))/2):
		l2 = l2 + range(l[2*i],l[2*i+1]+1)
	return l2

#list([Atom])*list(int)->list([Atom])
def selgoodatoms(lesAtomes,sel):
	at = []
	for i in sel :
		for j in lesAtomes:
			if j.numRes == i:
				at.append(j)
				break
	return at

#[Atom]*[Atom]->float
def distanceAtoms(a1,a2):
	return math.sqrt((a1.x-a2.x)**2 + (a1.y-a2.y)**2 + (a1.z-a2.z)**2)

#list([Atom])*list([Atom])*list(int)*list(int)->float
def rmsd(lesAtomes1,lesAtomes2,sel_p1,sel_p2) :
	at1 = selgoodatoms(lesAtomes1,sel_p1)
	at2 = selgoodatoms(lesAtomes2,sel_p2)
	s = 0
	n = len(at1)
	for i in range(n) :
		x = (at1[i].x - at2[i].x)**2
		y = (at1[i].y - at2[i].y)**2
		z = (at1[i].z - at2[i].z)**2
		s += x + y + z
	s = math.sqrt(s/n)
	return s


if __name__ == '__main__':
	if (len(sys.argv)<2):
		print "USAGE: ", sys.argv[0], "<pdb file1> <pdb file2> <sel file1> <sel file2>"
		sys.exit(1)

	filename1 = sys.argv[1]
	filename2 = sys.argv[2]
	filename3 = sys.argv[3]
	filename4 = sys.argv[4]
	_, _, _, _, lesAtomes1=readPDB(filename1, ["A"])
	_, _, _, _, lesAtomes2=readPDB(filename2, ["A"])
	sel1 = readsel(filename3)
	sel2 = readsel(filename4)

	print "Nb Atomes lus Molec.1, Molec.2:", len(lesAtomes1),",",len(lesAtomes2)
	print "Nb Atomes dans Select.1, Select.2:", len(sel1),",",len(sel2)
	print "RMSD de la selection: ", rmsd(lesAtomes1,lesAtomes2,sel1,sel2)



'''
3) La RMSD calculée sur la base des deux sélection rouge et bleu est 9.17399686312.
4) Sur PyMol, on remarque que 3pdz se superpose avec une partie de la molécule 1fcf. La valeur du RMSD correspondante est RMSD =    1.958 (266 to 266 atoms).
Le résidu 21 de 3pdz et le résidu 159 de 1cfc ne sont pas superposés. On remarque alors
que les alignements sur séquence et sur structure ménent à deux résultats différents.
On en déduit que la structure tertiaire n'est pas seulement conditionnée par la séquence
d'acide aminée mais aussi par d'autre facteurs (reste de la molécule chez 1fcf par exemple).
'''
