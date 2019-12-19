# -*- coding:utf-8 -*-

import re
import sys
import math
import numpy as np
import matplotlib.pyplot as pplot
import ForceField
from LecturePDB import *
from partieA import *

#list([Atom])*list([Atom])*list(int)*list(int)->list(list(float))
def createMatrix(lesAtomes1,lesAtomes2,sel_p1,sel_p2):
	at1 = selgoodatoms(lesAtomes1,sel_p1)
	at2 = selgoodatoms(lesAtomes2,sel_p2)
	m = []
	for a1 in at1:
		l = []
		for a2 in at2:
			 l.append(distanceAtoms(a1,a2))
		m.append(l)
	return m


#list(list(float))*list(list(float))->float
def contactDissimilarity(matrix1,matrix2):
	n = len(matrix1)
	m = len(matrix1[0])
	S = 0
	for i in range(n):
		for j in range(m):
			S += (matrix1[i][j] - matrix2[i][j])**2
	S = sqrt(S/(n*m))
	return S

if __name__ == '__main__':
	if (len(sys.argv)<2):
		print "USAGE: ", sys.argv[0], "<pdb file1> <sel file1> <sel file2>"
		sys.exit(1)

	filename1 = sys.argv[1]
	filename2 = sys.argv[2]
	filename3 = sys.argv[3]
	_, _, _, _, lesAtomes1=readPDB(filename1, ["A"])
	sel1 = readsel(filename2)
	sel2 = readsel(filename3)

	print "Nb Atomes lus Molec.1:", len(lesAtomes1)
	print "Nb Atomes dans Select.1, Select.2:", len(sel1),",",len(sel2)

	M = createMatrix(lesAtomes1,lesAtomes1,sel1,sel2)
	print len(M),",",len(M[0])

	pplot.imshow(M)
	pplot.show()

