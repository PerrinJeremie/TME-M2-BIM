# -*- coding:utf-8 -*-

import re
import sys
import math
import numpy as np
import matplotlib.pyplot as pplot
import ForceField

class Atom:
	def __init__(self, x, y, z, numRes=-1, resType="X", sse="X",bvalue=0,chain="-",atomname="CA",atomnum=-1):
		self.x = x
		self.y = y
		self.z = z
		self.numRes=numRes #numéro PDB du résidu
		self.resType=resType #code 3 lettres de l'aa
		self.sse=sse #Code 1 lettre de la Structure secondaire
		self.bvalue=bvalue
		self.chain = chain
		self.atomname=atomname
		self.atomnum=atomnum

#EXPDTA    X-RAY D
#str->str
def parse_expdta(ligne):
	re_espaceFin=re.compile("[ ]*$")
	return re_espaceFin.sub('',ligne[10:79])

#REMARK 2 RESOLUTION. 1.74 ANGSTROMS.
def parse_resol(ligne):
	resol=float(ligne[23:30])
	return resol

#NUMMDL 20
#str->int
def parse_nummdl(ligne):
	return int(ligne[10:14])

#DBREF 2JHQ A 1 226 UNP Q9KPK8 UNG_VIBCH 1    226
#str->str
def getChainDBREF(ligne):

	chain = ligne[12]
	if chain == " ":
		chain ="-"
	return chain

#ATOM   1272  CA  MET A 145      16.640  22.247  -9.383  1.00 56.65           C
#str->str*Atom ou None
def getAtom(ligne,allAtoms):
	re_whitespaces = re.compile("[ ]+")
	atomnum = int(ligne[7:11])
	resname = re_whitespaces.sub('',ligne[17:20])
	atomname = re_whitespaces.sub('',ligne[12:16])
	chain = ligne[21]
	if chain == " ":
		chain ="-"
	resnumber = int(ligne[22:26])
	x = float(ligne[30:38])
	y = float(ligne[38:46])
	z = float(ligne[46:54])
	if atomname in ['H1','H2','H3'] :
		atomname = 'H'
	if atomname == "OXT":
		atomname = 'O'
	if resname == "HIS":
		resname = "HID"
	#Autres champs
	#inscode = re_whitespaces.sub('',ligne[26])
	#atomnumber = int(ligne[6:11])
	#altloc = re_whitespaces.sub('',ligne[16])
	#occupancy = float(ligne[54:60])
	#tempfactor = float(ligne[60:66])
	#element = re_whitespaces.sub('',ligne[76:78])
	#charge = re_whitespaces.sub('',ligne[78:80])
	if (atomname=="CA" or allAtoms):
		return chain, Atom(x,y,z,resnumber, resname,chain=chain,atomname=atomname, atomnum=atomnum)
	else:
		return chain, None

#HELIX 1 HAGLYA 86 GLYA 94 1 9
#str->[int]
def getISeqNumHelix(ligne):
	lesPosH=[]
	re_whitespaces = re.compile("[ ]+")
	#initResName = re_whitespaces.sub('',ligne[15:18])
	#initICode = re_whitespaces.sub('',ligne[25])
	#endResName = re_whitespaces.sub('',ligne[27:30])
	#endICode = re_whitespaces.sub('',ligne[37])
	#Class = int(ligne[38:40])
	#Length = int(ligne[71:76])
	initChainID = re_whitespaces.sub('',ligne[19])
	initSeqNum = int(ligne[21:25])
	endChainID = re_whitespaces.sub('',ligne[31])
	endSeqNum = int(ligne[33:37])
	if initChainID != endChainID:
		print "Problème deux chaines différentes ? ", initChainID, endChainID
 		print ligne
 	else:
		for i in range(initSeqNum, endSeqNum+1):
			lesPosH.append(i)
	return lesPosH

#str->[int]
def getiSeqNumSheet(ligne):
	lesPosS=[]
	re_whitespaces = re.compile("[ ]+")
	#sheetID = re_whitespaces.sub('',ligne[11:14])
	#numStrands = int(ligne[14:16])
	#initResName = re_whitespaces.sub('',ligne[17:20])
	initChainID = re_whitespaces.sub('',ligne[21])
	initSeqNum = int(ligne[22:26])
	#initICode = re_whitespaces.sub('',ligne[26])
	endResName = re_whitespaces.sub('',ligne[28:31])
	endChainID = re_whitespaces.sub('',ligne[32])
	endSeqNum = int(ligne[33:37])
	#endICode = re_whitespaces.sub('',ligne[37])
	#sense = int(ligne[38:40])
 	if initChainID != endChainID:
 		print "Problème deux chaines différentes ? ", initChainID, endChainID
 		print ligne
 	else:
		for i in range(initSeqNum, endSeqNum+1):
			lesPosS.append(i)
	return lesPosS

#[Atom]*[int]*[int]->[Atom]
#Attention lesAtomes est modifiée
def ajouterSSEAtomes(lesAtomes, lesH, lesS):
	if(len(lesH)==0 and len(lesS)==0):
		return lesAtomes
	lesH.sort()
	lesS.sort()
	iAt=0
	iH=0
	iS=0
	#3)lesAtomes
	print
	while iAt<len(lesAtomes) and (iH<len(lesH) or iS<len(lesS)):
		if iH<len(lesH) and lesAtomes[iAt].numRes==lesH[iH]:
			lesAtomes[iAt].sse="H"
			iH+=1
		elif iS<len(lesS) and lesAtomes[iAt].numRes==lesS[iS]:
			lesAtomes[iAt].sse="S"
			iS+=1
		iAt+=1
	return lesAtomes

def getNumModel(ligne):
	return int(ligne[10:14])

#str->str*float*int*[str]*[Atom]
#exp, resol, nummdl, lesChaines, lesAtomes
def readPDB(nomFi, chaineVoulue=["-"], modelVoulu="1",allAtoms=False):
	try:
		f = open(nomFi, "r")
	except IOError:
		print "readPDB:: Fichier <%s> introuvable, arret du programme"%(nomFi)
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
	lesAtomes=[]
	for ligne in f:
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
		elif re_atom.search(ligne):
			c,a=getAtom(ligne,allAtoms)
			#print modelVoulu, numModel, c, chaineVoulue, a
			if (numModel==None):
				numModel=0
			if (a!=None and (numModel==0 or str(numModel)==modelVoulu) and c in chaineVoulue) :
				lesAtomes.append(a)
		elif re_helix.search(ligne):
			lesResHelix+=getISeqNumHelix(ligne)
		elif re_sheet.search(ligne):
			lesResHelix+=getiSeqNumSheet(ligne)
	lesAtomes=ajouterSSEAtomes(lesAtomes, lesResHelix, lesResSheet)

	return  exp, resol, nbmdl, lesChaines, lesAtomes



##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            PARTIE B :                              ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

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

'''
def main() :
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


'''
3) La RMSD calculée sur la base des deux sélection rouge et bleu est 9.17399686312.
4) Sur PyMol, on remarque que 3pdz se superpose avec une partie de la molécule 1fcf. La valeur du RMSD correspondante est RMSD =    1.958 (266 to 266 atoms).
Le résidu 21 de 3pdz et le résidu 159 de 1cfc ne sont pas superposés. On remarque alors
que les alignements sur séquence et sur structure ménent à deux résultats différents.
On en déduit que la structure tertiaire n'est pas seulement conditionnée par la séquence
d'acide aminée mais aussi par d'autre facteurs (reste de la molécule chez 1fcf par exemple).
'''
##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                              PARTIE C :                             ##
##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##


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
'''
def main() :
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
'''

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                                PARTIE D                              ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

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


'''
def main() :
	if (len(sys.argv)<2):
		print "USAGE: ", sys.argv[0], "<pdb infile> <pdb outfile>"
		sys.exit(1)

	filename1 = sys.argv[1]
	filename2 = sys.argv[2]
	_, _, _, _, lesAtomes1=readPDB(filename1, ["A","B"], allAtoms = True)
	print "Nb Atomes lus Molec.1:", len(lesAtomes1)

	_ = calculCV(lesAtomes1,20)
	writePDB(filename1,filename2,lesAtomes1,["A","B"])
'''

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                              PARTIE E                                ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

F = 332.0522

def champsForce(dist,Re,q,eps):
	epsij = math.sqrt(eps[0]*eps[1])
	Reij = Re[0]+Re[1]
	A = epsij*(Reij)**8
	B = 2*epsij*(Reij)**6
	return A/dist**8 - B/dist**6 + F*q[0]*q[1]/(20*dist)

def computeE(chainAtoms1, chainAtoms2):
	dcharge = ForceField.chargePDB()
	dvdw, depsilon = ForceField.epsilon_vdw_PDB()
	s = 0
	for a in chainAtoms1:
		epsi = depsilon[a.resType][a.atomname]
		Rie = dvdw[a.resType][a.atomname]
		qi = dcharge[a.resType][a.atomname]
		for b in chainAtoms2:
			Rij = distanceAtoms(a,b)
			epsj = depsilon[b.resType][b.atomname]
			Rje = dvdw[a.resType][a.atomname]
			qj = dcharge[a.resType][a.atomname]
			s += champsForce(Rij,[Rie,Rje],[qi,qj],[epsi,epsj])
	return s

"""def main() :
	if (len(sys.argv)<2):
		print "USAGE: ", sys.argv[0], "<pdb infile>"
		sys.exit(1)

	filename1 = sys.argv[1]
	_, _, _, _, lesAtomes=readPDB(filename1, ["A","B"], allAtoms = True)
	lesAtomesA = [a for a in lesAtomes if a.chain == "A"]
	lesAtomesB = [a for a in lesAtomes if a.chain == "B"]
	print "Nb Atomes lus chainA, chainB :", len(lesAtomesA), len(lesAtomesB)
	s = computeE(lesAtomesA,lesAtomesB)
	print(s)
	#_ = calculCV(lesAtomes1,20)
	#writePDB(filename1,filename2,lesAtomes1,["A","B"])"""

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                              PARTIE D                                ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

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


def main() :
	if (len(sys.argv)<2):
		print "USAGE: ", sys.argv[0], "<pdb infile>"
		sys.exit(1)

	filename1 = sys.argv[1]
	_, _, _, _, lesAtomes=readPDB(filename1,chaineVoulue=["A"])

	kirchhoff = kirchhoff_matrix(lesAtomes, 1, 7)

	diag, L = np.linalg.eigh(kirchhoff)
	writePDB_elastic(filename1,'output.pdb',lesAtomes,kirchhoff,L[:,1],chaineVoulue=["A"])

	#_ = calculCV(lesAtomes1,20)
	#writePDB(filename1,filename2,lesAtomes1,["A","B"])





##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                               THE END                                ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

if __name__ == '__main__':
	main()
