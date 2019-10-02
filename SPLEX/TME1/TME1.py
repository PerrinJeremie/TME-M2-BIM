#!/usr/bin/python3
# -*- coding: latin-1 -*-
				##############################
                                ##       TME1 - SPLEX       ##
				##############################

import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import statsmodels.sandbox.stats.multicomp as sm
import math

#########––––––––––––––––––––––––––––––––––––###########
##                IMPORTING DATA                      ##
#########––––––––––––––––––––––––––––––––––––###########

breast_data = pd.read_table("BreastDiagnostic.txt", sep =",")
mice_data = pd.ExcelFile("Data_Cortex_Nuclear.xls")
mice_data = mice_data.parse()
mice_data.fillna(method = "bfill", inplace=True)
mice_data.fillna(0, inplace=True)


################––––––––––––––––––––––––#################
##            PLOTTING CORRELATION MATRIX              ##
################––––––––––––––––––––––––#################

classes = breast_data.values[:,1]
bdf = breast_data.values[:,2:]
nobs,nmeas = bdf.shape
m = np.zeros((nmeas,nmeas,2))
for i in range(nmeas):
	for j in range(i,nmeas):
		(c,p) = stats.pearsonr(bdf[:,i],bdf[:,j])
		m[i,j,0] = c
		m[i,j,1] = p


'''
plt.imshow(m[:,:,0])
plt.colorbar()
plt.show()
'''

#############–––––––––––––––––––––––––––––##############
##              MANN-WHITNEY-WILCOXON                 ##
#############–––––––––––––––––––––––––––––##############

wilcc = np.zeros((nmeas))
Mft = []
Bft = []
for i in range(nobs):
	if classes[i] == 'M':
		Mft.append(bdf[i,:])
	else:
		Bft.append(bdf[i,:])
Mft = np.array(Mft)
Bft = np.array(Bft)

for k in range(nmeas):
	wilcc[k] = stats.mannwhitneyu(Mft[:,k],Bft[:,k],alternative='two-sided')[1]

print(wilcc)

#########–––––––––––––––––––––––––––––––––––––##########
##                    MULTIPLE TEST                   ##
#########–––––––––––––––––––––––––––––––––––––##########

'''
There is a choice to be made between two p-values correction techniques:
	- FDR correction
	- FWER correction
	The most stringent method being FWER and the least FDR. FDR is less stringent because on average it rejects a higher number of true null hypothesis but at the same time this allows it to capture more false null hypothesis.
We also note that controlling the FWER also corrects FDR, so it is more stringent.
'''
wilcc_corrected = sm.multipletests(wilcc,method="holm-sidak")
print(wilcc_corrected)


##############–––––––––––––––––––––––––––###############
##              COMPARING DISTRIBUTIONS               ##
##############–––––––––––––––––––––––––––###############
print(Mft)
Mft = Mft.astype(np.float) 
Bft = Bft.astype(np.float) 
st, pval = stats.ttest_ind(Mft,Bft)

print(sm.multipletests(pval,method="holm-sidak"))

plt.boxplot(Mft)
plt.boxplot(Bft)
plt.show()

'''


'''



#######––––––––––––––––––––––––––––––––––––––––––#######
##                       ADABOOST                     ##
#######––––––––––––––––––––––––––––––––––––––––––#######

'''
Ca ca a l'air chouette, il faut s'amuser
'''


'''
mdf = mice_data.values[:,1:-4]
nobs,nmeas = mdf.shape


print(mdf[0])

m = np.zeros((nmeas,nmeas,2))
for i in range(nmeas):
	for j in range(i,nmeas):
		(c,p) = stats.pearsonr(mdf[:,i],mdf[:,j])
		m[i,j,0] = c
		m[i,j,1] = p

plt.imshow(m[:,:,0])
plt.colorbar()
plt.show()
'''




