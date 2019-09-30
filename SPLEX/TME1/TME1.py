				##############################
                                ##       TME1 - SPLEX       ##
				##############################

import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import statsmodels.sandbox.stats.multicomp as sm

breast_data = pd.read_table("BreastDiagnostic.txt", sep =",")
mice_data = pd.ExcelFile("Data_Cortex_Nuclear.xls")
mice_data = mice_data.parse()
mice_data.fillna(method = "bfill", inplace=True)
mice_data.fillna(0, inplace=True)


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

wilcc = np.zeros((nmeas))
for k in range(nmeas):
	Mft = []
	Bft = []

	for i in range(nobs):
		if classes[i] == 'M':
			Mft.append(bdf[i,k])
		else:
			Bft.append(bdf[i,k])

	Mft = np.array(Mft)
	Bft = np.array(Bft)

	wilcc[k] = stats.wilcoxon(Mft,Bft,)[1]

print(wilcc)

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




