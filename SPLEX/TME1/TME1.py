#!/usr/bin/python3
# -*- coding: latin-1 -*-
				#######################################
                                ##           TME1 - SPLEX            ##
				## PERRIN Jérémie - PODLEJSKI Witold ##
				#######################################

import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt
import statsmodels.sandbox.stats.multicomp as sm
from sklearn.ensemble import AdaBoostClassifier
from sklearn.datasets import make_classification
import math

#########––––––––––––––––––––––––––––––––––––###########
##                IMPORTING DATA                      ##
#########––––––––––––––––––––––––––––––––––––###########

breast_data = pd.read_table("BreastDiagnostic.txt", sep =",")
mice_data = pd.ExcelFile("Data_Cortex_Nuclear.xls")
mice_data = mice_data.parse()
#mice_data.fillna(method = "bfill", inplace=True)
mice_data.fillna(mice_data.median(), inplace=True)


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

plt.imshow(m[:,:,0])
plt.colorbar()
plt.show()

'''
We find that some groups of variable are strongly positively correlated 
But very few are strongly negatively correlated.
A positive correlation means a change in one variable implies a change in the other variable in the same direction. A negative correlation means a change in one variable 
is often linked to a change of the other variable in the other direction.
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

print("##############")
print("## WILCOXON ##")
print("##############")
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
print("########################")
print("## WILCOXON CORRECTED ##")
print("########################")
print(wilcc_corrected[1])


##############–––––––––––––––––––––––––––###############
##              COMPARING DISTRIBUTIONS               ##
##############–––––––––––––––––––––––––––###############


Mft = Mft.astype(np.float) 
Bft = Bft.astype(np.float) 
st, pval = stats.ttest_ind(Mft,Bft)

print("###############")
print("## TTEST_IND ##")
print("###############")
print(sm.multipletests(pval,method="holm-sidak")[1])

	

plt.figure()
plt.boxplot([Mft[:,0],Bft[:,0]])
plt.title("Discriminating Variable")
plt.figure()
plt.boxplot([Mft[:,9],Bft[:,9]])
plt.title("Non Discriminating Variable")
plt.show()

'''
Some variable distributions are indistinguishable in between the two classes (col 9)
Other have different distributions (col 0)
Which is coherent with Student's T-Test.
'''

######––––––––––––––––––––––––––––––––––######
##               CONCLUSION                 ##
######––––––––––––––––––––––––––––––––––######

'''
We want to select variables (to use for example in the next task), to do this we test 
which variables we should use. TTEST_IND and WILCOXON both indicate some variables are not significant.
The adjustment method really depends on the costs of being too stringent or too laxist.
'''


#######–––––––––––––––––––––––––––––––––––––––––––––––––#######
##                       ADABOOST - TEST                     ##
#######–––––––––––––––––––––––––––––––––––––––––––––––––#######


classifier = AdaBoostClassifier(n_estimators=100, random_state=0)

X = np.concatenate((np.random.normal(loc=12, size=(500,1)) , np.random.normal(loc=10, size=(500,1))), axis=0)
y = np.concatenate((np.ones((500)), np.zeros((500))))


#######–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#######
##                       ADABOOST - Weak Learner Numbers                     ##
#######–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#######



t = []
for i in range(1,100,10):
	classifier = AdaBoostClassifier(n_estimators=i, random_state=0)
	classifier.fit(X,y)
	t.append(classifier.score(X,y))
fig = plt.figure()
plt.xlabel('Number of Weak Learners', fontsize=18)
plt.ylabel('Score of Classifier', fontsize=16)
plt.plot(t)
plt.show()

#######–––––––––––––––––––––––––––––––––––––––––––––––––––––––#######
##             ADABOOST - 2D with breast cancer data               ##
#######–––––––––––––––––––––––––––––––––––––––––––––––––––––––#######

classifier = AdaBoostClassifier(n_estimators=100, random_state=0)
X = np.array([[bdf[i,0],bdf[i,1]] for i in range(nobs)])
y = np.array([ (0 if a == 'M' else 1) for a in classes])


classifier.fit(X,y)

fig = plt.figure()
plot_colors = "br"
class_names = "AB"

# Plot the decision boundaries
x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
plot_step_x = (x_max - x_min)/100
plot_step_y = (y_max - y_min)/100
xx, yy = np.meshgrid(np.arange(x_min, x_max, plot_step_x),
                     np.arange(y_min, y_max, plot_step_y))

Z = classifier.predict(np.c_[xx.ravel(), yy.ravel()])
Z = Z.reshape(xx.shape)
cs = plt.contourf(xx, yy, Z, cmap=plt.cm.Paired)
plt.axis("tight")

# Plot the data points
for i, n, c in zip(range(2), class_names, plot_colors):
    idx = np.where(y == i)
    plt.scatter(X[idx, 0], X[idx, 1],
                c=c, cmap=plt.cm.Paired,
                s=20, edgecolor='k',
                label="Class %s" % n)

fig.suptitle('Classifying based on columns 0 and 1', fontsize=20)
plt.xlabel('Column 0', fontsize=18)
plt.ylabel('Column 1', fontsize=16)
plt.show()

'''
Adaboost is a non-linear classifier.
'''
