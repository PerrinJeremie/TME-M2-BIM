#!/usr/bin/python3
# -*- coding: latin-1 -*-

import matplotlib.pyplot as plt
from sklearn import cluster
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.cluster import AgglomerativeClustering
from sklearn.datasets import make_classification
from sklearn.datasets import make_blobs
from sklearn.datasets import make_moons

path = "../TME1/"

#########–––––––––––––––––––––––––––––##########
##       PRELIMINARIES - Data Creation        ##
#########–––––––––––––––––––––––––––––##########

# First simulated data set
#plt.title("Two informative features, one cluster per class", fontsize='small')
X1, Y1 = make_classification(n_samples=200, n_features=2, n_redundant=0, n_informative=2,n_clusters_per_class=1)
#plt.scatter(X1[:, 0], X1[:, 1], marker='o', c=Y1,s=25, edgecolor='k')
#plt.show()

# Second simulated data set
#plt.title("Three blobs", fontsize='small')
X2, Y2 = make_blobs(n_samples=200, n_features=2, centers=3)
#plt.scatter(X2[:, 0], X2[:, 1], marker='o', c=Y2, s=25, edgecolor='k')
#plt.show()

# Third simulated data set
#plt.title("Non-linearly separated data sets", fontsize='small')
X3, Y3 = make_moons(n_samples=200, shuffle=True, noise=None, random_state=None)
#plt.scatter(X3[:, 0], X3[:, 1], marker='o', c=Y3, s=25, edgecolor='k')
#plt.show()

#######–––––––––––––––––––––––––––––––––––######
##                   K-MEANS                  ##
#######–––––––––––––––––––––––––––––––––––######

km2 = KMeans(n_clusters=2, init='k-means++', max_iter=100, n_init=1)
km3 = KMeans(n_clusters=2, init='k-means++', max_iter=100, n_init=1)

km2.fit(X1)
plt.scatter(X1[:, 0], X1[:, 1], s=10, c=km2.labels_)
plt.show()


km3.fit(X2)
plt.scatter(X2[:, 0], X2[:, 1], s=10, c=km3.labels_)
plt.show()

km2.fit(X3)
plt.scatter(X3[:, 0], X3[:, 1], s=10, c=km2.labels_)
plt.show()

##############––––––––––––––––––––###############
##          HIERARCHICAL CLUSTERING            ##
##############––––––––––––––––––––###############

for linkage in (’ward’, ’average’, ’complete’):
	clustering2 = AgglomerativeClustering(linkage=linkage, n_clusters=2)
	clustering3 = AgglomerativeClustering(linkage=linkage, n_clusters=3)
	clustering2.fit(X1)
	plt.scatter(X2[:, 0], X2[:, 1], s=10, c=km3.labels_)
plt.show(
	clustering3.fit(X2)
	clustering2.fit(X3)




##############–––––––––––––––––––––###############
##             SPECTRAL CLUSTERING              ##
##############–––––––––––––––––––––###############


