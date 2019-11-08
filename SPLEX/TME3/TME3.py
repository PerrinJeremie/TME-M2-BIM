import pandas as pd
import graphviz
from sklearn import tree
from sklearn.ensemble import RandomForestClassifier

data_diabetes = pd.read_table(’patients_data.txt’,sep=’\t’,header=None)
classes_diabetes = pd.read_table(’patients_classes.txt’,sep=’\t’,header=None)
