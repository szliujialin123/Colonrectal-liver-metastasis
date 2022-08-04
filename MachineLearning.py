for name in dir():
 if not name.startswith('_'):
      del globals()[name]
del(name)

import pandas as pd
import numpy as np
from hyperopt import hp
import matplotlib.pyplot as plt
from sklearn.metrics import plot_roc_curve
from sklearn import datasets
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score
from sklearn.metrics import RocCurveDisplay
from sklearn.metrics import auc
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from xgboost import XGBClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

train = pd.read_csv("D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/MachineLearning/train.csv")
ar = np.array(train)
rownames = ar[:,0]
colnames = train.columns[1:]
ar = ar[:,1:]
train = pd.DataFrame(ar)
train.index = rownames
train.columns = colnames
train = np.transpose(train)

#

test = pd.read_csv("D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/MachineLearning/test.csv")
ar = np.array(test)
rownames = ar[:,0]
colnames = test.columns[1:]
ar = ar[:,1:]
test = pd.DataFrame(ar)
test.index = rownames
test.columns = colnames
test = np.transpose(test)

#

x = np.array(train)
rownames = list(train.index)
y = []
for row in rownames:
    if row[0:1] == 'M':
        y.append(1)
    if row[0:1] == 'P':
        y.append(0)
y = np.array(y)

x_test = np.array(test)
rownames = list(test.index)
y_test = []
for row in rownames:
    if row[6:] == 'LM':
        y_test.append(1)
    if row[6:] == 'CRC':
        y_test.append(0)
y_test = np.array(y_test)

x_train, x_validation, y_train, y_validation = train_test_split(x, y, test_size = 0.25)

cross_validation = StratifiedKFold(n_splits = 10, shuffle = True, random_state = 100)

#################################### SVC #######################################

tuned_parameters = [
    {"kernel": ["rbf"], "gamma": [1e-3, 1e-4], "C": [1, 10, 100, 1000]},
    {"kernel": ["linear"], "C": [1, 10, 100, 1000]},
]

clf = GridSearchCV(SVC(), param_grid = tuned_parameters, scoring = 'roc_auc',
                          cv = cross_validation)
clfSVC = clf.fit(x_train, y_train)

#################################### RF ########################################

tuned_parameters = [
    {"max_features": [2, 3, 4, 5, 6, 7, 8], 
    "n_estimators": [10, 25, 50, 75, 100, 200]}
]

clf = GridSearchCV(RandomForestClassifier(), param_grid = tuned_parameters, scoring = 'roc_auc',
                          cv = cross_validation)

clfRF = clf.fit(x_train, y_train)

################################## Xgboost #####################################

tuned_parameters = {'max_depth': [3, 4, 5, 6, 8, 10, 12, 15], 
          'colsample_bytree': [0.3, 0.4, 0.5, 0.7],
          'gamma': [0.0, 0.1, 0.2 , 0.3, 0.4],
          'min_child_weight': [1, 3, 5, 7]
         }
      
clf = GridSearchCV(XGBClassifier(use_label_encoder = False), param_grid = tuned_parameters, scoring = 'roc_auc',
                          cv = cross_validation)
                          
# The use of label encoder in XGBClassifier is deprecated and will be removed in a future release. To remove this warning, do the following: 
# 1) Pass option use_label_encoder=False when constructing XGBClassifier object;     
     
clfXGB = clf.fit(x_train, y_train, eval_metric = 'rmse') # Explicitly set eval_metric if you'd like to restore the old behavior.

################################ GaussianNB ####################################

tuned_parameters = {'var_smoothing': np.logspace(0, -9, num = 100)}

clf = GridSearchCV(GaussianNB(), param_grid = tuned_parameters, scoring = 'roc_auc',
                          cv = cross_validation)
                 
clfGaussianNB = clf.fit(x_train, y_train)

################################### LDA ########################################

tuned_parameters = {
  'solver': ['lsqr'],
  'shrinkage': [0, 1, 0.01]
}

clf = GridSearchCV(LinearDiscriminantAnalysis(), param_grid = tuned_parameters, scoring = 'roc_auc',
                          cv = cross_validation)

clfLDA = clf.fit(x_train, y_train)

################################# Visualize ####################################

Allclf = {'SVC':clfSVC, 'RF':clfRF, 'XGB':clfXGB, 'GaussianNB':clfGaussianNB, 'LDA':clfLDA}
fig, ax = plt.subplots()
for i in ['SVC', 'RF', 'XGB', 'GaussianNB', 'LDA']: 
    clf = Allclf[i]  
    
    if i == 'XGB':
        clf.fit(x_train, y_train, eval_metric = 'rmse')
    
    else:
        clf.fit(x_train, y_train)
    
    viz = RocCurveDisplay.from_estimator(
        clf,
        x_test,
        y_test,
        name = "ROC {}".format(i),
        alpha = 0.3,
        lw = 1,
        ax = ax,
    )

ax.plot([0, 1], [0, 1], linestyle = "--", lw = 2, color = "r", alpha = 0.8)

plt.show()

plt.savefig("D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/MachineLearning/AllclfROC.png")


clfSVC.best_params_
clfRF.best_params_
clfXGB.best_params_
clfGaussianNB.best_params_
clfLDA.best_params_
