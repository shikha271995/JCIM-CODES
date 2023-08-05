
##Data collection and Data Analysis

## Import libraries
import pandas as pd
import numpy  as np

x=pd.read_csv("dpp4.csv" , sep =";")
printx.head()
x=x[x["Smiles"] .notna()]
x.shape
x=x[x["Standard Value"] .notna()]
x.shape
x=x[x["Standard Units"] .notna()]
x.shape
x=x[['Molecule ChEMBL ID','Smiles', 'Standard Type','StandardValue','Standard Units','Document Year','Target Name']]
x.head()
x=x[x["Standard Units"] .str.contains("nM")]
x.head()
x.shape
x["Standard Units"] .value_counts()
y=x["Molecule ChEMBL ID"] .value_counts()
y
y.to_csv("duplicates.csv" , sep =",", index =None , header =None )
y
y.to_csv("dupilcates")
x1=x[x["Molecule ChEMBL ID"] .str.contains("CHEMBL1422")]
x1.sort_values("Standard Value",ascending = True )
x.head
x["Molecule ChEMBL ID"] .value_counts()
x["new_value"]=x[["Molecule ChEMBL ID","Smiles", "Standard Type","Standard Value","Standard Units","Document Year","Target Name"]]].groupby (["Molecule Chembl Id"])["Standard Value"].transform("mean")
x=x.sort_values("new_value" , ascending = True )
x.head
x=x.drop_duplicates("Molecule ChEMBL ID",keep ="first")
x.shape
x["Standard Type"] .value_counts()
x.to_csv("Dpp4 results of analysis")


Molecular Data analysis , Train and test split  , Tanimoto

import pandas as pd
import numpy as np


#  x=pd.read_csv("Dpp4 results of analysis.csv")
#active_RS=x.loc[x["Standard Value"]<=10000]
inactive_RS=x.loc[x["Standard Value"]>10000]

#active_RS["label"]=1
inactive_RS["label"]=0
#combined=pd.concat([active_RS,inactive_RS], axis=0)

#In[ ]:combined[["Smiles","label"]].to_csv("dpp4_binary_labelled_model_.s mi",  sep=" \t", header=None, index=None)

#dataset = combined[["Smiles", "label"]]
#dataset.columns
# import molvs
from molvs import Standardizer
from molvs import standardize_smiles

# x1=combined[["Smiles","Molecule ChEMBL  ID","label"]]

#  x1.shape
#  import rdkit
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import rdBase
rdBase.rdkitVersion
# Draw.MolsToGridImage([Chem.MolFromSmiles(smi) for smi in  x1.iloc[0:, 0]])

# [standardize_smiles(smi) for smi in x1.iloc[0:, 0]]
#allsmiles=[standardize_smiles(smi) for smi in x1.iloc[0:, 0]]
#allsmiles=pd.DataFrame(allsmiles)
# In[ ] :x2=pd.concat([allsmiles,x1[["Molecule ChEMBL ID","label"]]], axis=1)
#x2.columns="Smiles","ChEMBL ID","label"
#x2.to_csv("data2new.smi" , sep=" \t" , index=None, header=None )
#mol2=[Chem.MolFromSmiles(smi) for smi in allsmiles.iloc[0: , 0]]
# from rdkit.Chem import SaltRemover

remover = SaltRemover.SaltRemover(defnData="[Cl,Br]")
len(remover.salts)

#mol2=[remover.StripMol(mol) for mol in mol2]
#import rdkit.Chem as chem

#In[ ]:Draw.MolsToGr idImage([chem.MolFromSmiles(smi) for smi in x2["Smiles"].iloc[0:4821]])

# dataset.to_csv("datanew.smi", sep=" \t", index=None , header=None)
# from rdkit.Chem import Scaffolds

from rdkit.Chem import AllChem
import  rdkit.Chem as Chem

# from rdkit.Chem.Scaffolds import MurckoScaffold
t1=Chem.S milesMolSupplier("datanew.smi",  delimiter=" \t",titleLine=False)

# [MurckoScaffold.GetScaffoldForMol(mol)for mol in t1 if mol is not None]
#  mol2=[MurckoScaffold.GetScaffoldForMol(mol)for mol in t1 if mol is not  None ]
#  mol2=[Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol)) for mol in mol2 if mol2 is not None]
#  mol2
#x2
#ms=pd.DataFrame(mol2)
ms.head()
#ms.shape
#combined=pd.concat([x2, ms], axis=1)
#combined.head()
#combined.shape
#combined.columns="Smiles","ChEMBL ID","label","Murcko"
#    combined[["Murcko","ChEMBL ID"]].to_csv("data_murcko.smi" ,  sep=" \t", header=None ,index=None)
#
get_ipython().system('python/home/hemant/Downloads/NOOR_PYTHON/mayachemtools/maya/bin/RDKitClusterMolecules.py —infileParams smilesColumn,1,smilesNameColumn,2,smilesDelimiter,tab,smilesTitleLine,auto,sanitize, yes -i data_murcko.smi -o clusterfinal2.smi')

#import pandas as pd
#clus=pd.read_csv("clusterfinal2.smi" , sep=" ")
#clus
#clus["ClusterNumber"].value_counts()
#Draw.MolsToGridImage([Chem.MolFromSmiles(smi) for smi in clus.iloc[0:, 0] if smi is not None])

#clus
#X_train=clus.groupby("ClusterNumber").filter(lambda x: len(x) >=2)
#X_train.shape
#X_train.head()
#X_test=c lus.groupby("ClusterNumber").filter(lambda x: len(x) <2)
#X_test.shape
#X_test.head()
#combined.head()
#combined=combined.rename(columns={"ChEMBL ID":"Name"})
#X_train_ = pd.merge(X_train, combined, on="Name")
#X_train_
#X_train_=X_train_[["Smiles","Name","label"]]
#X_train_.to_csv("xtrain.smi", sep=" \t", index=None, header=None)
#X_train_.to_csv("xtrain.csv")
#X_test_ =pd.merge(X_test, combined, on="Name")
#X _test_
#X_test_[["Smiles","Name"]]
#X_test_.to_csv("xtest.smi", sep=" \t", index=None, header=None)
#X_test_.to_csv("xtest.csv")
#X_test_=X_test_[["Smiles","Name"]]
# Draw.MolsToGridImage([Chem.MolFromSmiles(smi) for smi in X_train_.iloc[0:, 0]])
# Xtrain=Chem.SmilesMolSupplier("xtrain.smi", delimiter=" \t", titleLine=False)
#   Xtest=Chem.SmilesMolSupplier("xtest.smi", delimiter=" \t", titleLine=False)
#  Xtrain_fp=[AllChem.GetMorganFingerprint(mol, 2) for mol in Xtrain if mol is not None]
#Xtrain_ids=[mol.GetProp('_Name') for mol in Xtrain if mol is not None]
#Xtrain_ids
#labels=np.asarray(Xtrain_ids)
#labels
#Xtrain_fp
#Xtrain_fp
#Xtest_fp=[AllChem.GetMorganFingerprint(mol, 2) for mol in Xtest if mol is not None]
#Xtest_ids=[mol.GetProp('_Name') for mol in Xtest if mol is not None]
#Xtest_ids
#Xtest_fp
#from rdkit import DataStructs
#from rdkit.DataStructs import TanimotoSimilarity
# In[ ] #pairwise similarity
tc=[]
for i, fp1 in enumerate (Xtrain_fp):
    for j, fp2 in enumerate (Xtest_fp):
        tc.append(round(TanimotoSimilarity(fp1,fp2),2))

#tc=pd.DataFrame(tc)
#tc
#
k,l=[],[]
for i, id1 in enumerate (Xtrain_ids):
    for j, id2 in enumerate (Xtest_ids):
        k.append(id1)
        l.append(id2)
        
#k=pd.DataFrame(k)
#k
#l=pd.DataFrame(l)
#l
#tc_combined=pd.concat([k,l,tc],axis=1)
#tc_combined
#tc_combined.columns="Xtrain","Xtest","TC"
#tc_combined
#tc_combined.sort_values("TC" , ascending=False)
#tc_combined.to_csv ("Tc_combined.csv")
Building models and Cross Validation

#import pandas as pd
import numpy as np
#x=pd.read_csv("Dpp4 results of analysis.csv")
#x.head()
#x.shape
#
active_RS=x.loc[x["Standard Value"]<=10000]
inactive_RS=x.loc[x["Standard Value"]>10000]
#
active_RS["label"]=1
inactive_RS["label"]=0

#active_RS.head()
#active_RS.shape
#inactive_RS.shape
#inactive_RS.head()
#combined=pd.concat([active_RS,inactive_RS], axis=0)
#combined.head()
#combined.shape
#
combined[["Smiles","label"]].to_csv("dpp4_binary_labelled_model_.smi", sep=" \t", header=None, index=None)
#dataset = combined[["Smiles", "label"]]
#dataset.columns
#dataset["label"].value_counts()
#dataset=dataset.reset_index(drop=True)
#dataset.head()
#
import rdkit
from rdkit import Chem, DataStructs, RD Config
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
#
FP, IDS = [], []
for i, mol in enumerate(supplier):
    FP.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024))
    IDS.append(mol.GetProp("_Name"))
    
#morgan_fp =np.asarray(FP, dtype=np.int32)
#all_ids = np.asarray(IDS, dtype=np.int32).reshape( -1,1)
#morgan_fp
#from sklearn.model_selection import train_test_split
#
X_train, X_test, y_train, y_test = train_test_spli t(combined["Smiles"],
combined["label"], random_state=42, test_size=0.2, shuffle=True,
stratify=combined["label"])
#  X_train.head()
#  train_RS = pd.concat([X_train,y_train], axis=1)
#train_RS.head()
#train_RS.shape

# test_RS = pd.concat([X_test,y_test], axis=1)
# test_RS .shape
active_train=train_RS.loc[train_RS["label"]==1]
inactive_train=train_RS.loc[train_RS["label"]==0]
active_test=test_RS.loc[test_RS["label"]==1]
inactive_test=test_RS.loc[ test_RS["label"]==0]
#
active_train.to_csv("dpp4_r_active_train.smi",sep=" \t",header=None, index=None)
inactive_train.to_csv("dpp4_r_inactive_train.smi",sep=" \t",header=None, index=None)
active_test.to_csv("dpp4_r_active_test.smi",sep=" \t",header= None, index=None)
inactive_test.to_csv("dpp4_r_inactive_test.smi",sep=" \t",header=None, index=None)
#
yc=pd.concat([train_RS, test_RS], axis=0)
#
yc.shape
#
train_RS=train_RS.reset_index(drop=True)
test_RS=test_RS.reset_index(drop=True)
#
train_RS.head()
#
smii = [Chem.MolFromSmiles(mol) for mol in train_RS["Smiles"].iloc[0:]if mol is not None]
#
train_RS_fp=[AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in smii if mol is not None]

# train_RS_fp
# morgan_train_RS_fp=np.asarray(train_RS_fp, dtype=np.int32)
# test_RS.head()

smiii = [Chem.MolFromSmiles(mol) for mol in test_RS["Smiles"].iloc[0:]if mol is not None]
test_RS_fp=[AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in smiii if mol is not None]

# morgan_test_RS_fp=np.asarray(test_RS_fp, dtype=np.int32)
#
from rdkit import DataStructs
from rdkit.DataStructs import TanimotoSimilarity

#
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import Lab elEncoder
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import ShuffleSplit


# rf= RandomForestClassifier(random_state=42)
# param_rf = {"n_estimators": [100,200,1000]}
# cv=ShuffleSplit(n_splits=10, test_size=0.2, random_state=42)
# rf_gs = GridSearchCV(rf,param_rf,cv=cv, scoring="accuracy", verbose=10)
# In[ ] :rf_gs.fit(morgan_train_RS_fp,train_RS["label"].ravel())
#  predicted=rf_gs.predict(morgan_test_RS_fp)
# get_ipython().system('conda install -c sepandhaghighi pycm -y')
# from pycm import *
#
print(ConfusionMatrix(actual_vector = test_RS["label"].ravel(), predict_vector=predicted))

#
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB

#
svcclf = SVC()
svcclf.fit(morgan_train_RS_fp,train_RS["label"].ravel())
predictedsv=svcclf.predict(morgan_test_RS_fp)
print(ConfusionMatrix(actual_vector = test_RS["label"].ravel(), predict_vector=predictedsv))

NBCLF=GaussianNB()
NBCLF.fit(morgan_train_RS_fp,train_RS["label"].ravel())
predictedNB=NBCLF.predict(morgan_test_RS_fp)
print(ConfusionMatrix(actual_vector = test_RS["label"].ravel(), predict_vector=predictedNB))
KNNCLF=KNeighborsClassifier()
KNNCLF.fit(morgan_train_RS_fp,train_RS["label"].ravel())

predictedkn=KNNCLF.predict(morgan_test_RS_fp)
print(ConfusionMatrix(actual_vector = test_RS["label"].ravel(), predict_vector=predictedkn))
Random Split Model using Random forest Algorithm


import pandas as pd
import numpy as np

#x=pd.read_csv("Dpp4 results of analysis.csv")
#x.head()
#x.shape

active_RS=x.loc[x["Standard Value"]<=10000]
inactive_RS=x.loc[x["Standard Value"]>10000]

active_RS["label"]=1
inactive_RS["label"]=0

# active_RS.head()
# active_RS.shape
# inactive_RS.shape
# inactive_RS.head()
# combined=pd.concat([active_RS,inactive_RS], axis=0)
# combined.head()
# combined.shape
#
combined[["Smiles","label"]].to_csv("dpp4_binary_labelled_model_.smi",
sep=" \t", header=None, index=None)

# dataset = combined[["Smiles", "label"]]
# dataset.columns
# dataset["label"].value_counts()
# In [ ]:dataset=dataset.reset_index(drop=True)
# dataset.head()

import rdkit
from rdkit import Chem, DataStructs, RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys

supplier=Chem.SmilesMolSupplier("dpp4_binary_labelled_model_.smi",
delimiter=" \t", titleLine=False)

FP, IDS = [], []
for i, mol in enumerate(supplier):
    FP.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024))
    IDS.append(mol. GetProp("_Name"))

# morgan_fp =np.asarray(FP, dtype=np.int32)
# all_ids = np.asarray(IDS, dtype=np.int32).reshape( -1,1)
# morgan_fp
# from sklearn.model_selection import train_test_split
#
X_train, X_test, y_train, y_test = train_test_split(combined["Smiles"],
combined["label"], random_state=42, test_size=0.2, shuffle=True,
stratify=combined["label"])

# X_train.head()
# train_RS = pd.concat([X_train,y_train], axis=1)
# train_RS.head()
#train_RS.shape
# test_RS = pd.concat([X_test,y_test], axis=1)
# test_RS .shape
#
active_train=train_RS.loc[train_RS["label"]==1]
inactive_train=train_RS.loc[train_RS["label"]==0]
active_test=test_RS.loc[test_RS["label"]==1]
inactive_test=test_RS.loc[test_RS["label"]==0]
#
active_train.to_csv("dpp4_r_active_train.smi",sep=" \t",header=None,index=None)
inactive_train.to_csv("dpp4_r_inactive_train.smi",sep=" \t",header=None,index=None)
active_test.to_csv("dpp4_r_active_ test.smi",sep=" \t",header=None,index=None)
inactive_test.to_csv("dpp4_r_inactive_test.smi",sep=" \t",header=None,index=None)

# yc=pd.concat([train_RS, test_RS], axis=0)
# yc.shape
#
train_RS=train_RS.reset_index(drop=True)
test_RS=test_RS.reset_index(drop=True)

# train_RS.head()
smii = [Chem.MolFromSmiles(mol) for mol in train_RS["Smiles"].iloc[0:]if mol is not None]

train_RS_fp=[AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in smii if mol is not None]

# train_RS_fp
# morgan_train_RS_fp=np.asarray(train_RS_fp, dtype=np.int32)
# test_RS.head()
#
smiii = [Chem.MolFromSmiles(mol) for mol in test_RS["Smiles"].iloc[0:]if mol is not None]
test_RS_fp=[AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in smiii if mol is not None]

# morgan_test_RS_fp=np.asarray(test_RS_fp, dtype=np.int32)
#
from rdkit import DataStructs
from rdkit.DataStructs import TanimotoSim ilarity

#pairwise similarity
tc=[]
for i, fp1 in enumerate (train_RS_fp):
    for j, fp2 in enumerate (test_RS_fp):
        tc.append(round(TanimotoSimilarity(fp1,fp2),2))

# tc=pd.DataFrame(tc)
# tc

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn. model_selection import cross_val_score, cross_val_predict
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import ShuffleSplit


# rf= RandomForestClassifier(random_state=42)
# param_rf = {"n_estimators": [100,200,1000]}
# cv=ShuffleSplit(n_splits=10, test_size=0.2, random_state=42)

rf_gs = GridSearchCV(rf,param_rf,cv=cv, scoring="accuracy", verbose=10)

rf_gs.fit(morgan_train_RS_fp,train_RS["label"].ravel())
# predicted=rf_gs.predict(morgan_test_RS_fp)
# get_ipython().system('conda install -c sepandhaghighi pycm -y')
# from pycm import *
print(ConfusionMatrix(actual_vect or = test_RS["label"].ravel(), predict_vector=predicted))


# rf_gs.best_params_
# rf_best = rf_gs.best_estimator_
#
import pickle
with open("rf_model_dpp4_binary.pkl", "wb") as f:
    pickle.dump(rf_best, f)

with open("rf_model_dpp4_binary.pkl", "rb") as f:
    rf_clf=pickle.load(f)

# x=pd.read_csv("fda.csv")
# x.head()
#
x[["zinc_id","smiles"]].to_csv("zinc_dataset.smi", sep=" \t", header=None, index=None)

zinc = x[["smiles", "zinc_id"]]

# zinc.columns
# zinc.head()

zinc.to_csv("zinc_dataset.smi", sep=" \t", header=None, index=None)

zincsupplier=Chem.SmilesMolSupplier("zinc_dataset.smi", delimiter=" \t",titleLine=False)

from rdkit  import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import rdBase

Draw.MolsToGridImage([Chem.MolFromSmiles(smi) for smi in zinc.iloc[0:10,0] if smi is not None])

fp1=[rdkit.Chem.AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits = 1024)
for mol in zincsupplier if mol is not None]
train_nf=np.asarray(fp1, dtype=np.int)

# train_nf
# ids=[mol.GetProp("_Name")for mol in zincsupplier if mol is not None]
# train_nf.shape
# labels=np.asarray(ids, dtype=np.object).reshape( -1,1)
# labels
# In[ ] zinc_predicted=rf_clf.predict(train_nf)
#
zinc_df=pd.DataFrame(zinc_predicted)
zinc_labels=pd.DataFrame(labels)

predicted_results=pd.concat([zinc_df, zinc_labels], axis=1)

# predicted_results
# predicted_results=pd.concat([zinc_df, zinc_labels], axis=1)
# predicted_results
# predicted_results.columns=["predicted","zinc_id"]

positives=predicted_results.loc[predicted_results["predicted"]==1]

# positives
# zinc_predicted_prob=rf_clf.predict_proba(train_nf)
# zinc_predicted_prob

zinc_df_prob=pd.DataFrame(zinc_predicted_prob)
zinc_labels_prob=pd.DataFrame(labels)

predicted_results_prob=pd.concat([zinc_df_prob, zinc_labels_prob], axis=1)
predicted_results_prob.columns=["prob0", "prob1","zinc_id"]

# predicted_results_prob

positives_pb=pred icted_results_prob.loc[predicted_results_prob["prob1"]>=0.8]

# positives_pb.shape


# positives_pb.to_csv("positives_53.csv", sep=",", index=None)

Murcko Scaffold Cluster Based model .

import pandas as pd
import numpy as np

x=pd.read_csv("Dpp4 results of analysis.csv")

active_RS=x.loc[x["Standard Value"]<=10000]
inactive_RS=x.loc[x["Standard Value"]>10000]

active_RS["label"]=1
inactive_RS["label"]=0

combined=pd.concat([active_RS,inactive_RS], axis=0)

combined[["Smiles","label"]].to_csv("dpp4_binary_labelled_model_.smi",
sep=" \t", header=None, index=None)

dataset = combined[["Smiles", "label"]]

import molvs
from m olvs import Standardizer
from molvs import standardize_smiles

x1=combined[["Smiles","Molecule ChEMBL ID","label"]]
x1.shape
import rdkit
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import rdBase
rdBase.rdkitVersion
Draw.MolsToGridImage([Chem.MolFromSmiles(smi) for smi in  x1.iloc[0:, 0]])
[standardize_smiles(smi) for smi in x1.iloc[0:,  0]]
allsmiles=[standardize_smiles(smi) for smi in x1.iloc[0:, 0]]
allsmiles=pd.DataFrame(allsmiles)
x2=pd.concat([allsmiles,x1[["Molecule ChEMBL ID","label"]]], axis=1)
x2.columns="Smiles","ChEMBL ID","label"
x2.to_csv("data2new.smi" , sep=" \t" , index=None, header=None )
mol2=[Chem.MolFromSmiles(smi) for smi in allsmiles.iloc[0:, 0]]

from rdkit.Chem import SaltRemover
remover = SaltRemover.SaltRemover(defnData ="[Cl,Br]")
len(remover.salts)
mol2=[remover.StripMol(mol) for mol in mol2]
import rdkit.Chem as chem
Draw.MolsToGri dImage([chem.MolFromSmiles(smi) for smi in
x2["Smiles"].iloc[0:4821]])
dataset.to_csv("datanew.smi", sep=" \t", index=None , header=None)

from rdkit.Chem import Scaffolds
from rdkit.Chem import AllChem
import rdkit.Chem as Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
t1=Chem.SmilesMolSupplier("datanew.smi", delimiter=" \t", titleLine=False)
t1

[MurckoScaffold.GetScaffoldForMol(mol)for mol in t1 if mol is not None]
mol2=[MurckoScaffold.GetScaffoldForMol(mol)for mol in t1 if mol is not None]
mol2=[Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol)) for mol in
mol2 if mol2 is not None]
mol2
x2
ms=pd.DataFrame(mol2)
ms.head()
combined=pd.concat([x2, ms], axis=1)
combined.head()
combined.columns="Smiles","ChEMBL ID","label","Murcko"
combined[["Murcko","ChEMBL ID"]].to_csv("data_murcko.smi" , sep=" \t",
header=None ,index=None)
get_ipython().system('python/home/hemant/Downloads/NOOR_PYTHON/mayachemtools/maya/bin/RDKitClusterMolecules.py --infileParamssmilesColumn,1,smilesNameColumn,2,smilesDelimiter,tab,smilesTitleLine,auto,sanitize,yes -i data_murcko.smi -o murckoclusterfinal3.smi')
import pandas as pd
clus=pd.read_csv("murckoclusterfinal3.smi" , sep=" ")
clus
clus["ClusterNumber"].value_counts()
Draw.MolsToGridImage([Chem.MolFromSmiles(smi) for smi in clus.iloc[0:, 0] if
smi is not None])
clus
X_train=clus.groupby("ClusterNumber").filter(lambda x: len(x) >=2)
X_train.head()
X_test=clus.groupby("ClusterNumber").filter(lambda x: len(x) <2)
X_test.head()
combined.head()
combined=combined.rename(columns={"ChEMBL ID":"Name"})
train_CB = pd.merge(X_train, combined, on="Name")
train_CB
train_CB=train_CB[["Smiles","Name","label"]]
train_CB
test_CB = pd.merge(X_test, combined, on="Name")
test_CB
test_CB=test_CB[["Smiles","Name","label"]]
test_CB
train_CB=train_CB.reset_index(drop=True)
test_CB=test_CB.reset_index(drop=True)
train_CB["label"].value_counts()
smiiCB = [Chem.MolFromSmiles(mol) for mol in train_CB["Smiles"].iloc[0:]if
mol is not None]
train_CB_fp=[AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for
mol in smiiCB if mol is not None]
train_CB_fp
morgan_train_CB_fp=np.asarray(train_CB_fp, dtype=np.int32)
morgan_train_CB_fp
morgan_train_CB_fp.shape
test_CB.head()
smiiiCB = [Chem.MolFromSmiles(mol) for mol  in test_CB["Smiles"].iloc[0:]if mol is not None]
test_CB_fp=[AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in smiiiCB if mol is not None]
morgan_test_CB_fp=np.asarray(test_CB_fp, dtype=np.int32)
morgan_test_CB_fp
morgan_test_CB_fp.shape
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
rf= RandomForestClassifier(random_state=42)
from sklearn.model_selection import ShuffleSplit
param_rf = {"n_estimators": [100,200,1000]}
cv=ShuffleSplit(n_splits=10, test_size=0.2, random_state=42)
rf_gs = Gr idSearchCV(rf,param_rf,cv=cv, scoring="accuracy", verbose=10)
train_CB["label"].value_counts()
test_CB["label"].value_counts()
conda install -c conda -forge imbalanced -learn
rf_gs.fit(morgan_train_CB_fp,train_CB["label"].ravel())
predictedCB=rf_gs.predict(morgan_test_CB_fp)
from pycm import *
print(ConfusionMatrix(actual_vector = test_CB["label"].ravel(), predict_vector =predictedCB))
rf_gs.best_params_
rf_best = rf_gs.best_estimator_

import pickle

with open("rf_model_dpp4_binary.pkl", "wb") as f:
    pickle.dump(rf_best, f)

with open("rf_model_dpp4_binary.pkl", "rb") as f:
    rf_clf=pickle.load(f)

x=pd.read_csv("fda.csv")

x.head()
x[["zinc_id","smiles"]].to_csv("zinc_dataset.smi", sep=" \t", header=None,
index=None)
zinc = x[["smiles", "zinc_id"]]
zinc.columns
zinc.head()
zinc.to_csv("zinc_dataset.smi", sep=" \t", header=None, index=None)
zincsupplier=Chem.SmilesMolSupplier("zinc_dataset.smi", delimiter=" \t",
titleLine=False)

Draw.MolsToGridImage([Chem.MolFromSmiles(smi) for smi in zinc.iloc[0:10,0] if smi is not None])

fp1=[rdkit.Chem.AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits = 1024)
for mol in zincsupplier if mol is not None]
train_nf=np.asarray(fp1, d type=np.int)

train_nf

ids=[mol.GetProp("_Name")for mol in zincsupplier if mol is not None]
train_nf.shape
labels=np.asarray(ids, dtype=np.object).reshape( -1,1)
labels
zinc_predicted=rf_clf.predict(train_nf)
zinc_df=pd.DataFrame(zinc_predicted)
zinc_labels=pd.DataFrame(labels)
predicted_results=pd.concat([zinc_df, zinc_labels], axis=1)
predicted_results.columns=["predicted", "zinc_id"]
predicted_results
positives=predicted_results.loc[predicted_results["predicted"]==1]
positives


#


zinc_predicted_prob=rf_clf.predict_proba(train_nf)


#


zinc_predicted_prob


#


zinc_df_prob=pd.DataFrame(zinc_predicted_prob)
zinc_labels_prob=pd.DataFrame(labels)


#


predicted_results_prob=pd.concat([zinc_df_prob, zinc_labels_prob], axis=1)


#


predicted_results_prob.columns=["prob0", "prob1","zinc_id"]


#


predicted_results_prob


#


positives_pb=predicted_results_prob.loc[predicted_results_prob["prob1"]>=0.8]


#


positives_pb.shape


#


positives_pb.head()


#


positives_pb.to_csv("positives_52.csv", sep=",", index=None)


#


from sklearn.feature_selection import SelectFromModel


#


rf_cv_scoreCB =
cross_val_score(rf_gs,morgan_train_CB_fp,train_CB["label"].ravel(), cv=cv,
n_jobs= -1, scoring="roc_auc")


#


rf_cv_scoreCB


#


rf_cv_scoreCB.mean()


SMOTE BASED MODEL


#


import pandas as pd
import numpy as np


#


x=pd.read_csv("Dpp4 results of analysis.csv")


#


active_RS=x.loc[x["Standard Value"]<=10000]
inactive_RS=x.loc[x["Standard Value"]>10000]


#


active_RS["label"]=1
inactive_RS["label"]=0


#


combined=pd.concat([active_RS,inactive_RS], axis=0)


#


combined[["Smiles","label"]].to_csv("dpp4_binary_labelled_model_.smi",
sep=" \t", header=None, index=None)


#


data set = combined[["Smiles", "label"]]


#


import molvs
from molvs import Standardizer
from molvs import standardize_smiles


#


x1=combined[["Smiles","Molecule ChEMBL ID","label"]]


#


x1.shape


#


import rdkit
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import rdBase
rdBase.rdkitVersion


#


Draw.MolsToGridImage([Chem.MolFromSmiles(smi) for smi in  x1.iloc[0:, 0]])


#


[standardize_smiles(smi) for smi in x1.iloc[0:, 0]]


#


allsmiles=[standardize_smiles(smi) for smi in x1.iloc[0:, 0]]
allsmiles=pd.DataFrame(allsmiles)
x2=pd.concat([allsmiles,x1[["Molecule ChEMBL ID","label"]]], axis=1)


#


x2.columns="Smiles","ChEMBL ID","label"


#


x2.to_csv("data2new.smi" , sep=" \t" , index=None, header=None )


#


mol2=[Chem.MolFromSmiles(smi) for smi in allsmiles.iloc[0:, 0]]


#


from rdkit.Chem import SaltRemover
remover = SaltRemover.SaltRemover(defnData="[Cl,Br]")
len(remover.salts)


#


mol2=[remover.StripMol(mol) for mol in mol2]


#


import rdkit.Chem as chem


#


Draw.MolsToGridImage([chem.MolFromSmiles(sm i) for smi in
x2["Smiles"].iloc[0:4821]])


#


dataset.to_csv("datanew.smi", sep=" \t", index=None , header=None)


#


from rdkit.Chem import Scaffolds
from rdkit.Chem import AllChem
import rdkit.Chem as Chem


#


from rdkit.Chem.Scaffolds import MurckoScaffold


#


t1=Chem.SmilesMolSupplier("datanew.smi", delimiter=" \t", titleLine=False)


#


t1


#


[MurckoScaffold.GetScaffoldForMol(mol)for mol in t1 if mol is not None]


#


mol2=[MurckoScaffold.GetScaffoldForMol(mol)for mol in t1 if mol is not None]


#


mol2=[Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol)) for mol in
mol2 if mol2 is not None]


#


mol2


#


x2


#


ms=pd.DataFrame(mol2)


#


ms.head()


#


combined=pd.concat([x2, ms], axis=1)


#


combined.head()


#


combined.columns="Smiles","ChEMBL ID","label","Murcko"


#


combined[["Murcko","ChEMBL ID"]].to_csv("data_murcko.smi" , sep=" \t",
header=None ,index=None)


#


get_ipython().system('python
/home/hemant/Downloads/NOOR_PYTHON/mayachemtools/maya/bin/RDKit
ClusterMolecules.py --infileParams
smilesColumn,1,smilesNameColumn,2,smilesDelimite r,tab,smilesTitleLine,auto,
sanitize,yes -i data_murcko.smi -o murckoclusterfinal3.smi')


#


import pandas as pd


#


clus=pd.read_csv("murckoclusterfinal3.smi" , sep=" ")


#


clus


#


clus["ClusterNumber"].value_counts()


#


Draw.MolsToGridImage([Chem.MolFromSmiles(smi) for smi in clus.iloc[0:, 0] if
smi is not None])


#


clus


#


X_train=clus.groupby("ClusterNumber").filter(lambda x: len(x) >=2)


#


X_train.head()


#


X_test=clus.groupby("ClusterNumber").filter(lambda x: len(x) <2)


#


X_test.head()


#


combined.head()


#


combined=combined.rename(columns={"ChEMBL ID":"Name"})


#


train_CB = pd.merge(X_train, combined, on="Name")


#


train_CB


#


train_CB=train_CB[["Smiles","Name","label"]]


#


train_CB


#


test_CB = pd.merge(X_test, combined, on="Name")


#


test_CB


#


test_CB=test_CB[["Smiles","Name","label"]]


#


test_CB


#


train_CB=train_CB.reset_index(drop=True)
test_CB=test_CB.reset_index(drop=True)


#


train_CB["label"].value_counts()


#


smiiCB = [Chem.MolFromSmiles(mol) for mol in train_CB["Smiles"].iloc[0:]if
mol is not None]


#


train_CB_fp=[AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for
mol in smiiCB if mol is not None]


#


train_CB_fp


#


morgan_train_CB_fp=np.asarray(train_CB_fp, dtype=np.int32)


#


morgan_train_CB_fp


#


morgan_train_CB_fp.shape


#


test_CB.head()


#


smiiiCB = [Chem.MolFromSmiles(mol) for mol in test_CB["Smiles"].iloc[0:]if
mol is not None]
test_CB_fp=[AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for
mol in smiiiCB if mol is  not None]


#


morgan_test_CB_fp=np.asarray(test_CB_fp, dtype=np.int32)


#


morgan_test_CB_fp


#


morgan_test_CB_fp.shape


#


from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix


#


rf= RandomForestClassifier(random_state=42)


#


from sklearn.model_selection import ShuffleSplit


#


param_rf = {"n_estimators": [100,2 00,1000]}


#


cv=ShuffleSplit(n_splits=10, test_size=0.2, random_state=42)


#


rf_gs = GridSearchCV(rf,param_rf,cv=cv, scoring="accuracy", verbose=10)


#


train_CB["label"].value_counts()


#


test_CB["label"].value_counts()


#


conda install -c conda -forge imbalanced -learn


#


import imblearn
from imblearn.over_sampling import SMOTE


#


smote = SMOTE()


#


Xs_train,
ys_train=smote.fit_resample(morgan_train_CB_fp,train_CB["label"].ravel())


#


Xs_train.shape


#


rf_gs.fit(Xs_train, ys_train)


#


predictedCBim=rf_gs.predict(morgan_test_CB_fp)


#


from pycm import *


#


print(ConfusionMatrix(actual_vector = test_CB["label"].ravel(), predict_vector
=predictedCBim))


#


rf_gs.best_params_


#


rf_best = rf_gs.best_estimator_


#


import pickle

with open("rf_model_dpp4_binary.pkl", "wb") as f:
    pickle.dump(rf_best, f)


#


with open("rf_model_dpp4_binary.pkl", "rb") as f:
    rf_clf=pickle.load(f)


#


x=pd.read_csv("fda.csv")


#


x.head()


#


x[["zinc_id","smiles"]].to_csv("zinc_dataset.smi", sep=" \t", header=None,
index=None)


#


zinc = x[["smiles", "zinc_id"]]


#


zinc.columns


#


zinc.head()


#


zinc.to_csv("zinc_dataset.smi", sep=" \t", header=None, index=None)


#


zincsupplier=Chem.SmilesMolSupplier("zinc_dataset.smi", delimiter=" \t",
titleLine=False)


#


Draw.MolsToGridImage([Chem.MolFromSmiles(smi) for smi in zinc.iloc[0:10,
0] if smi is not None])


#


fp1=[rdkit.Chem.AllChem.GetMorganFingerprint AsBitVect(mol, 2, nBits = 1024)
for mol in zincsupplier if mol is not None]
train_nf=np.asarray(fp1, dtype=np.int)


#


train_nf


#


ids=[mol.GetProp("_Name")for mol in zincsupplier if mol is not None]


#


train_nf.shape


#


labels=np.asarray(ids, dtype=np.object).reshape( -1,1)


#


labels


#


zinc_predicted=rf_clf.predict(train_nf)


#


zinc_df=pd.DataFrame(zinc_predicted)
zinc_labels=pd.DataFrame(labels)


#


predicted_results=pd.concat([zinc_df, zinc_labels], axis=1)


#


predicted_results


#


predicted_results.columns=["predicted","zinc_id"]


#


positives=predicted_results.loc[predicted_results["predicted"]==1]


#


positives


#


zinc_predicted_prob=rf_clf.predict_proba(train_nf)


#


zinc_predicted_prob


#


zinc_df_prob=pd.DataFrame(zinc_predicted_prob)
zinc_labels_prob=pd.DataFrame(labels)


#


predicted_results_prob=pd.concat([zinc_df_prob,  zinc_labels_prob], axis=1)


#


predicted_results_prob.columns=["prob0", "prob1","zinc_id"]


#


predicted_results_prob


#


positives_pb=predicted_results_prob.loc[predicted_results_prob["prob1"]>=0.
8]


#


positives_pb.shape


#


positives_pb.head()


#


positives_pb.to_csv("positives_8.csv", sep=",", index=None)


#


from sklearn.feature_selection import SelectFromModel


#


rf_cv_scoreCB =
cross_val_score(rf_gs,morgan_train_ CB_fp,train_CB["label"].ravel(), cv=cv,
n_jobs= -1, scoring="roc_auc")


#


rf_cv_scoreCB


#


rf_cv_scoreCB.mean()