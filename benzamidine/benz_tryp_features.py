def Get_dihedral_features_benzamidine():
 import os 
 import shutil
 import pandas as pd
 import numpy as np
 os.chdir('/homes/anuginueni/benzamidine_trajs')
 if(os.path.isdir('/homes/anuginueni/benzamidine_trajs/diheds')):  
   shutil.rmtree('/homes/anuginueni/benzamidine_trajs/diheds')
 from msmbuilder.dataset import dataset
 xyz = dataset( "/homes/anuginueni/benzamidine_trajs/*.xtc",topology='/homes/anuginueni/benzamidine_trajs/top.pdb',stride=1) 
 from msmbuilder.featurizer import DihedralFeaturizer        #for dihedrals          
 featurizer = DihedralFeaturizer(types=['phi', 'psi'])       #for dihedrals
 diheds = xyz.fit_transform_with(featurizer, '/homes/anuginueni/benzamidine_trajs/diheds/', fmt='dir-npy') #for dihedrals
 import mdtraj as md
 t=md.load( "/homes/anuginueni/benzamidine_trajs/trajectory-10.xtc",top='/homes/anuginueni/benzamidine_trajs/top.pdb',stride=5)
 des_feat=featurizer.describe_features(t)
 res = [ sub['resids'] for sub in des_feat ]
 print(str(res)) 
 return diheds

def Get_contacts_features_benzamidine():
 import os 
 import shutil
 os.chdir('/homes/anuginueni/benzamidine_trajs')
 if(os.path.isdir('/homes/anuginueni/benzamidine_trajs/contacts')):  
   shutil.rmtree('/homes/anuginueni/benzamidine_trajs/contacts')
 from msmbuilder.dataset import dataset
 xyz = dataset( "/homes/anuginueni/benzamidine_trajs/*.xtc",topology='/homes/anuginueni/benzamidine_trajs/top.pdb',stride=1) 
 from msmbuilder.featurizer import ContactFeaturizer        #for contacts          
 r1=list(range(1,223))
 index_list=[]
 for i in range(len(r1)):
   w=[]
   w.append(r1[i])
   w.append(224)
   index_list.append(w)
 print(len(index_list))
 print(index_list)
 featurizer = ContactFeaturizer(contacts=index_list,ignore_nonprotein=False)       #for contacts
 contacts = xyz.fit_transform_with(featurizer, '/homes/anuginueni/benzamidine_trajs/contacts/', fmt='dir-npy') #for contacts
 import mdtraj as md
 t=md.load( "/homes/anuginueni/benzamidine_trajs/trajectory-10.xtc",top='/homes/anuginueni/benzamidine_trajs/top.pdb',stride=5)
 des_feat=featurizer.describe_features(t)
 res = [ sub['resids'] for sub in des_feat ]
 print(str(res))
 return contacts


def Get_rawposition_features_benzamdine():
 import os 
 import shutil
 os.chdir('/homes/anuginueni/benzamidine_trajs')
 if(os.path.isdir('/homes/anuginueni/benzamidine_trajs/rawpositions')):  
   shutil.rmtree('/homes/anuginueni/benzamidine_trajs/rawpositions')
 from msmbuilder.dataset import dataset
 xyz = dataset( "/homes/anuginueni/benzamidine_trajs/*.xtc",topology='/homes/anuginueni/benzamidine_trajs/top.pdb',stride=1) 
 from msmbuilder.featurizer import RawPositionsFeaturizer        #for rawposition          
 featurizer = RawPositionsFeaturizer()       #for rawposition
 rawpositions = xyz.fit_transform_with(featurizer, '/homes/anuginueni/benzamidine_trajs/rawpositions/', fmt='dir-npy') #for rawpositions
 return rawpositions

def Get_combined_features_benzamidine():                                                     
  from msmbuilder.featurizer import DihedralFeaturizer
  from msmbuilder.featurizer import ContactFeaturizer                                           
  diheds= DihedralFeaturizer()
  contacts=ContactFeaturizer()
  features=[("di_fspeptide",diheds),("con_fspeptide",contacts)] 
  import os
  import shutil
  os.chdir('/homes/anuginueni/benzamidine_trajs')
  if(os.path.isdir('/homes/anuginueni/benzamidine_trajs/combined')):
   shutil.rmtree('/homes/anuginueni/benzamidine_trajs/combined')
  from msmbuilder.dataset import dataset
  xyz = dataset( "/homes/anuginueni/benzamidine_trajs/*.xtc",topology='/homes/anuginueni/benzamidine_trajs/top.pdb',stride=1)
  from msmbuilder.feature_selection import FeatureSelector 
   
  comb_features=FeatureSelector(features)
  co=xyz.fit_transform_with(comb_features, '/homes/anuginueni/benzamidine_trajs/combined/', fmt='dir-npy')
  return co



def Laplacian_score(diheds):
  import scipy.io
  import numpy
  import os
  #os.chdir('/home/anu/Downloads/scikit-feature-1.0.0')
  from skfeature.function.similarity_based import lap_score
  from skfeature.utility import construct_W
  from numpy import mean
  kwargs_W = {"metric": "euclidean", "neighbor_mode": "knn", "weight_mode": "heat_kernel", "k": 5, 't': 1}
  idx = []
  #change the path for every system to be run.
  #os.chdir('/home/anu/Downloads/traj_benz_trypsin/')
  for i in range(len(diheds)):
   X= diheds[i]
   W = construct_W.construct_W(X, **kwargs_W)
   score = lap_score.lap_score(X, W=W)
   idx.append(score)
  col_mean = mean(idx, axis =0)
  imp_features = numpy.argsort(col_mean)
  return col_mean,imp_features

def spec_score(diheds):
  import scipy.io
  import numpy
  from numpy import mean
  import os
  #os.chdir('/home/anu/Downloads/scikit-feature-1.0.0')
  from skfeature.function.similarity_based import SPEC


  idx = []
  #change the path for every system to be run.
  #os.chdir('/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein')
  for i in range(0,len(diheds),5):
   X= diheds[i]
   kwargs = {'style':0}
   score = SPEC.spec(X, **kwargs)
   print(score)
   idx.append(score)
  col_mean = mean(idx, axis =0)
  idx=SPEC.feature_ranking(col_mean,**kwargs)
  return col_mean,idx

def mcfs_score(diheds):
  import scipy.io
  import numpy
  from numpy import mean
  import os
  #os.chdir('/home/anu/Downloads/scikit-feature-1.0.0')
  from skfeature.function.sparse_learning_based import MCFS
  from skfeature.utility import construct_W
  from skfeature.utility import unsupervised_evaluation
  idx = []
  kwargs = {"metric": "euclidean", "neighborMode": "knn", "weightMode": "heatKernel", "k": 5, 't': 1}
  #change the path for every system to be run.
  #os.chdir('/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein')
  for i in range(0,len(diheds),5):
   X= diheds[i]
   W = construct_W.construct_W(X, **kwargs)
   score = MCFS.mcfs(X, n_selected_features=20, W=W, n_clusters=20)

   idx.append(score)
  col_mean = mean(idx, axis =0)
  imp_features=MCFS.feature_ranking(col_mean)
  return col_mean,imp_features

def ndfs_score(diheds):
  import scipy.io
  import numpy
  from numpy import mean
  import os
  #os.chdir('/home/anu/Downloads/scikit-feature-1.0.0')
  from skfeature.function.sparse_learning_based import NDFS
  from skfeature.utility import construct_W
  from skfeature.utility import unsupervised_evaluation
  from skfeature.utility.sparse_learning import feature_ranking
  idx = []
  kwargs = {"metric": "euclidean", "neighborMode": "knn", "weightMode": "heatKernel", "k": 5, 't': 1}
  #change the path for every system to be run.
  #os.chdir('/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein')
  for i in range(0,len(diheds),5):
   X= diheds[i]
   W = construct_W.construct_W(X, **kwargs)
   score = NDFS.ndfs(X, W=W, n_clusters=20)

   idx.append(score)
  col_mean = mean(idx, axis =0)
  imp_features=feature_ranking(col_mean)
  return col_mean,imp_features

def udfs_score(diheds):
  import scipy.io
  import numpy
  from numpy import mean
  import os
  #os.chdir('/home/anu/Downloads/scikit-feature-1.0.0')
  from skfeature.function.sparse_learning_based import UDFS
  from skfeature.utility.sparse_learning import feature_ranking
  idx = []
    # sort the feature scores in an ascending order according to the feature scores
  for i in range(0,len(diheds),5):
   X= diheds[i]
   score = UDFS.udfs(X, gamma=0.1, n_clusters=20)
   idx.append(score)
  col_mean = mean(idx, axis =0)
  imp_features= feature_ranking(col_mean)
  return col_mean,imp_features

