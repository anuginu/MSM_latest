def Get_dihedral_features_thrombin():
 import os 
 import shutil
 os.chdir('/homes/anuginueni/thrombin_trajs')
 if(os.path.isdir('/homes/anuginueni/thrombin_trajs/diheds')):  
   shutil.rmtree('/homes/anuginueni/thrombin_trajs/diheds')
 from msmbuilder.dataset import dataset
 xyz = dataset( "/homes/anuginueni/thrombin_trajs/*.xtc",topology='./top.pdb',stride=1) 
 from msmbuilder.featurizer import DihedralFeaturizer        #for dihedrals          
 featurizer = DihedralFeaturizer(types=['phi', 'psi'])       #for dihedrals
 diheds = xyz.fit_transform_with(featurizer, '/homes/anuginueni/thrombin_trajs/diheds/', fmt='dir-npy') #for dihedrals
 import mdtraj as md
 t=md.load( "/homes/anuginueni/thrombin_trajs/trajectory-10.xtc",top='/homes/anuginueni/thrombin_trajs/top.pdb',stride=5)
 des_feat=featurizer.describe_features(t)
 res = [ sub['resids'] for sub in des_feat ]
 print(str(res))
 return diheds


def Get_contacts_features_thrombin():
 import os 
 import shutil
 os.chdir('/homes/anuginueni/thrombin_trajs')
 if(os.path.isdir('./contacts')):  
   shutil.rmtree('./contacts')
 from msmbuilder.dataset import dataset
 xyz = dataset( "/homes/anuginueni/thrombin_trajs/*.xtc",topology='./top.pdb',stride=1) 
 from msmbuilder.featurizer import ContactFeaturizer        #for contacts          
 r1=list(range(1,277))
 index_list=[]
 for i in range(len(r1)):
   w=[]
   w.append(r1[i])
   w.append(278)
   index_list.append(w)
 print(len(index_list))
 print(index_list)
 featurizer = ContactFeaturizer(contacts=index_list,ignore_nonprotein=False)       #for contacts
 contacts = xyz.fit_transform_with(featurizer, 'contacts/', fmt='dir-npy') #for contacts
 import mdtraj as md
 t=md.load( "/homes/anuginueni/thrombin_trajs/trajectory-10.xtc",top='/homes/anuginueni/thrombin_trajs/top.pdb',stride=5)
 des_feat=featurizer.describe_features(t)
 res = [ sub['resids'] for sub in des_feat ]
 print(str(res))
 return contacts

def Get_rawposition_features_thrombin():
 import os 
 import shutil
 os.chdir('/homes/anuginueni/thrombin_trajs')
 if(os.path.isdir('./rawpositions')):  
   shutil.rmtree('./rawpositions')
 from msmbuilder.dataset import dataset
 xyz = dataset( "/homes/anuginueni/thrombin_trajs/*.xtc",topology='./top.pdb',stride=1) 
 from msmbuilder.featurizer import RawPositionsFeaturizer        #for contacts          
 featurizer = RawPositionsFeaturizer(raw_position_atom_index_pair)       #for contacts
 rawpositions = xyz.fit_transform_with(featurizer, 'rawpositions/', fmt='dir-npy') #for rawpositions
 return rawpositions

def Get_combined_features_thrombin():                                         
  from msmbuilder.featurizer import DihedralFeaturizer
  from msmbuilder.featurizer import ContactFeaturizer                            
  diheds= DihedralFeaturizer()
  contacts=ContactFeaturizer()
  features=[("di_thrombin",diheds),("con_thrombin",contacts)]
  import os
  import shutil
  os.chdir('/homes/anuginueni/thrombin_trajs')
  if(os.path.isdir('/homes/anuginueni/thrombin_trajs/combined')):
   shutil.rmtree('/homes/anuginueni/thrombin_trajs/combined')
  from msmbuilder.dataset import dataset
  xyz = dataset( "/homes/anuginueni/thrombin_trajs/*.xtc",topology='/homes/anuginueni/thrombin_trajs/top.pdb',stride=1)
  from msmbuilder.feature_selection import FeatureSelector

  comb_features=FeatureSelector(features)
  co=xyz.fit_transform_with(comb_features, '/homes/anuginueni/thrombin_trajs/combined/', fmt='dir-npy')
  return co

