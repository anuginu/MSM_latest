def Get_contacts_features_wwdomain():
 import os 
 import shutil
 os.chdir('/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein')
 if(os.path.isdir('/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein/contacts')):  
   shutil.rmtree('/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein/contacts')
 from msmbuilder.dataset import dataset
 xyz = dataset( "/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein/*.dcd",topology='/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein/gtt-1-protein.pdb',stride=1) 
 from msmbuilder.featurizer import ContactFeaturizer        #for contacts          
 featurizer = ContactFeaturizer()       #for contacts
 contacts = xyz.fit_transform_with(featurizer, '/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein/contacts', fmt='dir-npy') #for contacts
 return contacts
#BEN RESID IS 225
def Get_dihedral_features_wwdomain():
 import os 
 import shutil
 os.chdir('/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein')
 if(os.path.isdir('/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein/dihedrals')):  
   shutil.rmtree('/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein/dihedrals')
 from msmbuilder.dataset import dataset
 xyz = dataset( "/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein/*.dcd",topology='/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein/gtt-1-protein.pdb',stride=1) 
 from msmbuilder.featurizer import DihedralFeaturizer        #for contacts          
 featurizer =  DihedralFeaturizer(types=['phi', 'psi'])     #for contacts
 dihedrals = xyz.fit_transform_with(featurizer, '/home/anu/Downloads/DESRES-Trajectory_GTT-1-protein/GTT-1-protein/dihedrals', fmt='dir-npy') #for contacts
 return dihedrals
