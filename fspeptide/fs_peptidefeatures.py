#for Fs peptide
def Get_dihedral_features_fspeptide():                                                     
                                             
  from msmbuilder.featurizer import DihedralFeaturizer                   
  from msmbuilder.example_datasets import FsPeptide                     
  trajs = FsPeptide().get().trajectories 
  trajs = [t[::20] for t in trajs]                                
  featurizer = DihedralFeaturizer(types=['phi', 'psi'])  
  diheds = featurizer.fit_transform(trajs)                               
  import mdtraj as md
  t=md.load( "/homes/anuginueni/msmbuilder_data/fs_peptide/trajectory-10.xtc",top='/homes/anuginueni/msmbuilder_data/fs_peptide/fs-peptide.pdb',stride=5)
  des_feat=featurizer.describe_features(t)
  res = [ sub['resids'] for sub in des_feat ]
  print(str(res))
  return diheds               

def Get_contacts_features_fspeptide():
  from msmbuilder.featurizer import DihedralFeaturizer                   
  from msmbuilder.example_datasets import FsPeptide                     
  trajs = FsPeptide().get().trajectories 
  trajs = [t[::10] for t in trajs]   
  from msmbuilder.featurizer import ContactFeaturizer        #for contacts          
  featurizer = ContactFeaturizer()
  contacts = featurizer.fit_transform(trajs)
  import mdtraj as md
  t=md.load( "/homes/anuginueni/msmbuilder_data/fs_peptide/trajectory-10.xtc",top='/homes/anuginueni/msmbuilder_data/fs_peptide/fs-peptide.pdb',stride=5)
  des_feat=featurizer.describe_features(t)
  res = [ sub['resids'] for sub in des_feat ]
  print(str(res))
  return contacts

def Get_rawposition_features_fspeptide():
  from msmbuilder.featurizer import DihedralFeaturizer                   
  from msmbuilder.example_datasets import FsPeptide                     
  trajs = FsPeptide().get().trajectories 
  trajs = [t[::10] for t in trajs]   
  from msmbuilder.featurizer import RawPositionsFeaturizer
  featurizer = RawPositionsFeaturizer()
  rawpositions=featurizer.fit_transform(trajs)
  return rawpositions
