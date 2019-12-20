from benz_tryp_features import *
from ga_functions import *
from calculate_gmrq_score import *
dih=Get_dihedral_features_benzamidine()
#cont=Get_contacts_features_benzamidine()
#comb=Get_combined_features_benzamidine()

score_spec_dih=spec_score(dih)
score_lap_dih=Laplacian_score(dih)
score_mcfs_dih=mcfs_score(dih)
score_ndfs_dih=ndfs_score(dih)
score_udfs_dih=udfs_score(dih)
score_spec_feature=[]
for i in range(50,450,50):
  dih_gmrq=calculate_scor(score_spec_dih[1][0:i],dih)
  score_spec_feature.append(dih_gmrq)

print(score_spec_feature) 
f=open("benzamidine_dih_spec_ranking.txt","a+")  
f.write(str(score_spec_feature))
f.close()
score_lap_feature=[]
for i in range(50,450,50):
  dih_lap=calculate_scor(score_lap_dih[1][0:i],dih)
  score_lap_feature.append(dih_lap)

dih_lap=calculate_scor(score_lap_dih[1][0:888],dih)
score_lap_feature.append(dih_lap)
f=open("benzamdine_dih_lap_ranking.txt","a+")
f.write(str(score_lap_feature))
f.close()



score_mcfs_feature=[]
for i in range(50,450,50):
  dih_mcfs=calculate_scor(score_mcfs_dih[1][0:i],dih)

  score_mcfs_feature.append(dih_mcfs)
dih_mcfs=calculate_scor(score_mcfs_dih[1][0:888],dih)
score_mcfs_feature.append(dih_mcfs)
f=open("benzamidine_dih_mcfs_ranking.txt","a+")
f.write(str(score_mcfs_feature))
f.close()


score_ndfs_feature=[]
for i in range(50,450,50):
  dih_ndfs=calculate_scor(score_ndfs_dih[1][0:i],dih)

  score_ndfs_feature.append(dih_ndfs)
f=open("benzamidine_dih_ndfs_ranking.txt","a+")
dih_ndfs=calculate_scor(score_ndfs_dih[1][0:888],dih)
score_ndfs_feature.append(dih_ndfs)
f.write(str(score_ndfs_feature))
f.close()



score_udfs_feature=[]
for i in range(50,450,50):
  dih_udfs=calculate_scor(score_udfs_dih[1][0:i],dih)

  score_udfs_feature.append(dih_udfs)
f=open("benzamidine_dih_udfs_ranking.txt","a+")
dih_udfs=calculate_scor(score_udfs_dih[1][0:888],dih)
score_udfs_feature.append(dih_udfs)
f.write(str(score_udfs_feature))
f.close()

