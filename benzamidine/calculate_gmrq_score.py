def calculate_scor(feature_list,extracted_feat):
  import pandas as pd
  import numpy as np
  new_diheds = []    
  print("feature list is")                                                   
  print(feature_list)
  for i in range(0,len(extracted_feat)):                                        
       X= extracted_feat[i]                                                       
       selected_features = X[:, feature_list]                                    
       new_diheds.append(selected_features) 
  from msmbuilder.preprocessing import RobustScaler                                 
  scaler = RobustScaler()                                                           
  scaled_diheds = [] 
  scaled_diheds = scaler.fit_transform(new_diheds)   
  from msmbuilder.decomposition import tICA                                         
  tica_model = tICA(lag_time=1, n_components=5, kinetic_mapping=True)
  tica_model.fit(new_diheds)                                         
  tica_trajs = tica_model.transform(new_diheds)                      
  from msmbuilder.cluster import MiniBatchKMedoids                        
  clusterer = MiniBatchKMedoids(n_clusters=50, random_state=42)
  clustered_trajs = clusterer.fit_transform(tica_trajs)   
  from msmbuilder.msm import MarkovStateModel       
  msm = MarkovStateModel(lag_time=50,n_timescales=5)
  from sklearn.cross_validation import KFold            
  n_states = [2]          
  cv = KFold(len(clustered_trajs), n_folds=5)                  
  results = []                                                 
  for n in n_states:                                           
      msm.n_states_ = n                                           
      for fold, (train_index, test_index) in enumerate(cv):       
          train_data = [clustered_trajs[i] for i in train_index]
          test_data = [clustered_trajs[i] for i in test_index]
          msm.fit(train_data)                
          train_score = msm.score(train_data)
          test_score = msm.score(test_data)
          results.append({            
               'train_score': train_score,
               'test_score': test_score,
               'n_states': n,         
               'fold': fold}) 
  msm.fit(clustered_trajs)
  print(msm.timescales_)             
  results = pd.DataFrame(results)
  print(results)
  avgs = (results 
         .groupby('n_states')
         .aggregate(np.var)
         .drop('fold', axis=1))
  avgs_mean=(results.groupby('n_states').aggregate(np.median).drop('fold',axis=1))
  best_n = avgs_mean['test_score'].idxmax()
  best_score_mean = avgs_mean.loc[best_n, 'test_score']

  best_n = avgs['test_score'].idxmax()
  best_score = avgs.loc[best_n, 'test_score']
  return best_score_mean
