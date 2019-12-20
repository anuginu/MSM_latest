import numpy as np

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

def initial_population(imp_features):
  print("total number of features is ")
  print(len(imp_features))
  print("enter number of features")
  k=input()
  k=int(k)
  population_dihedral = []
  m=0
  while m < len(imp_features)-k:
   population_dihedral.append(imp_features[m:m+k])
   m=m+1
  end=len(imp_features)
  population_dihedral.append(imp_features[m:end])
  return population_dihedral


def initial_population_without_laplacian(imp_features):
  population_dihedral = [] 
  m=0
  features=list(range(0,len(imp_features)))
  while m < len(imp_features)-20:
   population_dihedral.append(features[m:m+20])
   m=m+1
  end=len(imp_features)
  population_dihedral.append(features[m:end]) 
  return population_dihedral



def calculate_fitness(population_dihedral,diheds,score_global,i,lock):
   import pandas as pd
   import numpy as np
   pop_index=i
   new_diheds = [] 
                                                        
   for i in range(0,len(diheds)):                                        
       X= diheds[i]                                                       
       selected_features = X[:, population_dihedral]                                    
       new_diheds.append(selected_features) 
   from msmbuilder.preprocessing import RobustScaler                                 
   scaler = RobustScaler()                                                           
   scaled_diheds = scaler.fit_transform(new_diheds)
   scaled_diheds=new_diheds                                  
   from msmbuilder.decomposition import tICA                                         
   tica_model = tICA(lag_time=2, n_components=5)
   tica_model.fit(scaled_diheds)                                         
   tica_trajs = tica_model.transform(scaled_diheds)                      
   from msmbuilder.cluster import MiniBatchKMeans                        
   clusterer = MiniBatchKMeans(n_clusters=200, random_state=42)

   clustered_trajs = clusterer.fit_transform(tica_trajs)   
   from msmbuilder.msm import MarkovStateModel       
   msm = MarkovStateModel(lag_time=50,n_timescales=5)                
   #msm.fit_transform(clustered_trajs)
   from sklearn.cross_validation import KFold            
   n_states = [4]           
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
          time_score = msm.timescales_[0]
          time_test_score= time_score+test_score
          print(time_score)
          print(test_score)
          av_score=time_test_score/2
          results.append({            
               'train_score': train_score,
               'test_score': test_score,
               'time_score': time_score,
               'av_score': av_score,
               'n_states': n,         
               'fold': fold}) 
          print(msm.timescales_)     
   results = pd.DataFrame(results)
   avgs = (results 
         .groupby('n_states')
         .aggregate(np.median)
         .drop('fold', axis=1))
   best_nt = avgs['test_score'].idxmax()
   best_n = avgs['av_score'].idxmax()
   best_score = avgs.loc[best_n, 'av_score']
   best_scorent=avgs.loc[best_nt,'test_score']
   print(best_scorent)
   lock.acquire()
   score_global.update({pop_index: best_scorent})
   lock.release()

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


def select_parents_roulette(score,population_dihedral):
    import random
    sum_fitness=0
    sum_temp=0
    parents = []
    score_keys=list(score.keys())
    for i in score_keys:
      sum_fitness=sum_fitness+score.get(i)
    for i in range(0,50):
      num=random.randint(1,sum_fitness+1)
      sum_temp=0
      for i in score_keys:
         sum_temp=sum_temp+score.get(i)
         if(sum_temp>=num):
           parents.append(population_dihedral[i])
           break

    return parents

def select_parents_rank_based(score,population_dihedral,cross_probability):
    import random
    import operator
    parents=[]
    num_parents=cross_probability*len(population_dihedral)
    num_parents=int(num_parents)
    if(num_parents%2==1):
     num_parents=num_parents+1
    scored_population=dict(sorted(score.items(), key=operator.itemgetter(1)))
    scored_pop=list(scored_population.keys())
    pop_new=[]
    for i in scored_pop:
      pop_new.append(population_dihedral[i])
    #rank=list(range(len(scored_pop),0,-1))
    rank=list(range(0,len(scored_pop)))
    rank_prob=[]
    value=0
    sum_rank_prob=0
    for i in rank:
     value=i/len(scored_pop)
     rank_prob.append(value)
    #sum of the rank_probability
    sum_rank_prob=0
    for i in range(len(rank_prob)):
      sum_rank_prob=sum_rank_prob+rank_prob[i]
    sum_rank_prob=int(sum_rank_prob)
    for j in range(0,num_parents):
      num=random.randrange(0,sum_rank_prob+1)
      sum_temp=0
      for i in range(len(rank_prob)):
         sum_temp=sum_temp+rank_prob[i]
         if(sum_temp>=num):
           parents.append(pop_new[i])
           break
    
    return parents


def parents_binarize(parents,imp_features):
    import numpy as np
    parents_binary=[]            
    for j in range(len(parents)):                         
     parents_binary.append(np.zeros(len(imp_features)))
     for i in  list(parents[j]):    
      parents_binary[j][i]=int(1)  
    return parents_binary


def crossover(parents,population_dihedral):       
    import numpy
    offspring = [] 
    #to find the chromosome with the minimum length
    min_length=[]
    for i in range(len(parents)):                    
      min_length.append(len(population_dihedral[0]))
    minimum=min(min_length)
   
    
    if(len(parents)%2==0):
     end = len(parents)-1
    else:
     end = len(parents)-1
    for k in range(0,end,2):
        #rand_ind= int(random.randrange(0,minimum-2) )
        crossover_point =(int)(minimum/2)
        parent1_idx = k                     
        parent2_idx = k+1             
        if(len(list(set(list(parents[parent1_idx][0:crossover_point])+list(parents[parent2_idx][crossover_point:]))))==len(parents[0])):       
          offspring.append(list(set(list(parents[parent1_idx][0:crossover_point])+list(parents[parent2_idx][crossover_point:]))))
        if(len(list(set(list(parents[parent2_idx][0:crossover_point])+list(parents[parent1_idx][crossover_point:]))))==len(parents[0])):
          offspring.append(list(set(list(parents[parent2_idx][0:crossover_point])+list(parents[parent1_idx][crossover_point:]))))
    return offspring






def crossover_binary(parents):       
    import numpy
    offspring = [] 
    #to find the chromosome with the minimum length
    min_length=[]
    for i in range(len(parents)):                    
      min_length.append(len(population_dihedral[0]))
    minimum=min(min_length)
   
    
    if(len(parents)%2==0):
     end = len(parents)-1
    else:
     end = len(parents)-2
    for k in range(0,end,2):
        rand_ind= int(random.randrange(0,minimum-2) )
        crossover_point = rand_ind
        parent1_idx = k                     
        parent2_idx = k+1 
        crossover_point = numpy.uint8(minimum/2)     
        if(len(list(parents[parent1_idx][0:crossover_point])+list(parents[parent2_idx][crossover_point:]))==len(parents[0])):       
         offspring.append(list(parents[parent1_idx][0:crossover_point])+list(parents[parent2_idx][crossover_point:]))
        if(len(list(parents[parent2_idx][0:crossover_point])+list(parents[parent1_idx][crossover_point:]))==len(parents[0])):
         offspring.append(list(parents[parent2_idx][0:crossover_point])+list(parents[parent1_idx][crossover_point:]))
    return offspring

def mutation_binary_offspring_withweightage(offsprings,rate,count_mutation,col_mean,imp_features):
    import random
    import numpy
    count=count_mutation

    for i in range(0,len(offsprings)):
     
     
     if count >0 : 
      index_0=[i for i,x in enumerate(offsprings[i]) if x==0] 
      index_1=[i for i,x in enumerate(offsprings[i]) if x==1]
      rank_index_1=[]
      for m in index_1:
         rank_index_1.append(col_mean[m])
      min_index_to_mutate=index_1[numpy.argmin(rank_index_1)]
      
      rank_index_0=[]
      for m in index_0:
         rank_index_0.append(col_mean[m])
      max_index_to_mutate=index_0[numpy.argmax(rank_index_0)]
      
      if min_index_to_mutate < len(offsprings[i]) :
        offsprings[i][min_index_to_mutate]=0
      if max_index_to_mutate < len(offsprings[i]):  
        offsprings[i][max_index_to_mutate]=1
      count=count-1
     else:
      break 
    return offsprings[0:count_mutation]  











def mutation_binary_offspring(offsprings,rate,count_mutation):
    import random
    count=count_mutation
    for i in range(0,len(offsprings)):
     if count >0 : 
      for j in range(rate):
       rand_ind= random.randrange(0,len(offsprings[i])-1)
       if offsprings[i][rand_ind]==1:
         offsprings[i][rand_ind]=0
         index_0=[i for i,x in enumerate(offsprings[i]) if x==0] 
         index_0 = [x for x in index_0 if x != rand_ind]
         #rand_ind_0=rank_next(offsprings[i]
         rand_ind_0=random.randrange(0,len(index_0)-1)
         offsprings[i][index_0[rand_ind_0]]=1
       else:
         offsprings[i][rand_ind]=1
         index_1=[i for i,x in enumerate(offsprings[i]) if x==1] 
         index_1 = [x for x in index_1 if x != rand_ind]
         #rand_ind_1=rank_next(offsprings[i]
         rand_ind_1=random.randrange(0,len(index_1)-1)
         offsprings[i][index_1[rand_ind_1]]=0
      count=count-1
     else:
      break 
      
    return offsprings[0:count_mutation]





def binary_to_pop_dih(offsprings):
    offspring_from_binary=[]
    for i in range(len(offsprings)):
       each_offs=[]
       for j in range(len(offsprings[i])):
          if(offsprings[i][j]==1):
             each_offs.append(j)
       offspring_from_binary.append(each_offs)
    return offspring_from_binary 

