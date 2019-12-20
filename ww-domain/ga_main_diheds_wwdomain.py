from  ga_functions import *
from ww_domain_feature import *
import numpy as np
def mutation_binary_offspring_withweightage(offsprings,rate,count_mutation,col_mean,imp_features):
    import random
    import numpy as np
    count=count_mutation

    for i in range(0,len(offsprings)):
     
     
     if count >0 : 
      index_0=[i for i,x in enumerate(offsprings[i]) if x==0] 
      index_1=[i for i,x in enumerate(offsprings[i]) if x==1]
      rank_index_1=[]
      for m in index_1:
         rank_index_1.append(col_mean[m])
      min_index_to_mutate=index_1[np.argmin(rank_index_1)]
      
      rank_index_0=[]
      for m in index_0:
         rank_index_0.append(col_mean[m])
      max_index_to_mutate=index_0[np.argmax(rank_index_0)]
      
      if min_index_to_mutate < len(offsprings[i]) :
        offsprings[i][min_index_to_mutate]=0
      if max_index_to_mutate < len(offsprings[i]):  
        offsprings[i][max_index_to_mutate]=1
      count=count-1
     else:
      break 
    return offsprings[0:count_mutation]  







def main_modified(generations):
 import numpy as np
 from msmbuilder.preprocessing import RobustScaler
 import time
 import pickle
 import os
 import multiprocessing
 os.environ["OMP_NUM_THREADS"] = "1"
 import operator
 from multiprocessing import Pool
 from operator import itemgetter
 diheds=Get_dihedral_features_wwdomain()
 scaler = RobustScaler()                                                           
 scaled_feature = scaler.fit_transform(diheds) 
 Val=Laplacian_score(scaled_feature) # output of imp_features and col_mean of the laplacian score of each  dihedral
 col_mean=Val[0]
 imp_features=Val[1]
 current_gen = 0
 for_each_gen_score =[]
 population_each_gen=[]
 population_dihedral=[]
 population_dihedral=initial_population(imp_features)
 cross_probability=0.8
 num_parents=(int)(cross_probability*len(population_dihedral))
 population_dihedral_duplicate=[]
 numberOfThreads = multiprocessing.cpu_count()
 f = open("benzamidine_diheds_ga_score"+str(generations)+".txt", "a")
 while current_gen < generations:
   manager = multiprocessing.Manager()
   score = manager.dict()
   processes = []
   lock = multiprocessing.Lock()
   for i in range(len(population_dihedral)):
         p = multiprocessing.Process(target=calculate_fitness, args=(population_dihedral[i],scaled_feature,score,i,lock))
         processes.append(p)
 #starttime = time.time()
   for i in chunks(processes,numberOfThreads): #chunks is a function : has to be defined
      p_count=0   
      for process in i:
         process.start() 
         p_count=p_count+1
      print("the started process are"+str(p_count))
      for process in i:

         process.join()
         p_count=p_count-1
      print("the joined process are"+str(p_count))
      for process in i:
         process.terminate()
         p_count=p_count+1
      print("the terminated process are"+str(p_count)) 
   scored_population={}
   scored_population=dict(sorted(score.items(), key=operator.itemgetter(1)))
   for_each_gen_score.append(scored_population)
   population_each_gen.append(population_dihedral)
   scored_population_list=list(scored_population.keys())
   parents=[]
   parents = select_parents_rank_based(scored_population,population_dihedral,cross_probability)
   offsprings_1=[]
   offsprings_1=crossover(parents,population_dihedral)
   
   parents_binary=[]
   parents_binary=parents_binarize(parents,imp_features)
   offsprings_2_binary=[]
   count_mutation=len(population_dihedral)-len(offsprings_1)
   offsprings_2_binary=mutation_binary_offspring(parents_binary,4,count_mutation)
#,col_mean,imp_features)
   offsprings_2=[]
   offsprings_2=binary_to_pop_dih(offsprings_2_binary)
   for i in range(len(offsprings_2)):
     offsprings_2[i]=np.asarray(offsprings_2[i])
   for i in range(len(offsprings_1)):
     offsprings_1[i]=np.asarray(offsprings_1[i])
   offsprings=[]
   offsprings=offsprings_1+offsprings_2
  # offsprings.append(population_dihedral[scored_population_list[len(scored_population_list)-1]])
   #offsprings.append(population_dihedral[scored_population_list[len(scored_population_list)-2]])
   population_dihedral=[]
   population_dihedral=offsprings
   current_gen = current_gen+1
 print(for_each_gen_score,file=f)
 f.close()
 return for_each_gen_score,population_each_gen,scaled_feature,imp_features  
