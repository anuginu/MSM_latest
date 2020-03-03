# MSM_latest
To find the optimum residues that can be used as collective variable to construct the most efficient MSM.
Different folders- Benzamidine, Fspeptide, Thrombin-pip,villin and ww-domain contains the scripts to obtain the collective variable for the corresponding bio-molecules. 
To run the script in ipython, following are the steps to be performed 
1) %run ga_functions.py
2) %run fs_peptidefeatures.py / %run benz_trypfeatures.py/%run thrombin_features.py etc

To obtain contact features are collective variable 

3a) %run ga_main_contacts_fspeptide.py/%run ga_main_contacts_benzamidine.py/%run ga_main_contacts_thrombin.py etc

To obtain dihedral features are collective variable 

3b) %run ga_main_diheds_fspeptide.py/%run ga_main_diheds_benzamidine.py/%run ga_main_diheds_thrombin.py etc

Then run 
4) main_modified(#generations) where #generations is the number of iterations to be run for the GA.  
