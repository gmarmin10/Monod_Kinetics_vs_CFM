# Monod_Kinetics_vs_CFM
Quantifying key interactions between nutrient data and growth of phytoplankton using the Monod mathematical formula and the model, CFM-Phyto


Please download all the files in order for the main codes to run. The main codes here are titled, 1)"Monod_Github" and 2)"CFM_Github" 3)"macro_allocation". 

  1) "Monod_Github" produces two figures. Figure 1 plots the data against the Monod curve that is modeled by the Monod mathematical model. Figure 2 depicts the guess of maximum growth rate and K for each iteration. 
  2) "CFM_Github" produces 6 figures. Figure 1 plots the data against the model produced by the Cell flux model of phytoplankton. Figure 2 plots the allocation of nitrogen to carbon to groups of macromolecules with increasing nitrate concentration. Figure 3 plots the macromolecular allocation of carbon with increasing nitrate concentration. Figures 4 and 5 show the guesses for each parameter for each iteration the model runs. 
  3) "macro_allocation" plots the molecular allocation using CFM. 


For this purpose, an example dataset "Lee2019A.csv" was used and is also located within this folder to be downloaded. However, any dataset can be used. To change the dataset in the "Monod_Github" file, use line 15 and simply change the file name. Similarly, in the "CFM_Github" file, use line 167 and change the file name. 


If it is required to tune the model to the data, Lines 178-220 in the CFM code, are used to adjust the initial guesses of aN and Yn. In the Monod code, lines 25-66 can be used to do the same. To change the light intensity in the cell flux model, use line 49 and change the parameter labeled "I". 
