'''
Created on Nov 7, 2016

@author: Keisuke
'''



def Qc_essential_computation(D,Chl_const, Chl_D, Nchl_const, Nchl_D, Nphoto_const, Nphoto_D, Nbiosynth_D,\
        Nprotein_const, Nprotein_D, Nrna_const, Nrna_D, Nrna_D2, Ypthylakoid_chl, YrnaC_N,\
        CNprotein, YpgC_P, Cdna_const, Cconst_protein, Cessential):
    
    #=======================
    #N computation for i+1
    #=======================
    Chl       = Chl_const      + Chl_D*D
    Nchl      = Nchl_const     + Nchl_D*D               #(molN chl cell-1) Chlorophyll N concentration
    Nphoto    = Nphoto_const   + Nphoto_D*D        #(molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nbiosynth =                  Nbiosynth_D*D
    Nprotein  = Nprotein_const + Nprotein_D*D
    Nrna      = Nrna_const     + Nrna_D*D + Nrna_D2*D*D
    
    #=========================
    #Phosphorus related
    #=========================
    Pthylakoid = Chl*Ypthylakoid_chl          #(molP cell-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
    
    #=============================
    #preparing for Qc computation
    #=============================
    Crna = Nrna*YrnaC_N       #(molC cell-1) carbon in RNA (195-16)
    Cchl = Chl                    #(molC cell-1) carbon in chlorophyll (195-16)
    Cphoto = Nphoto*CNprotein     #(molC cell-1) carbon in photosystem protein (195-16)
    Cbiosynth = Nbiosynth*CNprotein   #(molC cell-1) carbon in biosynthesis protein (195-16)
    CthylakoidPG = Pthylakoid*YpgC_P           #(molC cell-1) carbon in PG (phosphatidyl glycerol) in thylakoid membranes
    
    #=============================
    #Computing Qc_essential
    #=============================
    Qc_essential = Crna + Cdna_const + Cchl + Cphoto + Cbiosynth\
                        +Cconst_protein + CthylakoidPG + Cessential
                        
    return Qc_essential