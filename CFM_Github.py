'''
From a800_05_12_07_SameQc.py
'''

from pylab import *
from af001_energy_calculation import *
from Solver_2D_O import *
from Solver_3D_O import *
from Qc_essential_computation import *
from FigSetting2 import *
import random
import time
import numpy
from matplotlib.pyplot import figure, show, xlabel, ylabel, stackplot, plot, title

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Parameter sets (preparation)
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

m=3.79146798299876E-19                  #(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
Pmax=0.00320513285659728

OT=0.00863364097132997

Nstore_max=2.91679384515998E-15         #(molN cell-1) Constant protein pool in nitrogen (193-25)
Cnrna_variable=6212.59249917364         #(s) Constant for Variable part of RNA (193-26)
Ypthylakoid_chl=0.0281633095303638      #((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
Pconst_other=5.44534485638617E-17               #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
Qp_max=25.26/(3.097e16)                                                                              #(molP cell-1) total phosphorus content in the cell (193-26)
Cessential=1.51786753491048E-15         #(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%
Nconst_protein = 4.45336898828389E-15   #(molN cell-1) Constant protein pool in nitrogen (193-25)

#==============================

E3=evalue()
E=E3.E
Qc=1.00*10**(-12)/12                    #(molC/cell) biomass C per cell (196-18)(average of N and P limited cases from Healey 1985)
YchlN_C=4/55

#Conversion parameters================
CNprotein=4.49                          #(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
YcyanoC_N=2                             #(molC molN) C/N molar ratio of cyanophycin
YpgC_P=40                               #(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)
Nunit=1/Qc*14*10**6/(12*10**3)          #((ug N / mgC)/(molN cell-1) unit conversion term (164-20)
Punit=1/Qc*30.97*10**6/(12*10**3)       #((ug P / mgC)/(molP cell-1) unit conversion term (164-20)

#Photosynthesis================

I=20                                    #(umolE m-2 s-1) light intensity
Pchl=Pmax*(1-exp(-OT*I))                #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
Pchl=(14/24)*Pchl
#=========================================
#DNA and RNA-> Kei 193-28
#=========================================              

Molar_mass_DNA_AT_average=307.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_DNA_CG_average=307.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_RNA_AT_average=316.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_RNA_CG_average=323.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")

#------------------------------------
#E coli
#------------------------------------
CG_Ecoli=0.506                          #(dimensionless) from [http://www.ncbi.nlm.nih.gov/genome/167 (accessed 06/18/2016)]
AT_Ecoli=1-CG_Ecoli                     #(dimensionless) 

Molar_mass_DNA_Ecoli=Molar_mass_DNA_AT_average*CG_Ecoli+Molar_mass_DNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of DNA unit
Molar_mass_RNA_Ecoli=Molar_mass_RNA_AT_average*CG_Ecoli+Molar_mass_RNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of RNA unit

RNA_DNA_mass_ratio=17.844/6.5239        #(ug/ug) from values ad D=0 "07 Bremer and Dennis 1996 data plot.xlsx"
RNA_DNA_molar_ratio=RNA_DNA_mass_ratio/Molar_mass_RNA_Ecoli*Molar_mass_DNA_Ecoli               #(mol mol-1)

#---------------------------------------------
#Stoichiometric parameters for DNA and RNA
#---------------------------------------------

CG=0.563                                #GC% not CG but I started with CG so I stick with it; it does not matter as "AT GC".   [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
YnucacidP_N=1/(3.5*(1-CG)+4*CG)         #(molP molN-1) P/N molar ratio of RNA (193-26) values (193-28) excel file "08 N to P ratio in DNA and RNA.xlsx"

YdnaC_N=3.5*(1-CG)+2.5*CG               #(molC molN-1) C/N molar ratio of dna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
YrnaC_N=3.25*(1-CG)+2.5*CG              #(molC molN-1) C/N molar ratio of rna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)

DNAmb=2.1269                            #(Mb) Megabase pair of synechococcus DNA in mega (million) base pairs [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
Avogadro=6.022*10**23                   #(molecules mol-1) Avogadro constant
Pdna_const=DNAmb*2*10**6/Avogadro       #(molP cell-1) Constant part of DNA in phosphorus 
Prna_const=Pdna_const*RNA_DNA_molar_ratio       #(molP cell-1) Constant part of RNA in phosphorus
#* Make sure to multiply by 2 as they are base PAIRs"
Ndna_const=Pdna_const/YnucacidP_N       #(molN cell-1) Constant part of DNA in nitrogen
Nrna_const=Ndna_const*RNA_DNA_molar_ratio       #(molN cell-1) Constatn part of RNA in phosphorus
Ndna=Ndna_const                         #(molN cell-1) DNA in nitrogen (here assuming constant)

Converter = 4.34728279914354E-10/3.56099164557551    #For fixed ratios of two parameters: Cnbiosynth and Ynphoto_chl. The value based on original CFM-Phyto

def function(aNO3,Ynphoto_chl,NO3):
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Main calculation 199-21~
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    Vn=aNO3*NO3                             #(mol N cell-1 s-1) nitrogen uptake per cell 
    Cnbiosynth= Ynphoto_chl*Converter       #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
    A=((1+E)*Qc*Ynphoto_chl)/Pchl+Cnbiosynth
    B=Nconst_protein+(m*Ynphoto_chl)/Pchl
    G=((1+E)*Qc*YchlN_C)/Pchl
    H=(m*YchlN_C)/Pchl
    I=A*Cnrna_variable
    J=A+B*Cnrna_variable+G
    K=B+Nrna_const+Ndna+H
    L = ((1 + E)*Qc*Ypthylakoid_chl)/Pchl
    M = (m*Ypthylakoid_chl)/Pchl
    
    aN=I
    bN=J
    cN=K
    dN=-Vn
    
    DN=solver_3D(aN,bN,cN,dN)
    D=DN.X
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Obtaining output values (Similar to the previous steady state one)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #==============================
    #Carbon related: Solver setting up (for C) (Kei 200-33)
    #==============================
    #Preparation---------------------------------------
    #Chlorophyll related
    Chl_const = m/Pchl                              # (molC chl cell-1) cN[i]hlrophyll concentration (193-25)
    Chl_D = (1 + E)*Qc/Pchl                           # (molC chl cell-1) cN[i]hlrophyll concentration (193-25)
    
    #Nitrogen related
    Nchl_D = Chl_D*YchlN_C                          # (molN chl cell-1) Chlorophyll N concentration
    Nphoto_D = Chl_D*Ynphoto_chl                    # (molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nchl_const = Chl_const*YchlN_C                  # (molN chl cell-1) Chlorophyll N concentration
    Nphoto_const = Chl_const*Ynphoto_chl            # (molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nbiosynth_D = Cnbiosynth                        # (molN cell-1) various part of biosynthesis related protein in N (193-37)
    Nprotein_D = Nphoto_D + Nbiosynth_D             # (molN cell-1) All the proteins in N (193-26)
    Nprotein_const = Nphoto_const + Nconst_protein  # (molN cell-1) All the proteins in N (193-26)
    Nrna_D = Nprotein_const*Cnrna_variable          # (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
    Nrna_D2 = Nprotein_D*Cnrna_variable             # (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
    
    #Constant carbon parameters------------------
    Cconst_protein = Nconst_protein*CNprotein       #(molC cell-1) carbon in other protein assumed constant (195-16)
    Cdna_const = Ndna_const*YdnaC_N                 #(molC cell-1) carbon in constant part of DNA (195-16)   
    
    #Calculating factors and Qc_essential-------------
    
    Qc_D2 = A*Cnrna_variable*YrnaC_N
    Qc_D = (1+E)*Qc/Pchl + A*CNprotein + L*YpgC_P + B*Cnrna_variable*YrnaC_N
    Qc_const = m/Pchl + B*CNprotein + M*YpgC_P + Nrna_const*YrnaC_N + Ndna_const*YdnaC_N + Cessential
    
    Qc_essential = Qc_essential_computation(D,Chl_const, Chl_D, Nchl_const, Nchl_D, Nphoto_const, Nphoto_D, Nbiosynth_D,\
                                                  Nprotein_const, Nprotein_D, Nrna_const, Nrna_D, Nrna_D2, Ypthylakoid_chl, YrnaC_N,\
                                                  CNprotein, YpgC_P, Cdna_const, Cconst_protein, Cessential)  
    
    i0 = Qc_essential >= Qc                         #conditions where Mu max applies
    i1 = isnan(Qc_essential)                        #Including this avoide when Qc_essential is nan with high Vn
    i = i0 + i1
    D[i] = 2*(Qc - Qc_const)/(Qc_D + sqrt(Qc_D*Qc_D + 4*Qc_D2*(Qc - Qc_const)))
    
    return D
# function ending here

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# AW method: Largely copied from J001 00 00
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
data = genfromtxt("Lee2019A.csv",delimiter=',').T
NO3 = data[0]                                       #(uM) NO3 concentration
Mu_d = data[1]                                      #(d-1) Growth rate

#Removing negative values
NO3[Mu_d<0] = nan
Mu_d[Mu_d<0] = nan

#Standard deviation of samples
sigma = nanstd(Mu_d)*sqrt(size(Mu_d)-sum(isnan(Mu_d)))/sqrt(size(Mu_d)-sum(isnan(Mu_d))-1)

# Initial values
aN = 0.05*Qc/86400                                  #aNO3 #(mol N cell-1 s-1 / umolN L-1) affinity of NO3 (initial value from d007_06_00)
Yn = 3.56099164557551                               #Ynphoto_chl     #((molN cell-1)/(molC chl cell-1)) the stoichiometric ratio for cell photosynthetic enzyme (Rubisco etc.) nitrogen to chlorophyll (193-25)

#Array preparation
repeat_times = arange(100000+1)
n = zeros(size(repeat_times))*nan
aN_array = copy(n)
Yn_array = copy(n)
Nc_array = copy(n)

#change ratio per iteration
a = 0.5
a_aN = a*aN
a_Yn = a*Yn

#For step 0
Mu0 = function(aN,Yn,NO3)
X2 = nansum((Mu_d - Mu0)**2/(2*sigma**2))
P = exp(-X2)
Pbest = P

t0 = time.time()
for i in repeat_times:

    #checking if the parameter value come within the realistic range
    while True:
        aN1 = aN + a_aN*random.uniform(-1,1)
        if aN1>0 and aN1<10e-19:
            break
    
    while True:
        Yn1 = Yn + a_Yn*random.uniform(-1,1)
        if Yn1>0 and Yn1<10:
            break
    
    #Candidate value calculation
    Mu1 = function(aN1,Yn1,NO3)*86400
    X2 = nansum((Mu_d - Mu1)**2/(2*sigma**2))
    P1 = exp(-X2)
    Pratio = P1/P
    r01 = random.uniform(0,1)

    #testing the quality of Pratio
    if r01 < Pratio:
        P = P1
        aN = aN1
        Yn = Yn1
        
        #Pbest update
        if P>Pbest:
            Pbest = P
            aNbest = aN
            Ynbest = Yn
    
    #Recording values
    aN_array[i] = aN
    Yn_array[i] = Yn
    
    #Time printing
    if mod(i,10000) == 0:
        print(i,round(time.time()-t0,2),'(s)')
        t0 = time.time()

#for loop end here
#===================

#Printing best values
print('Pbest=', Pbest,', aNbest=', aNbest, ', Ynbest=', Ynbest)
#Running the model for plotting
NO3_high_res = arange(0,100,1)
Mu = function(aNbest,Ynbest,NO3_high_res)
D = Mu

Ynphoto_chl = Ynbest
Cnbiosynth = Ynphoto_chl*Converter

#Preparation---------------------------------------
#Chlorophyll related
Chl_const = m/Pchl                              # (molC chl cell-1) cN[i]hlrophyll concentration (193-25)
Chl_D = (1 + E)*Qc/Pchl                         # (molC chl cell-1) cN[i]hlrophyll concentration (193-25)

#Nitrogen related
Nchl_D = Chl_D*YchlN_C                          # (molN chl cell-1) Chlorophyll N concentration
Nphoto_D = Chl_D*Ynphoto_chl                    # (molN cell-1) Photosynthesis related protein nitrogen (193-25)
Nchl_const = Chl_const*YchlN_C                  # (molN chl cell-1) Chlorophyll N concentration
Nphoto_const = Chl_const*Ynphoto_chl            # (molN cell-1) Photosynthesis related protein nitrogen (193-25)
Nbiosynth_D = Cnbiosynth                        # (molN cell-1) various part of biosynthesis related protein in N (193-37)
Nprotein_D = Nphoto_D + Nbiosynth_D             # (molN cell-1) All the proteins in N (193-26)
Nprotein_const = Nphoto_const + Nconst_protein  # (molN cell-1) All the proteins in N (193-26)
Nrna_D = Nprotein_const*Cnrna_variable          # (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
Nrna_D2 = Nprotein_D*Cnrna_variable             # (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)

#Constant carbon parameters------------------
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# for plotting
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
Chl       = Chl_const      + Chl_D*D
#=================
#Nitrogen related
#=================
Nchl      = Nchl_const     + Nchl_D*D          #(molN chl cell-1) Chlorophyll N concentration
Nphoto    = Nphoto_const   + Nphoto_D*D        #(molN cell-1) Photosynthesis related protein nitrogen (193-25)
Nbiosynth =                  Nbiosynth_D*D
Nprotein  = Nprotein_const + Nprotein_D*D      #(mol N cell-1) all the protein
Nrna      = Nrna_const     + Nrna_D*D + Nrna_D2*D*D
Nessential= Nchl + Nphoto + Nbiosynth + Nconst_protein + Nrna + Ndna
Qn_max    = Nessential + Nstore_max

#=========================
#Phosphorus related
#=========================
Pdna=Ndna*YnucacidP_N                          #(mol P cell-1) DNA in phosphorus
Pthylakoid=Chl*Ypthylakoid_chl                 #(molP cell-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
Prna=Nrna*YnucacidP_N                          #(molP cell-1) Phosphorus in RNA 
Pessential=Pthylakoid+Prna+Pdna+Pconst_other   #(molP cell-1)

#Constant carbon parameters------------------
Cconst_protein = Nconst_protein*CNprotein      #(molC cell-1) carbon in other protein assumed constant (195-16)
Cdna_const = Ndna_const*YdnaC_N                #(molC cell-1) carbon in constant part of DNA (195-16)   
#==================================
#Unite conversion
#==================================
Numbertoarray=ones(size(NO3_high_res))         #(dimensionless) Number to array converter
#----------------------------------
#For N/C
#----------------------------------
Nconst_protein_plot=Nconst_protein*Nunit*Numbertoarray      #(ug N/ mgC) constant protein pool in nitrogen (193-25)(193-33)
Nphoto_plot=Nphoto*Nunit                                    #(ug N/ mgC) Photosynthesis related protein nitrogen (193-25)(193-33)
Nbiosynth_plot=Nbiosynth*Nunit                              #(ug N/ mgC) biosynthesis related protein in N (193-37)
Ndna_const_plot=Ndna_const*Nunit*Numbertoarray              #(ug N/ mgC) Nitrogen in constant part of DNA
Nrna_plot=Nrna*Nunit                                        #(ug N/ mgC) Nitrogen in constant part of RNA
Nchl_plot=Nchl*Nunit                                        #(ug N/ mgC) Chlorophyll nitrogen (actually almost negligiable) (193-33)

#==============================
#Carbon related 
#==============================
Cchl = Chl                          #(molC cell-1) carbon in chlorophyll (195-16)
Cphoto = Nphoto*CNprotein           #(molC cell-1) carbon in photosystem protein (195-16)
Cbiosynth = Nbiosynth*CNprotein     #(molC cell-1) carbon in biosynthesis protein (195-16)
Crna = Nrna*YrnaC_N                 #(molC cell-1) carbon RNA (195-16)
CthylakoidPG = Pthylakoid*YpgC_P    #(molC cell-1) carbon in PG (phosphatidyl glycerol) in thylakoid membranes


Cother = Qc - Cphoto - Cbiosynth - Cconst_protein - Cchl\
        - Crna - Cdna_const\
        - Cessential - CthylakoidPG 

#=======================================
#Unit conversions for C
#=======================================

percentorratio=100                  #100: percent, 1:ratio
Cphoto_plot=Cphoto/Qc*percentorratio           
Cbiosynth_plot=Cbiosynth/Qc*percentorratio
Cconst_protein_plot=Cconst_protein/Qc*percentorratio*Numbertoarray
Cchl_plot=Cchl/Qc*percentorratio
Crna=Crna/Qc*percentorratio*Numbertoarray
Cdna_const_plot=Cdna_const/Qc*percentorratio*Numbertoarray
Cother_plot=Cother/Qc*percentorratio
Cessential_plot=Cessential/Qc*percentorratio*Numbertoarray
CthylakoidPG_plot=CthylakoidPG/Qc*percentorratio


#======================
# Color setting
#======================

Color_DNA_const='orange'
Color_RNA_const='yellow'
Color_protein_const='blue'
Color_photo='yellow'
Color_RNA_variable='red'
Color_DNA_variable='blue'
Color_protein_biosynthesis='pink'
Color_chl='#66FF66'
Color_other='#006400'
Color_other='#008000'
Color_P_const='#CCFFFF'
Color_DNA='black'
Color_RNA='red'
Color_Nstore='purple'
Color_Pstore='#D0CECE'
Color_Cessential='brown'
Color_thylakoid='#FFD966'

#======================

StackPlotColorsC=(Color_Cessential,Color_protein_const,Color_photo,Color_protein_biosynthesis,Color_DNA,Color_RNA,Color_chl,Color_thylakoid,\
                  Color_Nstore, Color_other)
StackPlotColorsN=(Color_protein_const,Color_photo,Color_protein_biosynthesis,Color_DNA,Color_RNA,Color_chl,Color_Nstore)
StackPlotColorsP=(Color_P_const,Color_thylakoid,Color_DNA,Color_RNA,Color_Pstore)

Xlabel = 'NO$_{3}^{-}$ ($\mu$mol L$^{-1}$)'

figure(1)
plot(NO3,Mu_d,'o')
plot(NO3_high_res,D*86400)
xlabel('NO$_{3}^{-}$ ($\mu$M)')
ylabel('$\mathit{\mu}$ (d$^{-1}$)')

figure(2)
stackplot(NO3_high_res,Nconst_protein_plot,Nphoto_plot,Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,colors=StackPlotColorsN)
title("N:C")
xlabel(Xlabel)
ylabel("N:C ($\\mu$g N / mg C)")

figure(4)
stackplot(NO3_high_res,Cessential_plot,Cconst_protein_plot,Cphoto_plot,Cbiosynth_plot,Cdna_const_plot,\
          Crna,Cchl_plot,CthylakoidPG_plot,Cother_plot,colors=StackPlotColorsC)
xlabel(Xlabel, fontsize=30)                       #copied from 73
ylabel('C allocation ($\%$)', fontsize=30)  #copied from 73

figure(5)
plot(repeat_times,aN_array)
ylabel('aN')
xlabel('repeat time')

figure(6)
plot(repeat_times,Yn_array)
ylabel('Yn')
xlabel('repeat time')

show()
    
    