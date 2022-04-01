'''
Created on Feb 14, 2022

@author: keiin
'''

from pylab import *
from FigSetting2 import *
from sf import *
import random
import time
from matplotlib.pyplot import figure, show, xlabel,ylabel, subplot, plot
import statistics

data = genfromtxt("Lee2019A.csv",delimiter=',').T
NO3 = data[0]                                       #(uM) NO3 concentration
Mu_d = data[1]                                      #(d-1) Growth rate

#Removing negative values
NO3[Mu_d<0] = nan
Mu_d[Mu_d<0] = nan

sigma = nanstd(Mu_d)*sqrt(size(Mu_d)-sum(isnan(Mu_d)))/sqrt(size(Mu_d)-sum(isnan(Mu_d))-1)

#Initial values
K = 20                                              #(uM) Initial Kno3: half saturation constant for NO3
Mumax = 0.4                                         #(d-1) Maximum growth rate 

#Array preparation
repeat_times = arange(100000+1)
n = zeros(size(repeat_times))*nan
K_array = copy(n)
Mumax_array = copy(n) 

#Defining a function
def function(K,Mumax,NO3):
    Mu = Mumax*(NO3/(K+NO3))
    return Mu

#change ratio per iteration
a = 0.2
a_K = a*K
a_Mumax = a*Mumax

#For step 0
Mu0 = function(K,Mumax,NO3)
X2 = nansum((Mu_d - Mu0)**2/(2*sigma**2))
P = exp(-X2)
Pbest = P
Kbest = K
Mumaxbest = Mumax

t0 = time.time()
for i in repeat_times:
    
    #checking if the parameter value come within the realistic range
    while True:
        K1 = K + a_K*random.uniform(-1,1)
        if K1>0 and K1<200:
            break
    
    while True:
        Mumax1 = Mumax + a_Mumax*random.uniform(-1,1)
        if Mumax1>0 and Mumax1<60/24:
            break
    
    #Candidate value calulation
    Mu1 = function(K1,Mumax1,NO3)
    X2 = nansum((Mu_d - Mu1)**2/(2*sigma**2))
    P1 = exp(-X2)
    Pratio = P1/P
    r01 = random.uniform(0,1)   
    
    #testing the quality of Pratio
    if r01 < Pratio:
        P = P1
        K = K1
        Mumax = Mumax1
        
        #Pbest update
        if P>Pbest:
            Pbest = P
            Kbest = K
            Mumaxbest = Mumax
    
    #Recording values
    K_array[i] = K
    Mumax_array[i] = Mumax

    #Time printing
    if mod(i,10000) == 0:
        print(i,round(time.time()-t0,2),'(s)')
        t0 = time.time()

#for loop end here
#===================

#Printing best values
print('Pbest=',Pbest,', Kbest =',Kbest,', Mumaxbest=',Mumax)
#Running the model for plotting
NO3_high_res = arange(0,130,1)
Mu = function(Kbest,Mumaxbest,NO3_high_res)

#Plotting
figure(1)
plot(NO3,Mu_d,'o')
plot(NO3_high_res,Mu,)
xlabel('NO$_{3}^{-}$ ($\mu$M)')
ylabel('$\mathit{\mu}$ (d$^{-1}$)')

figure(2)
subplot(2,1,1)
plot(repeat_times,K_array)
xlabel('Repeat times')
ylabel('K ($\mu$mol L$^{-1}$)')

subplot(2,1,2)
plot(repeat_times,Mumax_array)
xlabel('Repaet times')
ylabel('Mumax (d$^(-1)$')

show()




