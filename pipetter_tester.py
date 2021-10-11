# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 12:19:09 2021

@author: lopedan
"""

#######0)General packages useful#################

import matplotlib.pyplot as plt  #for simplicity, to not write matplotlib.pyplot
        #everytime we want to plot something

#from scipy.stats import norm               ##norm.fit() fit to gaussian
import numpy as np
    #np contain linspaces as np.linspace(a,b,N)
import pandas as pd

##########



def Pipetting_tester(pipette_volume, mass_measurements, delta_m = .00001):
    '''
    Python function to test your pipetting skills
    pippete volume in mL
        *Inputs:
        .pipette_volume: volume of the pippete [mL]
        .mass_measurements: arrays (not numpy) with the measurements from the balance [g].
        .delta_m = error of the balance. 
        
    *Outputs:
        .Dictionary with lot of things: slope and its error, intercept and its error, correlation
        	coefficient
    '''
    mass_measurements = np.array(mass_measurements)     #conversion to numpy array
    N = len(mass_measurements)                                        #number of measurements
    
    ###### 1) Calcs ###################
    m_inc =  np.array( [mass_measurements[i+1] - mass_measurements[i] 
                                 for i in range(0,N-1)] ) #[g] Mass increment
    
    delta_m_inc = np.sqrt(2) * delta_m     #[g] error of the mass increment


    #The mean and its error of the mass increment is:
    
    m_mean = sum(m_inc) / len(m_inc)   #[g] Mean value of the mass increment
    delta_m_mean = np.sqrt( len(m_inc) ) * delta_m_inc / len(m_inc) #[g] error of
                                #the value

    #The standar deviation is:
            # sigma = sqrt(1/N * sum_i(x_i - <x>) ) 
    sigma = np.sqrt( 1/len(m_inc) * sum( (m_inc - m_mean)**2 ) )  #[g]


    #This guys computes the systematic and random error as (ranodm is associated to you,
    #while systematic to the instrument):

    syst_error = 100* np.abs(m_mean - pipette_volume) /pipette_volume   #[%] Systematic error. 
            #Note that pipette volume in mL = mass in g
    rand_error = sigma / m_mean * 100                           #[%] Random error

    
    ######### 2) Bar making ####################
    #To make a bar plot with the error, first I have to create an histogram, store
    #its variables, and use them for the bar plot
    
    counts, bin_edges = np.histogram(m_inc, bins = 15)
    
    #The bin center can be computed easily, provided that the previous line returns 
    #all the borders of the bars (so, if there is 3 bars, 4 borders):
    
    bin_center = np.array([ (bin_edges[i+1] + bin_edges[i]) /2 for 
                       i in range(0,len(bin_edges)-1 )] )       #center of the bins
    
    
    ############# 3) Plot ######################
    plt.figure(figsize=(10,8))  #width, heigh 6.4*4.8 inches by default
    plt.bar(bin_center, counts, width=bin_edges[1] - bin_edges[0], 
        yerr = delta_m_inc, edgecolor="black")
    plt.title("Pipeting skills. In theory $\Delta$(m) = " + str(pipette_volume) + "g", fontsize=20)           #title
    plt.xlabel("$\Delta$(m) [g]", fontsize=14)                        #xlabel
    plt.ylabel("Counts", fontsize=14)              #ylabel
    # Set size of tick labels.
    plt.tick_params(axis='both', labelsize=14)              #size of axis
    plt.grid(True) 
    plt.savefig('Pipeting_skills.png', format='png') 
    
    
    ########## 4) Results ###########
    #Randome error should be less than 0.2%, systematic depend on the quantity:
        
    print("Random error: " + str(rand_error) + "%. Ideally <= 0.2%")       #random error printing
    
    #Since the random error depends on the quantity, ahve to do if statements:
    if pipette_volume == 1: #1ml
        print("Systematic error: " + str(syst_error) + "%. Ideally <= 0.6%")       #random error printing
   
    elif pipette_volume == .5: #.5ml
        print("Systematic error: " + str(syst_error) + "%. Ideally <= 0.1%")       #random error printing
            
    elif pipette_volume == .1: #.1ml
        print("Systematic error: " + str(syst_error) + "%. Ideally <= 3%")       #random error printing
           
            
           
    #################5) Return of values############
        #the values will be returned in a dictionary indicating what is each
        #value
    values = {'Delta_m[g]' : m_inc, 'delta(Delta_m[g])' : delta_m_inc , 
              '<Delta_m>[g]' : m_mean, 'delta(<Delta_m>[g])' : delta_m_mean,
              'sigma[g]' : sigma, 'Systematic_error[%]' : syst_error, 'Random_error[%]' : rand_error}
    
    return values