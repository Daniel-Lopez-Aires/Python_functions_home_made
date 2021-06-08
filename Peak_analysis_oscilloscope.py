#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 14:37:16 2021

@author: dla

function to process the data from the histograms for B.1.4. experiment.

*Inputs:
        .Voltages, Times = variables that contains all the datasets, sepparating the
            datasets by columns (np.arrays)
        .Column_index = colum that contains the data we want to analyse
        .Index that identifies the values (voltage and times), which comes from a variable
            storing all the values
        . Indexes that limits the peak
        .delta_V, delta_t = error of the measurements, from the oscilloscope
        
*Outputs:
        .Peak value (V)
        .Baseline (V)
        .delta_V and delta_t
        .Integral of the peak and its error
        .Number of elements of the peak
"""

#######0) General packages useful###########

import matplotlib.pyplot as plt  #for simplicity, to not write matplotlib.pyplot
        #everytime we want to plot something
import scipy.optimize              #to do the fit. doing only import scipy sometimes
                                #gives an error, so have to do this
import numpy as np          #np contain linspaces as np.linspace(a,b,N)
from scipy import stats as stats     #to find the most commom element in a lis
####



def Peak_analysis_oscillo(Voltages, Times, column_index, index_min_peak, index_max_peak, delta_V, delta_t):
    
    #TO find the peak I could easily do (absolute value to be able to detect positive peaks (raw signal) and 
    		#negative signals (from pre) ):
    peak = max( np.absolute(Voltages[:,column_index]) )                           #[V] Peak value
    index_peak = np.where( np.absolute(Voltages[:,column_index]) == peak )        #index

    #using the index of the peak I see the interval by looking at the .csv file
    #and the plot

    voltage_peak_neg = Voltages[:,column_index][index_min_peak-1:index_max_peak-1]  
    time_peak = Times[:,column_index][index_min_peak-1:index_max_peak-1]
    voltage_peak = np.array(np.absolute(voltage_peak_neg))   #abs() because this
        #peak contains both > and ,0 values, so to add them in order to count
        #them, I have to put everything in the positive value
        
    len_peak =  len(Voltages[:,column_index][index_min_peak-1:index_max_peak-1])           #len of the peak   

    #Baseline.
        #to choose the baseline, the 1st approach was to do the mean between channel 0 and the
        #channel where the peak start to appear. But this is abd idea becaue the error is very high,
        #so the 2nd approach will simply be choose the most frquent value

    #len_baseline_signal_points = len(Voltages[column_index][0:index_min_peak-1])  #length of the voltages used
        #to compute the baseline. Will be neccesary for error calcs
    #baseline_signal = sum(Voltages[column_index][0:index_min_peak-1]) / len_baseline_signal_points    
        #[V] baseline, to compute the peak amplitude. I average between the initial value and
        #when the peak starts
        

    baseline = stats.mode(Voltages[:,column_index])[0]       #[V] baseline voltage, the most common value
                    #[0]contains the value, [1] the frequency

        
    #sum_voltage = sum(-voltage_peak + baseline)         #[V] total voltage of the peak, corrected with
                        #the baseline, so that this is the real voltage created!! Baseline-voltage because
                        #voltage of the amplitude is the greatest. DO NOT NEED IT ANYMORE!!!!

    integral = np.trapz(voltage_peak, time_peak)        #[V*t] area under the peak
    delta_integral = 2 * np.sqrt( delta_V / (max(voltage_peak) - min(voltage_peak) ) 
                             + delta_t/ (max(time_peak) - min(time_peak) ) )    #overstimation of the
                #error of the integral. This is the error of the 
                #area of the rectangle (Vmax-Vmin)*(tmax-tmin)


    #len_baseline_signal_points_stored.append(len_baseline_signal_points)

    #Plot (debug)
    #plt.plot([1e6 * x for x in Times[-1]],[1e3 * x for x in Voltages[-1]], 'bo-')   
                                                    #-1 chooses last element, which is the
    plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
    plt.plot( 1e6 * Times[:,column_index], 1e3 * Voltages[:,column_index], 'b.-')    
    plt.plot( 1e6 * time_peak , 1e3 * voltage_peak_neg, 'r.-')   
        #widht so that each bar touches each other!
    plt.title("Waveform", fontsize=22)           #title
    plt.xlabel("time (us)", fontsize=14)                        #xlabel
    plt.ylabel("voltage (mV)", fontsize=14)              #ylabel
    plt.legend(['data', 'peak interval'], fontsize=14) 
    # Set size of tick labels.
    plt.tick_params(axis='both', labelsize=14)              #size of axis
    plt.grid(True) 
    #plt.savefig('Signal_LYSO_sum.png', format='png') 
    
    
   #2) Return of values########################
   #the values will be returned in a dictionary indicating what is each
   #value
    values = {'peak[V]' : peak, 'baseline[V]' : baseline, 
              'N_peak' : len_peak, 
              '\Delta(V[V])' : delta_V, '\Delta(t[s])' : delta_t,
              'integral[V*s]' : integral, '\Delta(integral[V*s])' : delta_integral
              }
    return values
