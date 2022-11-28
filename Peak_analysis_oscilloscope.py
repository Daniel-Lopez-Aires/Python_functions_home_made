#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 14:37:16 2021

@author: Daniel López Aires // danlopair@gmail.com

"""

################### 0) General useful packages ##################

import matplotlib.pyplot as plt  #for simplicity, to not write matplotlib.pyplot
        #everytime we want to plot something
import numpy as np          #np contain linspaces as np.linspace(a,b,N)
from scipy import stats    #to find the most commom element in a lis
##################################



def Peak_analysis_oscillo(Voltage, Time, signal_type, index_min_peak, index_max_peak, delta_V, delta_t):
    '''
    This function is to do the Pulse Shape Analysis (PSA) of a electric pulse (nuclear detectors)
from the oscilloscope.

    *Inputs:
        .Voltage, Time = variables (vectors) that contain the data (np.array)
        .Signal_type = 'pre' or 'raw'. This indicated whether the peak is positive ('raw')
        	or negative ('pre'), which is neccessary to compute the amplitude
        .Indexes: indexes that limits the peak, eye spotted (will be plotted to confirm the selection)
        .delta_V, delta_t = error of the measurements, from the oscilloscope

    *Outputs:
        .Peak value (V)
        .Baseline (V)
        .delta_V and delta_t
        .Integral of the peak and its error
        .Number of elements of the peak
        .Rise and decay time of the peak, and their error (same for both)
    '''

    ############### 1) Finding the peak (max/min value) ################
    
    	#to find the peak I could easily do (absolute value to be able to detect positive peaks (raw signal) and 
    		#negative signals (from pre) ):
    if signal_type == 'pre':                   #signal from pre ==> negative peak
                
        peak = min( Voltage )                           #[V] Peak value
        peak_abs = max( np.absolute(Voltage) )          #[V] |Peak value|
    
    elif signal_type == 'raw':                  #signal from raw ==> positive peak
        peak = max( Voltage )             #[V] |Peak value|
        peak_abs = max( np.absolute(Voltage) )          #[V] |Peak value|
        
    
    index_peak = np.where( np.absolute(Voltage) == peak_abs )[0][0]        #index
            #index_peak is a tuple, which contains a single np.aray, so with [0][0] you get the value
    #using the index of the peak I see the interval by looking at the .csv file
    #and the plot

    voltage_peak_neg = Voltage[index_min_peak-1:index_max_peak-1]  
    time_peak = Time[index_min_peak-1:index_max_peak-1]
    voltage_peak = np.absolute(voltage_peak_neg)   #abs() because this
        #peak contains both > and ,0 values, so to add them in order to count
        #them, I have to put everything in the positive value
        
    len_peak =  len(Voltage[index_min_peak-1:index_max_peak-1])           #len of the peak   


    ################### 2) Finding the baseline #####################
    
        #to choose the baseline, the 1st approach was to do the mean between channel 0 and the
        #channel where the peak start to appear. But this is a ad idea becaue the error is very high,
        #so the 2nd approach will simply be choose the most frquent value

    #len_baseline_signal_points = len(Voltages[column_index][0:index_min_peak-1])  #length of the voltages used
        #to compute the baseline. Will be neccesary for error calcs
    #baseline_signal = sum(Voltages[column_index][0:index_min_peak-1]) / len_baseline_signal_points    
        #[V] baseline, to compute the peak amplitude. I average between the initial value and
        #when the peak starts
        

    baseline = stats.mode(Voltage)[0][0]       #[V] baseline voltage, the most common value
                    #[0]contains the values (array), [1] the frequency (array), so another [0]
                    #is needed to get the value
        
    #sum_voltage = sum(-voltage_peak + baseline)         #[V] total voltage of the peak, corrected with
                        #the baseline, so that this is the real voltage created!! Baseline-voltage because
                        #voltage of the amplitude is the greatest. DO NOT NEED IT ANYMORE!!!!


    ####################### 3) Computing the integral of the peak ###################
    #The area under the peak (integral) can be computed numerically using the trapz method.
    #NEED TO THINK OF A WAY TO STIMATE ITS ERROR!!! I tried to put the error of the area of a rectange:
    # (Vmax-Vmin)*(tmax-tmin).


    integral = np.trapz(voltage_peak, time_peak)        #[V*t] area under the peak
    delta_integral = 2 * np.sqrt( delta_V / (max(voltage_peak) - min(voltage_peak) ) 
                             + delta_t/ (max(time_peak) - min(time_peak) ) )    #overstimation of the
                #error of the integral. This is the error of the 
                #area of the rectangle (Vmax-Vmin)*(tmax-tmin)


    #len_baseline_signal_points_stored.append(len_baseline_signal_points)


    ############### 4) Computing the rise and decay time ######################
    
    t_decay = Time[index_max_peak] - Time[index_peak]  
    			#time when the peak ends - time when the max value is reached
    t_rise = -Time[index_min_peak] + Time[index_peak]  #peak - start
		# -time when the peak starts + time when the max value is reached
    
    delta_t_rise_decay = np.sqrt(2) * delta_t  #error of both t_decay and t_rise
    
    
    ############# 5) Computing amplitude ##################

    amplitude = Amplitude(peak, baseline, delta_V)['amplitude[V]' ]
    delta_amplitude = Amplitude(peak, baseline, delta_V)['\Delta(amplitude[V])']
    

    ######################## 6)Plot #####################
    
    plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
    plt.plot( 1e6 * Time, 1e3 * Voltage, 'b.-')    
    plt.plot( 1e6 * time_peak , 1e3 * voltage_peak_neg, 'r.-')   
    plt.title("Waveform", fontsize=22)           #title
    plt.xlabel("time (us)", fontsize=14)                        #xlabel
    plt.ylabel("voltage (mV)", fontsize=14)              #ylabel
    plt.legend(['data', 'peak interval'], fontsize=14) 
    # Set size of tick labels.
    plt.tick_params(axis='both', labelsize=14)              #size of axis
    plt.grid(True) 
    #plt.savefig('Signal_LYSO_sum.png', format='png') 
    
    
   ########### 7) Return of values ############################
   
   #the values will be returned in a dictionary indicating what is each
   #value
    values = {'|peak|[V]' : peak_abs, 'baseline[V]' : baseline, 
              'N_peak' : len_peak, 'index_peak': index_peak,
              '\Delta(V[V])' : delta_V, '\Delta(t[s])' : delta_t,
              'integral[V*s]' : integral, '\Delta(integral[V*s])' : delta_integral,
              't_decay[s]' : t_decay, '\Delta(t_decay[s])' : delta_t_rise_decay,
              't_rise[s]' : t_rise, '\Delta(t_rise[s])' : delta_t_rise_decay,
              'amplitude[V]' : amplitude, '\Delta(amplitude[V])' : delta_amplitude}
              
    return values




###########################################################################
###########################################################################
#Function to compute the amplitude


def Amplitude(peak,baseline, delta_V):
    """Funciton to compute the amplitude of a peak from the oscilloscope
        *Inputs:
            .Peak: peak value [V] (with its sign)
            .Baseline: baseline value [V] (with its sign)
        *Outputs
            .Amplitude value [V]
            .\Delta(amplitude[V])
    """
    
    if peak > 0 and baseline >= 0:             #peak and baseline >0
        amp = peak-baseline
            
    elif peak < 0 and baseline <= 0:             #peak and baseline <0
        amp = np.absolute(peak) - np.absolute(baseline)
            
    elif peak > 0 and baseline <= 0:              #peak>0, amplitude <0
        amp = peak-baseline 
            
    elif peak < 0 and baseline >= 0:             #peak<0, baseline >0
        amp = np.absolute(peak) + baseline
    
    delta_amp = delta_V * np.sqrt(2)   #[V] error of the amplitude = peak +/- baseline
    
    #return of values
    
    values = {'amplitude[V]' : amp, '\Delta(amplitude[V])' : delta_amp}
    return values
