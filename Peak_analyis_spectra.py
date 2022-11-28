# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 11:05:28 2022

@author: lopedan
"""

#%% ################## 0) General useful packages ##################

import matplotlib.pyplot as plt  #for simplicity, to not write matplotlib.pyplot
        #everyx we want to plot something
import numpy as np          #np contain linspaces as np.linspace(a,b,N)
from scipy import stats    #to find the most commom element in a lis
import pandas as pd          #pandas, fundamental data science library
import sys                   #to import functions from other folders!!
sys.path.insert(0, '//net1.cec.eu.int/jrc-services/KRU-Users/lopedan/Desktop/PhD_Alemania-Sevilla_Residuos_nuclearesMariaVillaAlfageme/Python/Functions')   
                                    #path where I have the functions
import Fits

##################################



def Peak_analysis(x, y, index_min_peak, index_max_peak, index_df =0, signal_type = 'raw'):
    '''
    This function is to do the Pulse Shape Analysis (PSA) of a electric pulse (nuclear detectors)
from the oscilloscope.
    *Inputs:
        .x,y: x and y variables (list) that contains the data (np.array or pandas)
        .Signal_type = 'pre' or 'raw'. This indicated whether the peak is positive ('raw')
        	or negative ('pre'), which is neccessary to compute the amplitude
        .Indexes: indexes that limits the peak, eye spotted (will be plotted to confirm the selection)
    *Outputs:
        .Peak value (V)
        .Baseline (V)
        .delta_V and delta_t
        .Integral of the peak and its error
        .Number of elements of the peak
        .Rise and decay time of the peak, and their error (same for both)
    '''
    
    
################################################################################
    #%%############## 1) Formatting input ################
    '''
    If pandas is given, convert it into numpy, what we need (to avoid problems with
        indexes, since when you get a section of a df, the index remains, so the 1st value
        is not [0], but with np.array it is)
    '''
    if type(x) is not  np.ndarray:
        x = np.array(x)         #convert x into a np array
        
    if type(y) != np.ndarray:
        y = np.array(y)         #convert y into a np array
        
    #Note that both is not and != work similar. same as == and is 
    
    
    '''
     The peak data is:
    '''   
    y_peak = y[index_min_peak : index_max_peak]  
    x_peak = x[index_min_peak : index_max_peak]
    #y_peak = np.absolute(y_peak_neg)   #abs() because this
        #peak contains both > and ,0 values, so to add them in order to count
        #them, I have to put everything in the positive value
        
    len_peak =  len( y[index_min_peak-1:index_max_peak-1] )           #len of the peak   
    
    
    ################################################################################
    #%%############## 2) Finding the peak (max/min value) ################
    
    	#to find the peak I could easily do (absolute value to be able to detect positive peaks (raw signal) and 
    		#negative signals (from pre) ):
    if signal_type == 'pre':                   #signal from pre ==> negative peak
                
        peak = min( y )                           #[V] Peak value
        peak_abs = max( np.absolute(y) )          #[V] |Peak value|
    
    elif signal_type == 'raw':                  #signal from raw ==> positive peak
        peak = max( y )             #[V] |Peak value|
        peak_abs = max( np.absolute(y) )          #[V] |Peak value|
        
    
    index_peak = np.where( np.absolute(y) == peak_abs )[0][0]        #index
            #index_peak is a tuple, which contains a single np.aray, so with [0][0] you get the value
    #using the index of the peak I see the interval by looking at the .csv file
    #and the plot




################################################################################
    #%%################## 3) Finding the baseline #####################
    
        #to choose the baseline, the 1st approach was to do the mean between channel 0 and the
        #channel where the peak start to appear. But this is a ad idea becaue the error is very high,
        #so the 2nd approach will simply be choose the most frquent value

    #len_baseline_signal_points = len(ys[column_index][0:index_min_peak-1])  #length of the ys used
        #to compute the baseline. Will be neccesary for error calcs
    #baseline_signal = sum(ys[column_index][0:index_min_peak-1]) / len_baseline_signal_points    
        #[V] baseline, to compute the peak amplitude. I average between the initial value and
        #when the peak starts
        

    baseline = stats.mode(y)[0][0]       #[V] baseline y, the most common value
                    #[0]contains the values (array), [1] the frequency (array), so another [0]
                    #is needed to get the value
        
    #sum_y = sum(-y_peak + baseline)         #[V] total y of the peak, corrected with
                        #the baseline, so that this is the real y created!! Baseline-y because
                        #y of the amplitude is the greatest. DO NOT NEED IT ANYMORE!!!!



################################################################################
    #%% ###################### 4) Computing the integral of the peak ###################
    #The area under the peak (integral) can be computed numerically using the trapz method.
    #NEED TO THINK OF A WAY TO STIMATE ITS ERROR!!! I tried to put the error of the area of a rectange:
    # (Vmax-Vmin)*(tmax-tmin).


    integral = np.trapz(y_peak, x_peak)        #[V*t] area under the peak
    #delta_integral = 2 * np.sqrt( delta_V / (max(y_peak) - min(y_peak) ) 
                            # + delta_t/ (max(time_peak) - min(time_peak) ) )    #overstimation of the
                #error of the integral. This is the error of the 
                #area of the rectangle (Vmax-Vmin)*(tmax-tmin)


    #len_baseline_signal_points_stored.append(len_baseline_signal_points)


 ################################################################################   
    #%% ############ 5) Computing amplitude ##################

    amplitude = Amplitude(peak, baseline)['amplitude[V]' ]
    delta_amplitude = Amplitude(peak, baseline)['\Delta(amplitude[V])']
    
    
    
################################################################################    
    #%%########### 6) Gaussian fit ########################
    
    fit = Fits.Gaussian_fit(x_peak, y_peak) #dataframe containin the gausian fit
                #it also plot it. Note its index is 0, by defect
    
    
################################################################################    
    #%%####################### 7)Plot #####################
    
    plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
    plt.plot( x, y, 'b.-')    
    plt.plot( x_peak , y_peak, 'r.-')   
    plt.title("Peak", fontsize=22)           #title
    plt.xlabel("x", fontsize=14)                        #xlabel
    plt.ylabel("y ", fontsize=14)              #ylabel
    plt.legend(['data', 'peak interval'], fontsize=14) 
    # Set size of tick labels.
    plt.tick_params(axis='both', labelsize=14)              #size of axis
    plt.grid(True) 
    #plt.savefig('Signal_LYSO_sum.png', format='png') 
    
    
    
################################################################################    
   #%% ########## 8) Return of values ############################
   
   #the values will be returned in a dictionary indicating what is each
   #value
    aux = {'|peak|' : peak_abs, 'baseline' : baseline, 
              'N_peak' : len_peak, 'index_peak': index_peak,
              'integral' : integral, 
              'amplitude' : amplitude, '\Delta(amplitude)' : delta_amplitude,
              #Gaussian fit data
              'heigh' : fit['heigh'][0], '\Delta(heigh)' : fit['\Delta(heigh)'][0],
              'mean' : fit['mean'][0], '\Delta(mean)' : fit['\Delta(mean)'][0],
              'sigma' : fit['sigma'][0], '\Delta(sigma)' : fit['\Delta(sigma)'][0],
              'FWHM' : fit['FWHM'][0], '\Delta(FWHM)' : fit['\Delta(FWHM)'][0],
              'R[%]' : fit['R[%]'][0], '\Delta(R[%])' : fit['\Delta(R[%])'][0],
              }         #values to return
              
    values= pd.DataFrame(aux, index = [index_df] )  #df to return
    return values







################################################################################
################################################################################
#%% ##########################FUNCTIONS ################################################
###########################################################################
#Function to compute the amplitude


def Amplitude(peak,baseline, Delta_y = 0):
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
    
    delta_amp = Delta_y * np.sqrt(2)   #[V] error of the amplitude = peak +/- baseline
    
    #return of values
    
    values = {'amplitude[V]' : amp, '\Delta(amplitude[V])' : delta_amp}
    return values