#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 10:20:36 2021

@author: dla

This function is to read .csv data from COMPASS (filtered data, which is the one 
with the good data), in order to capture the wave, or what we think it should be 
the wave. 

*Inputs:
        .name = filename in string format, eg: 'blablabla.csv'
        .row (optional) = the row to be used to extract the data. The .csv contains many rows, 
        	each row is a measurement. They are very similar between each other, 
        	so in general 0 is fine, that's why is the default value
        
*Outputs:
        .Channel = channel of the peak in the histogram
        .Board channel = board channel, the ch of the digitizer used (usually 15)
        .Samples = vector data with the waveform in strange units.
"""


#######0)General packages useful#################

import matplotlib.pyplot as plt  #for simplicity, to not write matplotlib.pyplot
        #everytime we want to plot something

#from scipy.stats import norm               ##norm.fit() fit to gaussian
import numpy as np
    #np contain linspaces as np.linspace(a,b,N)
import pandas as pd
        
from plotly.graph_objs import Bar, Layout
from plotly import offline
from scipy import stats as stats     #to find the most commom element in a list

import sys                   #to import functions from other folders!!
sys.path.insert(0, '/home/dla/Python/Functions_homemade')   #path where I have the functions
##########



#
############1) Function definition####################

def ReadCsvCOMPASS(name, row=0 ):
    
    load_csv = pd.read_csv('DataF_CH15@V1725S_646_run.csv', sep=';', header = None, skiprows=1)
    #skiprow = 1 to skip the first row, the header, which is extremely strange.
    #header:    #Board = np.array([])
                #Channel = np.array([])
                #Timetag =  np.array([])
                #Energy =  np.array([])
                #Energyshort =  np.array([])
                #Flags =  np.array([])
                #Samples =  np.array([])

    #Each sample is in one row. The samples goes from column 6 (cuntings starts at 0)
    #to the end

    size = load_csv.shape                               #size of the csv
    
    sample = np.array( [] )                            #storing of the wave samples

    for i in range(6,size[1]):
        sample = np.append(sample, load_csv[i][row])
        #storing the sample value of the given row, but for the different columns
    
    board_channel = load_csv[1][row]                 #Board channel
    timetag = load_csv[2][row]			     #[ps] Timestamp of the sample
    energy = load_csv[3][row]                       #Energy, in channels    
    
   ########### 6) Return of values ############################
   
   #the values will be returned in a dictionary indicating what is each
   #value
    values = {'Board_ch' : board_channel, 'E (ch)' : energy, 
              'Sample': sample, 'Timetag' : timetag
              }
              
    return values
