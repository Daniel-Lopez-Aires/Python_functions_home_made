#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 12:29:46 2021

@author: dla

Function to read histograms in .txt from CAEN's software SiPMkit_controlsoftware.
It needs the name in the proper format: "blablabla.txt". It also needs the measure time
to compute the count rate.
"""

#0) GENERAL PACKAGES, USEFULS
    
import numpy as np


def Read_hist_txt_CAEN(name, time):

    
    #1) 
    
    with open(name) as file_object:
            
        lines = file_object.readlines()
        print('the number of lines of the document is',len(lines))
        #This contains strings (have to be converted to numbers using int()
        #and \n, so the \n (salto de linea) have to be removed
        
        ADC_channels = np.array([] )                   #mid variable, to be used to count
        counts = np.array([] )                     #mid variable, to be used to count
        
        for i in range(len(lines)):        #get the values for each line    
            ADC_channels = np.append(ADC_channels, float(lines[i].split()[0]) )  
                    #store 1st number of the column
            counts = np.append(counts, float(lines[i].split()[1]) )            
                                        #store 2nd number of the column 
        
    return  np.array([ADC_channels, counts, counts/time]) #1st row time, 2nd row voltage


#Debug
#aaaa = Read_hist_txt_CAEN('Cs_137_1_elevacion_BGO_histo.txt',200)
