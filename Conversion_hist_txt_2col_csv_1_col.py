#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 08:10:46 2021

@author: dla

Conversion from .txt with 2 columns (ch and counts) to .csv with 1 columns (counts)
"""


import numpy as np
    #np contain linspaces as np.linspace(a,b,N)

import sys                   #to import functions from other folders!!
sys.path.insert(0, '/home/dla/Python/Functions_homemade')   #path where I have the functions

import Read_hist_txt


def Convert_txt_hist_2col_csv_1col(name,time, new_name):
    
    '''
    INPUTS:
            *name = file name, in the format 'whatever.txt'
            *time = measurement time
            *new_name is the .csv name, should be similar to the .txt, but with .csv. To be
    automatized..
    
    '''
    #name, time, new_name = 'Co60_BGO_3Thre_5000_8_histo_histo.txt', 300, 'Co60_BGO_3Thre_5000_8_histo_histo.csv'
                                #debug
    load = Read_hist_txt.Read_hist_txt_CAEN(name, time) 
   
    #Storing of the values
    ADC_channels = np.array( load[0] )
    counts_st = np.array( load[1] )
    rate_st = np.array( load[2] ) 
    
    #writting of the new file
    
    filename = new_name
    with open(filename, 'w') as file_object:
        for i in range(0,len(counts_st) ):
            
            if i< max(range(0, len(counts_st))): #for this case, a line break \n  needs to be added
                
                file_object.write(str(counts_st[i]) )
                file_object.write('\n')
            
            else:
                file_object.write(str(counts_st[i]) )
                
#Conversion:
    #Co done (testing)
    #Convert_txt_hist_2col_csv_1col('Ba133_BGO_3Thre_5000_8_histo_histo.txt',300, 'Ba133_BGO_3Thre_5000_8_histo_histo.csv')
    #Convert_txt_hist_2col_csv_1col('Cs137_BGO_3Thre_5000_8_histo_histo.txt',300, 'Cs137_BGO_3Thre_5000_8_histo_histo.csv')
    #Convert_txt_hist_2col_csv_1col('Na22_BGO_3Thre_5000_8_histo_histo.txt',300, 'Na22_BGO_3Thre_5000_8_histo_histo.csv')
    