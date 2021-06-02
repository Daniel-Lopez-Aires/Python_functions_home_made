#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 09:52:43 2021

@author: dla

function to read .csv files from the emulator (1 column). It needs the name in the string format, eg,
"blablabla.CSV", and with it the csv is read, and the column (counts) is
stored in a numpy array. 

"""
    #0) GENERAL PACKAGES, USEFULS
    
import csv
import numpy as np
    
def Read_csv_emulator(name):

    #name = 'Co_60_120s_gain_0_3.csv' 			#Debug
    
    #1) 
    
    with open(name) as file_object:            #LYSO, raw   24/5 
   
        reader = csv.reader(file_object) #reader object assoaciated with the 
        #filerow_count = sum(1 for row in reader)            #number of rows
        #header_row = next(reader) #next return the next line. Since we only call it
        #once, we only get the 1st line
        #print(header_row)
        #n_columns = len(header_row)     #number of columns
    
        #Storing of the voltage (y) and time (x)
        counts = np.array([])
        for row in reader:
            counts = np.append(counts, float(row[0]) )
            
            
    return np.array([counts]) #1st row time, 2nd row voltage


#debug:
#aaaa = Read_csv_oscilloscope('Co_60_120s_gain_0_3.csv')    
