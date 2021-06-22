#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 11:39:55 2021

@author: dla

function to read .csv files from the oscilloscope. It needs the name in the string format, eg,
"blablabla.CSV", and with it the csv is read, and the 2 columns are
stored in a numpy array. 1st row contains the time, and 2nd row contains the voltage
"""

#0) GENERAL PACKAGES, USEFULS
    
import csv
import numpy as np


def Read_csv_oscilloscope(name):
    
    #1) 
    
    with open(name) as file_object:            #LYSO, raw   24/5 
   
        reader = csv.reader(file_object) #reader object assoaciated with the 
        #filerow_count = sum(1 for row in reader)            #number of rows
        #header_row = next(reader) #next return the next line. Since we only call it
        #once, we only get the 1st line
        #print(header_row)
        #n_columns = len(header_row)     #number of columns
    
        #Storing of the voltage (y) and time (x)
        time_help = np.array([])
        voltage_help = np.array([])
        for row in reader:
            voltage_help = np.append(voltage_help, float(row[-1 -1]))     #voltage; -1 = last , so -1 -1 is the
                            #second to last
            time_help = np.append(time_help, float(row[-1 -2]))     #time
            
            
    return np.array([time_help, voltage_help]) #1st row time, 2nd row voltage


#debug:
#aaaa = Read_csv_oscilloscope('TEK0002_LYSO_raw_24_5_both.CSV')           
