#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 09:23:48 2021

@author: Daniel López Aires// danlopair@gmail.com

"""


#######0) General packages useful#############33
import numpy as np

#To import things from ROOT, previously you have to source thisroot.sh from 
    #the  command line! So, before opening anaconda to open spyder, you must do:
                #cd root/bin
                #source thisroot.sh
from ROOT import TCanvas, TFormula, TF1, TH1F, TTree, TFile, TArrayS
from ROOT import gROOT, TPaveText, TGraphErrors, TGraph
from ROOT import TPad, TPaveLabel, TTreeReader

#############


#%%  ###############################################################
######## 1) Function to read .root files containing a single channels######
####################################################################


def ReadRootSingleCOMPASS(name):

    """
This function is to read .root data from COMPASS (filtered data, which is the 
one  with the good data) if storing each channel in a single .root.
.root better than .csv because it is approx an order of magnitude lighter 
(10MB vs 100MB). Note that to obtain the wave, each channel have to be saved
in a .root, and the .root weights an order of magnitude more (8MB vs 200kB)

*Inputs:
        .name = filename in string format, eg: 
                'DataF_CH14@V1725S_646_run.root'
        
*Outputs:
        .Channel = channel of the peak in the histogram
        .Board channel = board channel, the ch of the digitizer used (usually 15)
        .waveform = vector data with the waveform in strange units.
	.time = vector data for the waveform plot
	.Flags, Timestamp = vectors given by the.root, but useless for us
	
#Biblio: https://root.cern.ch/doc/master/pyroot002__TTreeAsMatrix_8py.html
        https://root-forum.cern.ch/t/unable-to-read-a-ttree-with-leaves-and-a-branch-tarrays-with-pyroot/45578
        """ 


    data_root = TFile(name)                          #reading the .root file
    tree=data_root.Data_F                   #getting the Tree
    #Tree.Print()                            #print of the Tree
    n = tree.GetEntries()            #number of entries of the Tree (numbers)

#Whatch out, that Tree contains leaves (Energy, Sample, etc), but it also 
#contain a branche, that inside it have more leaves.


######### 1.1)Read only specific branches #############
#There is a simple function called AsMatrix, but since it will be removed,
#will not use it. eg of that function:
        #E = Tree.AsMatrix(columns=["Energy"])
#Other version to read the values , similar to the one used in C== is:
    
    E = np.array( [] )                    #store the Energy
    board = np.array( [] )                #store the Board
    ch_digi = np.array( [] )              #store the Channel
    timestamp = np.array( [] )            #store the Timestamp
    flags = np.array( [] )      #store the flags
    fN = np.array( [] )         #store fn, the number of points
                #of the waveform
    waveform = np.array( [] )    #store of the waveform 
               
    for event in tree:  #for each event, store the desired values
        E = np.append(E, tree.Energy) 
        board = np.append(board, tree.Board)
        ch_digi = np.append(ch_digi, tree.Channel) 
        timestamp = np.append(timestamp, tree.Timestamp)
        flags = np.append(flags, tree.Flags)
        fN = np.append(fN, tree.Samples.fN)

#I will not loop through all the event because it will give 7321000 numbers,
#which is too much, and the pc takes a lot, so will choose the last event, which
#is simply using the variable event, since it will contain the last iteration of the
#previous loop. 


    for i in range(0,int(fN[-1]) ): #i goes from 0 to 1000=fN[j], for all j
        waveform = np.append(waveform,event.Samples.fArray[i])      
        #print(event.Samples.fArray[i])     
    

 # for elem in event.Samples.fArray: #for last event of the previous loop
 #        print(elem)
 #        Waveform = np.append(Waveform, elem)
#This loop do not end!!!  


    #the x axis for the waveform plot is determined by the ADC sampling time, which
    #is 4ns in our case [CAEN's COMPASS manual], so the X data will then be (the Y data
    #is the sample):      
    time = np.linspace(0, int(fN[-1])*4, int(fN[-1]) )          #[ns] time for
            #the waveform sample    
        

   ########### 1.2) Return of values ############################
   
   #the values will be returned in a dictionary indicating what is each
   #value
    values = {'Board_ch' : board, 'E[ch]' : E, 
              'Voltage_wave[ch]': waveform, 'Timestamp' : timestamp,
              'Ch_digi' : ch_digi, 'Flags': flags,
              'Time_wave[ns]': time
              }
              
    return values          

#%%  ###############################################################
####### 2) Function to read .root files containing all the channels######
####################################################################

def ReadRootFullCOMPASS(name):

    """
This function is to read .root data from COMPASS (filtered data, which is the 
one  with the good data) if storing all the channels in a single .root.
.root better than .csv because it is approx an order of magnitude lighter 
(10MB vs 100MB). Note that to obtain the wave, each channel have to be saved
in a .root, and the .root weights an order of magnitude more (8MB vs 200kB)

*Inputs:
        .name = filename in string format, eg: 
                'HcompassF_run_20210630_095907.root'
        
*Outputs:
        .Counts = counts of all the energy histograms, by columns:
                counts of hist_0; counts of hist-1; etc
        .Channels = vector array with the channels, from 0 to the max
        .n_Channels = number of channels
	
        """ 

    data_root = TFile(name)                          #reading the .root file
    E_folder = data_root.Energy        #folder that contains the 
                        #Energy histograms

    #Now, the 16 (from 0 to 15) histograms will be obtained, getting the
    #energy histogram data, the one we are interested in (for the moment).
    #To avoid importing 1 hist, and then copy paste that 16 times, will create
    #a loop. To do so, I need to change, in the loop, the name of the hist
    #to be obtained, which can be done easily:
    
    
    aux_1 = '_F_EnergyCH'                   #the hist name, for the loop,1
    aux_2 = '@V1725S_646'                   #the hist name, for the loop,2
    
    c = np.array( [] )          #Store of the Counts of the energy hist  
    
    for i in range(0,16): #loop for each hist, from 0 to 15
        
        H = E_folder.Get(aux_1+str(i)+aux_2)  #load of histogram
        n_ch = H.GetNcells()    #4097, the number of channels in the hist
        ch = np.linspace(1,n_ch, n_ch)      #channel linspace, for plotting   
    
        #To get the bin value of each bin, we need to do a loop, and for each 
        #iteration, use GetBinContent():
    

        
        #Since I do not know how to automatically store things in columns,
        #I could do the following, which works, first use append, then 
        #columnstack
        
        if i==0:        #1st way to store things (clear variable), 
                            #so using append
            for j in range(0, n_ch):   #loop through all the bins
            
                c = np.append(c,H.GetBinContent(j))
                
        else:   #rest of the cases, using columnstack
            c_aux = np.array( [] ) #auxiliar variable
            
            for j in range(0, n_ch):   #loop through all the bins
                c_aux = np.append(c_aux, H.GetBinContent(j))   
            
            c= np.column_stack((c, c_aux ))      #store in a column; once
            #we have computed the full new column with c_aux, store it

   ########### 2.2) Return of values ############################
   
   #the values will be returned in a dictionary indicating what is each
   #value
    values = {'Counts' : c, 'Channels' : ch, 
              'n_Channels': n_ch
              }
              
    return values   
