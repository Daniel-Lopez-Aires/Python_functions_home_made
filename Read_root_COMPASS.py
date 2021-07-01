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
        .n_events = number of events of each hist
        	
#Biblio: https://root.cern.ch/doc/master/pyroot002__TTreeAsMatrix_8py.html
        https://root-forum.cern.ch/t/unable-to-read-a-ttree-with-leaves-and-a-branch-tarrays-with-pyroot/45578
        """ 


    data_root = TFile(name)                          #reading the .root file
    tree=data_root.Data_F                   #getting the Tree
    #Tree.Print()                            #print of the Tree
    n = tree.GetEntries()            #number of events (entries) of the Tree (numbers)

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

#How to get the waveform (fArray)? Seeing the TBrwoser, it has n*1000 entries, which would
#mean that for each event, it stores 1000 numbers (1000=fN). I can not store all the
#n*1000 entries becasue the pc takes too long, so will only store the 1000 numbers
#for a single event, the last. This would is simply mean using the variable event,
# since it will contain the last iteration of the previous loop, and save the values.
#Note that:
	#- size(event.Samples.fArray)>= n*1000, so that each event has more than 1000 values,
	#  which is extremely range and do not match the hypothesis explained above. 
	#  If choosing the 1st thousand values, from 0 to 1000, the
	#  waveform is obtained, while if choosing other range, say from 1000 to 2000, etc, 
	#  weird results are obtained. Chosssing the value number n*1000 breaks spyder. 


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
    n_events = len(E) 				#number of events of the E hist   

   ########### 1.2) Return of values ############################
   
   #the values will be returned in a dictionary indicating what is each
   #value
    values = {'Board_ch' : board, 'E[ch]' : E, 
              'Voltage_wave[ch]': waveform, 'Timestamp' : timestamp,
              'Ch_digi' : ch_digi, 'Flags': flags,
              'Time_wave[ns]': time, 'n_events' : n_events
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
        .E[ch] = vector array with the channels, from 0 to the max
        .n_Channels = number of channels
        .n_events = number of events of each hist
	
        """ 

    data_root = TFile(name)                          #reading the .root file
    E_folder = data_root.Energy        #folder that contains the 
                        #Energy histograms

    #Now, the 16 (from 0 to 15) histograms will be obtained, getting the
    #energy histogram data, the one we are interested in (for the moment).
    #To avoid importing 1 hist, and then copy paste that 16 times, will create
    #a loop. To do so, I need to change, in the loop, the name of the hist
    #to be obtained, which can be done easily.
    #Note that:
    	#. There is a lot of bins, more than 2100000000, while
    	#the number of channels is 4096. No idea why this happens, but since I am
    	#not interested in such high energies (the bin center of the bin 2100000000
    	#is 2099999999.5), and the content is 0, of course.
    	#I will only choose the channels close to 4096 ch.
    
    
    aux_1 = '_F_EnergyCH'                   #the hist name, for the loop,1
    aux_2 = '@V1725S_646'                   #the hist name, for the loop,2
    
    #Initializaiton
    c = np.array( [] )          #Store of the Counts of the energy hist  
    n_events = np.array( [] )	#store the number of events (entries) on each hist
    ch_b = np.array( [] )         #store of the bin center of the energy hist
    
    for i in range(0,16): #loop for each hist, from 0 to 15
        
        H = E_folder.Get(aux_1+str(i)+aux_2)  #load of histogram
        n_ch = H.GetNcells()    #4097, the number of channels in the hist     
        n_events = np.append(n_events, H.GetEntries()) #number of events on the 
    	#histogram. Each hist could have its own n_entries.
    	
        #To get the bin value of each bin, we need to do a loop, and for each 
        #iteration, use GetBinContent(). The bin center can be obtained with
        #GetBinCenter()
 
        #Since I do not know how to automatically store things in columns,
        #I could do the following, which works, first use append, then 
        #columnstack
        
        if i==0:        #1st way to store things (clear variable), 
                            #so using append
            for j in range(0, n_ch):   #loop through all the bins
            
                c = np.append(c,H.GetBinContent(j))
                ch_b = np.append(ch_b,H.GetBinCenter(j))
        else:   #rest of the cases, using columnstack
            c_aux = np.array( [] ) #auxiliar variable to store the counts
            ch_aux = np.array( [] ) #auxiliar variable to store the ch
            
            for j in range(0, n_ch):   #loop through all the bins
                c_aux = np.append(c_aux, H.GetBinContent(j))   
                ch_aux = np.append(ch_aux, H.GetBinCenter(j))  
                
            c= np.column_stack((c, c_aux ))      #store in a column; once
            #we have computed the full new column with c_aux, store in a 
            #columns
            ch_b= np.column_stack((ch_b, ch_aux ))      #store in a column;     	
    
    ch = np.linspace(1,n_ch, n_ch)      #channel linspace, for plotting 
   
   
   ########### 2.2) Return of values ############################
   
   #the values will be returned in a dictionary indicating what is each
   #value
    values = {'Counts' : c, 'E[ch]' : ch, 
              'n_Channels': n_ch, 'n_events' : n_events,
              'Bin_center' : ch_b
              }
              
    return values   
