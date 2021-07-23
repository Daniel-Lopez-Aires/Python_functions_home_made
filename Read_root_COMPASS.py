#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 09:23:48 2021

@author: Daniel López Aires// danlopair@gmail.com

"""


#######0) General packages useful#############33
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time     #to measure the time
#To import things from ROOT, previously you have to source thisroot.sh from 
    #the  command line! So, before opening anaconda to open spyder, you must do:
                #cd root/bin
                #source thisroot.sh
from ROOT import TCanvas, TFormula, TF1, TH1F, TH2F, TTree, TFile, TArrayS
from ROOT import gROOT, TPaveText, TGraphErrors, TGraph
from ROOT import TPad, TPaveLabel, TTreeReader

#############


#%%  ###############################################################
### 1) Function to read .root files containing a single channel (waveform included) ###
####################################################################


def ReadRootSingleCOMPASS(name, Waveform_saving = False):

    """
This function is to read .root data from COMPASS (filtered data, which is the 
one  with the good data) if 
	i) storing each channel in a single .root.
	ii) a single .root for all channels, but loading the file with a Tree 

Those files contains a Tree, with leafs and a TArrayS with the waveform

.root better than .csv because it is approx an order of magnitude lighter 
(10MB vs 100MB). Note that to obtain the wave, each channel have to be saved
in a .root, and the .root weights an order of magnitude more (8MB vs 200kB)

The data can be readed in 2 ways, using a function that will be removed called 
AsMatrix(), which is very fast (comparable to C++), and the other way, which
is the one found in the bilbio (example). The first takes about 0.5s, and the 
0ther 50s ==> 10 times shorter. So, will use it as long as it exists. When it
ceases to exist, the funciton automatically will sswitch to the other (try
    except blocks)

*Inputs:
        .name = filename in string format, eg: 
                'DataF_CH14@V1725S_646_run.root' / 'SDataF_run.root'
                                                    (S from time sorted)
        .Waveform_saving = if True, the outputs return the waveform. Default
            value = False
        
*Outputs:
        .A dictionary with 2 dafaframes:
            .Dataframe with:
                -Channel = channel of the peak in the histogram
                -Board channel = board channel, the ch of the digitizer used (usually 15)
                -Flags, Timestamp = vectors given by the.root, but useless for us
            .Dataframe with:  ONLY IF WAVEFORM_SAVING = TRUE
                - waveform = vector data with the waveform in strange units.
                -time = vector data for the waveform plot
        .Plot (python) of the 1D spectra and the 2D spectrum. Root plot, TH2F
                    computed but not plotted because it slows down a lot the PC
                    to play with that plot

#Biblio: https://root.cern.ch/doc/master/pyroot002__TTreeAsMatrix_8py.html
        https://root-forum.cern.ch/t/unable-to-read-a-ttree-with-leaves-and-a-branch-tarrays-with-pyroot/45578
        """ 


    data_root = TFile(name)                          #reading the .root file
    tree = data_root.Data_F                   #getting the Tree
    #Tree.Print()                            #print of the Tree
    n = tree.GetEntries()            #number of events (entries) of the Tree (numbers)

#Whatch out, that Tree contains leaves (Energy, Sample, etc), but it also 
#contain a branche, that inside it have more leaves.


######### 1.1)Read only specific branches #############
#There is a simple function called AsMatrix, which is incredibly faster, like
#20 times faster. However, it will be removed on root 6.26 (On July 2021, last
#version is root 6.24/00). So, will use it, but will implement a Try except bock,
#so that when that function do not exist, will run the other (more time 
#consuming) way, which is more similar to the C++ way.

    try:    #short way toe xtract the data, but which will be deleted, so that
        #when this do no exist, the except block will be executed
    
        dat = tree.AsMatrix(exclude=["Samples"])
            #load of all the data except the waveform.
        E = dat[:,3]
        board = dat[:,2]
        ch_digi = dat[:,0]
        timestamp = dat[:,1]
        flags = dat[:,5]
        
        #this way, the wave can not be loaded, so have to load the wave the old
        #way:
            
        if Waveform_saving:     #If Waveform_saving is True, then retrieve
            fN = np.array( [] )         #store fn, the number of points
                #of the waveform            
                
            for event in tree:  #for each event, store fN
                fN = np.append(fN, tree.Samples.fN)        
                
                
    except:                 #long way to extract the data
    
        E = np.array( [] )                    #store the Energy
        board = np.array( [] )                #store the Board
        ch_digi = np.array( [] )              #store the Channel
        timestamp = np.array( [] )            #store the Timestamp
        flags = np.array( [] )      #store the flags
        fN = np.array( [] )         #store fn, the number of points
                #of the waveform
               
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
    
    if Waveform_saving:     #If Waveform_saving is True, then retrieve
            #the waveform from the .root
            
        waveform = np.array( [] )    #store of the waveform 
    
    
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
   #The values will be returned in a dictionary. To return the values, pandas 
   #dataframe will be used.

    try:        #if using the AsMatrix version to load
        df_ch_timestamp = pd.DataFrame(data= dat,
                                   columns=['Ch digitizer', 'Timestamp[ps]', 'E[ch]', 'Board_ch', 'Flags'])
    
    except: #if loadind the data with the long way (loops)
        df_ch_timestamp = pd.DataFrame(data=np.array( [ch_digi,timestamp, E, board, flags] ).T,
                                   columns=['Ch digitizer', 'Timestamp[ps]', 'E[ch]', 'Board_ch', 'Flags'])
    
    
    if Waveform_saving:     #If waveform_saving is True, return the waveform too
        df_wave = pd.DataFrame(data=np.array( [time,waveform] ).T, columns=['Time[ns]', 'Voltage[ch]'])
        #the values will be returned in a dictionary indicating what is each
            #value
        values = {'Waveform' : df_wave,
              'Hist' : df_ch_timestamp}
        
    else:   #Waveform_saving = false, do not store it          
            values = {
              'Hist' : df_ch_timestamp}        
    
    return values          

#%%  ###############################################################
#### 2) Function to read .root files containing all the channels (only hist) ######
####################################################################

def ReadRootHistCOMPASS(name):

    """
This function is to read .root data from COMPASS (filtered data, which is the 
one  with the good data) if storing all the channels in a single .root. This reads
the .root that contains histograms, inside folders.

Actually, I do not use this, the other .root is the relevant one, but in any case here it
is this function just in case :))


.root better than .csv because it is approx an order of magnitude lighter 
(10MB vs 100MB). Note that to obtain the wave, each channel have to be saved
in a .root, and the .root weights an order of magnitude more (8MB vs 200kB).

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
    
    #Initialization
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




#%%  ###############################################################
####### 3) Function to make coincidences from the .root with all the channels######
####################################################################

def Coincidences_2ch_root(name, ch_A, ch_B, n_channels = 4096, gate = 3e5, 
                          save = True, debug = False, nbins = 70):

    """
This function is to make coincidences between 2 channels of MARS. It needs the
data from a single .root containng all the channels Time sorted, i.e., to load the
file 'SDataF_run.root'.

To do the coincidence, event by event, have to:
    1) Choose an event
    2) Check if the next one is from the other Channel or not
    3) If yes, coincidences are possible. Check if the time interval
           betweeen those events are small enough.
    4) Store the single energy values if the time interval is small enough,
            which will be the coincidences energies
            
The time interval to compare with its the gate, since once the gate is opened, 
the ADC records signals, so that, if the time interval between 2 measurements is
greater than the gate, those measurements were taken with 2 different gates and not
the same, and hence they are not in coincidence. Counterwise, if for a single gate
2 measurements were taken, those measurements are in coincidence.

@WATCH OUT:
    .To do the comparisons, I choose one event, and then I compare it with the next one,
        and I vary the original event. So, say event 1 is in coincidence with event 2.
        In the nex loop, it will choose event 2 and try t compare it with event 3.
        They won't be in coincidence, event 3 is from another gate (remember the gate
        is optimized to fit tightly the wave to avoid piling up events). 
        
            FIXED!!!

*Inputs:
        .name = filename in string format, eg:
                'SDataF_run.root'      (S from time sorted)
        .ch_A, ch_B = channels of the digitizer to be used to make the coincidences
        .gate [ps] = the time interval of the gate, which will be
            used to check whether 2 events are in coincidence or not. Default gate by
            COMPASS = 300ns = 3e5 ps
        .n_channels = number of channels of the ADC. Default: 4096
            NOTE that this could also be obtained from the other root, the one
            with the histograms, but to avoid loading 2 histos, could give it 
            as an input (in the load function from that file, it can be seen how
             the number of channels can be obtained). 
        .save = if True the plots with the single energies and the coincidences
            are saved. Default = True
        .debug: if debug = True, it plots the single energy that have 
            coincidences to tet the results obtained
        .nbins = number of bins for the single hist plots! DEfault value = 70
        
*Outputs:
        .Dictionary with:
            - Single energies values of both channels
            - Pandas dataframe containg the coincidence energies for 
            both channels (this is simply a subset of the single energies)
            - Dictionary with dataframes with the data from the .root file
        .Plots of single E spectras and 2D spectra in subplots.
        .Run time of the data loading and the coincidence making


@@@@@@@@@@ TO DO:
    1) TH2F plot do not generated if running this on terminal.

        """

#################0) Initialization ############
    t_begin = time.time()       #to measure time

    data  = ReadRootSingleCOMPASS(name)                     #data loading
    t_end_load = time.time()

    t_load = t_end_load - t_begin   #[s] time spent loading

    n_events = len(data['Hist'])   #number of events. Rows go from 0 to
            #n_events - 1


#################1) Single E extraction############
#This is very easy, just check the channel, and depending on that store the energy
    t_begin_coin = time.time()
#Initialization
    E_A = np.array( [] )                #Single energies of digi channel A
    E_B = np.array( [] )                #Single energies of digi channel B


    for i in range(0,n_events):  #events goes from 0 to n_events-1

        if data['Hist']['Ch digitizer'][i] == ch_A:  #If the event is ch 14
        
            E_A = np.append(E_A, data['Hist']['E[ch]'][i] )  #storing of the
                    #single energy
        elif data['Hist']['Ch digitizer'][i] == ch_B:  #If the event is ch 15
        
            E_B = np.append(E_B, data['Hist']['E[ch]'][i] )  #storing of the
                    #single energy  



#################2) Coincidences############
#To do the coincidence, event by event, have to:
    #1) Choose an event
    #2) Check if the next one is from the other Channel or not
    #3) If yes, coincidences are possible. Check if the time interval
           #betweeen those events are small enough.
    #4) Store the single energy values if the time interval is small enough,
            #which will be the coincidences energies

#A creation of a Th2F plot object (root) will also be implemented, since it is
#more faster than the matplotlib plot

    hist_2D = TH2F('2Dhist', '2Dhist', 4096, 0, 4096, 4096, 0, 4096)    #root hist


#Initialization
    E_A_c = np.array( [] )      #Energies of the ch A in coincidence with B
    E_B_c = np.array( [] )      #Energies of the ch B in coincidence with A
    
    
    for i in range(0,n_events-2):           #loop through all events
        
        #Time interval between the event i and the next one
        delta_t = data['Hist']['Timestamp[ps]'][i + 1] - data['Hist']['Timestamp[ps]'][i] 
                    #[ns] time interval between the events
                    
        if delta_t < gate: #if True, we have coincidences
            #now we have to find the channels of the events, to store it, or not
            #if the events are from the same channel
            
            if data['Hist']['Ch digitizer'][i] == ch_A:  #Check if the event is ch A
            
                if data['Hist']['Ch digitizer'][i + 1]== ch_B: #Check if the 
                #following row of the data correspond to the other 
                #channel (B)==> coincidence possible            
            
                    E_A_c = np.append(E_A_c, data['Hist']['E[ch]'][i] )
                    E_B_c = np.append(E_B_c, data['Hist']['E[ch]'][i+1] ) 
                    
                    hist_2D.Fill( data['Hist']['E[ch]'][i], 
                                 data['Hist']['E[ch]'][i+1] ) #fill of the hist
            
            else: #data['Hist']['Ch digitizer'][i] == ch_B:
                
                if data['Hist']['Ch digitizer'][i + 1]== ch_A: #Check if the 
                #following row of the data correspond to the other 
                #channel (A)==> coincidence possible  
                    E_A_c = np.append(E_A_c, data['Hist']['E[ch]'][i+1] )
                    E_B_c = np.append(E_B_c, data['Hist']['E[ch]'][i] )
                    
                    hist_2D.Fill( data['Hist']['E[ch]'][i+1], 
                                 data['Hist']['E[ch]'][i] ) #fill of the hist 


    time_end_coin = time.time()
    t_coin = time_end_coin - t_begin_coin       #[s] time spent in the comparison
    
    #The coincidences will be stores in a dataframe
    df_aux = pd.DataFrame(data=np.array( [E_A_c, E_B_c] ).T, 
                              columns=['E_coinc_ch'+str(ch_A), 
                                        'E_coinc_ch'+str(ch_B)] )
        #dataframe that contains the coincidence values.
        
    #To counts how many times do each pair appear, one can do:
    df_coinc_E = df_aux.groupby(df_aux.columns.tolist(),as_index=False).size()
    df_coinc_E.rename(columns = {'size':'Counts'}, inplace = True) #rename of
            #the new column to give it a good name
    
    #duplicate = df_coinc_E.duplicated()
   
    
   
    ########3) Plot ##############################3
    #Here both the single spectra and the 2D spectra will be plotted. 

    plt.figure(figsize=(16,12))  #width, heigh 6.4*4.8 inches by default
    #plt.suptitle("Spectra of the LED driver varying its amplitude", fontsize=22, wrap=True)           #title

    #1D spectra, ch A

    plt.subplot(1, 2, 1)
    plt.hist(E_A, bins = nbins)
    plt.title("Spectrum ch "+ str(ch_A), fontsize=22)           #title
    plt.xlabel("ADC Channels", fontsize=14)                        #xlabel
    plt.ylabel("Counts", fontsize=14)              #ylabel
    # Set size of tick labels.
    plt.tick_params(axis='both', labelsize=14)              #size of axis
    plt.grid(True) 
    #plt.xlim(0,n_channels)                       #limits of x axis

    
#1D spectra, ch B

    plt.subplot(1, 2, 2)
    plt.hist(E_B, bins = nbins)
    plt.title("Spectrum ch "+ str(ch_B), fontsize=22)           #title
    plt.xlabel("ADC Channels", fontsize=14)                        #xlabel
    plt.ylabel("Counts", fontsize=14)              #ylabel
    # Set size of tick labels.
    plt.tick_params(axis='both', labelsize=14)              #size of axis
    plt.grid(True) 
    #plt.xlim(0,n_channels)                       #limits of x axis
    
    if save:     #to save the plot    
        plt.savefig('Spectras_single_ch'+ str(ch_A) +
                '_ch'+ str(ch_B)+ '.png', format='png')        



    
    #2D spectra
    t_beggin_py_plot = time.time()
    
    plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
    plt.scatter(df_coinc_E['E_coinc_ch8'], df_coinc_E['E_coinc_ch11'], 
                df_coinc_E['Counts'], c=df_coinc_E.Counts)
    plt.title("2D spectrum, ch "+str(ch_A) + " and "+ str(ch_B) + " coincidence", fontsize=22, wrap=True)           #title
    plt.xlabel("E(ch) [ch" + str(ch_A) +"]", fontsize=14)                        #xlabel
    plt.ylabel("E(ch) [ch" + str(ch_B) +"]", fontsize=14)             #ylabel
    plt.tick_params(axis='both', labelsize=14)              #size of axis
    plt.grid(True) 
    cbar = plt.colorbar()                                      #Show the colorbar    
    cbar.set_label('Counts', fontsize=14)
    plt.xlim(0,n_channels)                       #limits of x axis
    plt.ylim(0,n_channels)                       #limits of y axis
        #both axis goes from 0 to the number of channels (each one to one digi 
            #channel)
    #plt.axis('scaled')
    
    t_plot_python_2D = time.time()- t_beggin_py_plot
    
    if save:     #to save the plot    
        plt.savefig('2D_spectrum_coinc_ch'+ str(ch_A) +
                '_ch'+ str(ch_B)+ '.png', format='png')      

    
    #########Debug plots###3
    #The plot of the energies values of each digi channels that is in 
    #coincidence with other energy value from the other channel will be plotted
    #as a debug method.
    if debug:
        
        plt.figure(figsize=(18,8))  #width, heigh 6.4*4.8 inches by default
        plt.suptitle("Spectra of the singles that have coincidences (subset of all the singles)",
                    fontsize=22, wrap=True)           #title


        plt.subplot(1, 2, 1)
        plt.hist(E_A_c, bins = nbins)
        plt.title("Spectrum (subset) ch " + str(ch_A), fontsize=20)           #title
        plt.xlabel("ADC Channels", fontsize=14)                        #xlabel
        plt.ylabel("Counts", fontsize=14)              #ylabel
    # Set size of tick labels.
        plt.tick_params(axis='both', labelsize=14)              #size of axis
        plt.grid(True) 
    #plt.xlim(0,n_channels)                       #limits of x axis
    

        plt.subplot(1, 2, 2)
        plt.hist(E_B_c, bins = nbins)
        plt.title("Spectrum (subset) ch "+ str(ch_B), fontsize=20)           #title
        plt.xlabel("ADC Channels", fontsize=14)                        #xlabel
        plt.ylabel("Counts", fontsize=14)              #ylabel
    # Set size of tick labels.
        plt.tick_params(axis='both', labelsize=14)              #size of axis
        plt.grid(True) 
    #plt.xlim(0,n_channels)                       #limits of x axis

        if save:     #to save the plot    
            plt.savefig('Spectras_debug_ch'+ str(ch_A) +
                '_ch'+ str(ch_B)+ '.png', format='png')      
    
    


    ########### 4) Return of values ############################
    #The values will be returned in a dictionary. To return the values, 
    #pandas dataframe will be used. Since E_A and E_B do not need to have the same
    #length, creating np.arrays with both of them do not work properly, so that
    #they will be stored via a dictionary

    #df_single_E = pd.DataFrame(data= np.array( [E_A, E_B] ).T,
                                    #columns=['E_single_ch', 'E_single_ch' ] )
                                   
    #the values will be returned in a dictionary indicating what is each
    #value
    values = {'E_single_ch'+str(ch_A) : E_A, 'E_single_ch'+str(ch_B) : E_B,
              'E_coinc' : df_coinc_E, 'Data_root' : data,
              'Run_time_load' : t_load , 'Run_time_coinc' : t_coin,
              #'Run_time_plot_TH2F' : t_root_plot, 
              'Run_time_plot_2D_python' : t_plot_python_2D}


    ############ 5) Plots, for command line run ############
    plt.show()  #call it at the end, to print all the plots
    #hist_2D.DrawCopy('colz')   #To plot the TH2F. Simpley Draw does not work
        #Since this slows down a lot th epc, although the generation is 
        #incredibly fast (1e-5s), will comment it, and use the python plot, 
        #which is virtually the same ;)
          
    return values   
    
    
    
   #####################GARBAGE###############
   
   
   
   #### 1) ALTERNATIVE MEHTOD TO MAKE HISTOGRAMS
   
   #To plot the single spectra with python, since we do not have counts, we have to do
    #the following. 
    
    #u, inv = np.unique(E_A, return_inverse=True)
    #counts = np.bincount(inv)

    #plt.subplot(1, 2, 1)
    #plt.bar(u, counts, edgecolor="black")
    #plt.title("Spectrum ch "+ str(ch_A), fontsize=22)           #title
    #plt.xlabel("ADC Channels", fontsize=14)                        #xlabel
    #plt.ylabel("Counts", fontsize=14)              #ylabel
    	# Set size of tick labels.
    #plt.tick_params(axis='both', labelsize=14)              #size of axis
    #plt.grid(True) 
    #plt.xlim(0,n_channels)                       #limits of x axis


