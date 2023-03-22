# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 10:29:58 2023

@author: lopedan

This script will contain all the functions I create to read files from experimental measurements,
say TGA, XRD, etc (the ones that needed ofc)
"""

######################################
#%% ######### 0) General packages ###########
######################################


import matplotlib.pyplot as plt  #for simplicity, to not write matplotlib.pyplot
        #everytime we want to plot something
#from scipy.stats import norm               ##norm.fit() fit to gaussian
import numpy as np
    #np contain linspaces as np.linspace(a,b,N)
import pandas as pd
import sys                   #to import functions from other folders!!
sys.path.insert(0, '//net1.cec.eu.int/jrc-services/KRU-Users/lopedan/Desktop/PhD_Residuos_nucleares/Python/Functions')   
                                    #path where I have the functions
import Fits, Peak_analyis_spectra


#############################################################



######################################
#%% ######### ICPMS excel reader #############
#####################################

def Read_ICPMS_excel (name):
    '''
    Function that will read the excel file from ICPMS adn will return df with the relevant
    information, for easier handling /plotting. Note the excel should be a bit preprocessed:
            1) Including sample preparation
            2) Computing D_f * cps in a new sheet called 'Df_cps'
    
    Maybe that could be automatized? note that requires computing stuff from different sheets,
    and the D_f position could differ from file to file, so maybe more challenging that simply
    wworking with the excels a bit (Eww)
    
    *Inputs:
        .name: string with the name of the excel, without the .xlsx
        
    *Outputs:
        .several df with the Df_cps, %rsd, cps. Note that if I return X outputs, if I want
        to obtain a variable per output, in the script I should call X variables, like:
            a, b, .. = Read_ICPMS_excel(name)
            
            
    TO DO:
        .Include way to substract blank if desired (indicating which one is blank, etc)  
        .Include plotting, also sorting out what happens with the automatization of indexes (-5, -4,
                                                                                             etc)
        '''
    
    
    ########### 1) Raw load ###########
    '''
    Can be done easily with pandas. Since th excel sheet containing the cps and the excel file only
    differs in the .xlsx we can define the excel sheet name with the name given as input:
    '''
    excel_name = name + '.xlsx'
    
    #Load
    Dat_cps = pd.read_excel(excel_name, name, header = [1])
        #header 1 means take row 1 to give names to the columns
        #That contains the cps and cps*dil factor
        
    Dat_sig = pd.read_excel(excel_name, '%rsd', header = [1])
        #This contains the sigma values (measured 3 times, automatically computed average)
        
    Dat_cpsDf = pd.read_excel(excel_name, 'Df_cps', header = [1])
    
    '''
    Note there, the 1st row, the isotopes row, have no name, since stefaan do the excel in the 
    way he do it, so we need to set it manually:
    
    '''
    Dat_cps.rename(columns = {'Unnamed: 0' : 'Isotopes'}, inplace = True)   #changing column
                            #name from Unnamed to Isotopes
    Dat_sig.rename(columns = {'Unnamed: 0' : 'Isotopes'}, inplace = True)   
    Dat_cpsDf.rename(columns = {'Unnamed: 0' : 'Isotopes'}, inplace = True)  

    ############### 2) Clean df, 1 ############
    '''
    After the raw load, we can clean that a bit, creating a handful df, not the preovious, which
    are literally the excel in a df. That is, only collecting the relevant columns and putting them
    in a df, for further analysis (plotting, etc).
    
    Nevertheless, here could be relevant the fact of removing the blank or not, and for that you should
    say if you have a blank and if you want to delete it.
    
    The 1st cleaning is removing the 1st 4 rows, which contain bullshit, so we could do it with the 
    .drop method. Note the index are no longer used, simply erased.
    '''
    df_cps = Dat_cps.drop([0,1,2,3])    #Removing rows 01,2,3 (their index)
    df_cpsDf = Dat_cpsDf.drop([0,1,2,3])    #Removing rows 01,2,3 (their index)
    df_sig = Dat_sig.drop([0,1,2,3])    #Removing rows 01,2,3 (their index)
    
    '''
    A further step could be the blank removal, a bit trickier, but possibly could be though. 
    To do it in the future when I need it...
    '''
    
    
    ############## Ouptuts
    raw_df = {'cps' : Dat_cps, 'sigma' : Dat_sig, 'cpsDf' : Dat_cpsDf}    #raw stuff, the excel
                    #essentially, for debug
                    
    return df_cps, df_sig, df_cpsDf, raw_df
    
    
    
######################################
#%% ########## TGA reader ############ 
######################################


def Read_TGA (name):
    '''
    Function that reads the .txt file from TGA/DTA (HeHeHe), returning a dataframe with the
    relevant data (omitting intro).

    '''
    
    with open(name) as file_obj:
        lines = file_obj.readlines()
        #print('the number of lines of the document is',len(lines))
        '''
        From row 1 to 38 there is info about the run. Line 39 has the header, and then the data. Sio
        with a loop we could get the values
        '''
        T = np.array([])       #[°] storage of T values
        t = np.array([])       #[min] storage of t values
        DTA  = np.array([])       #[uV/mg] storage of TDA
        m  = np.array([])       #[%] storage of mass values
        sens  = np.array([])       #[uV/mW] storage of sensitivity values
        seg  = np.array([])       # storage of segment values

        #Now the loop to store the values:
        for i in range(39, len(lines)-1):   
                #-1 because last line is empty!
            T = np.append(T, float(lines[i].split(';')[0]) )
            t = np.append(t, float(lines[i].split(';')[1]) )
            DTA = np.append(DTA, float(lines[i].split(';')[2]) )
            m = np.append(m, float(lines[i].split(';')[3]) )
            sens = np.append(sens, float(lines[i].split(';')[4]) )
            seg = np.append(seg, float(lines[i].split(';')[5]) )

        #Now lets create the dictionary that will be converted to a df:
            dic = {'T[°]': T, 't[min]' : t, 'DTA[uV/mg]' : DTA, 'm[%]' : m, 
           'sens[uV/mW]' : sens, 'seg' : seg}
            df = pd.DataFrame.from_dict( dic )          #final dataframe
            
            
        #Finally, the return of the relevant values
        return df
    
    
 
######################################
#%% ########## XRD reader ############ 
######################################    
 
def Read_XRD_WB (name):
    '''
    Function that reads the .dat file from XRD in F141 (W Bonani), returning a df
    with the relevant info.
    
    Note the files give Q, I and Delta I (its error). Q is the module of the scattering vector,
    the vector that goes from the initial momentum vector to the final, being 2Theta the angle that
    covers. Then it can be proved that
    [https://physics.stackexchange.com/questions/123297/why-is-scattering-vector-vecq-called-vector-of-momentum-transfer#123300]

    Q = 4pi/lambda * sin (theta)
    
    But we want 2Theta so we will convert that by doing:
        Theta = asin (Q * lambda /4pi)
    Taking care since the trigonometric operations are in radians, so we need to convert the angle
    to degrees
    '''
    
    with open(name) as f_obj:
        lines = f_obj.readlines()       # to read all the lines, 1 by 1, and store them
        #
        '''
    The next step is to sepparate the lines and get the values. 1st 2 lines are info,
    so could be discarded. THe rest, if we split by ' ' (space) we have several blanks ('') 
    and the 3 numbers (the sepparator, the amount of spaces is variable, so we need to do that).
    With if statements we can select the numbers (non blanks) store the 3 of them, and we know 
    that in order they are Q, ExpIntensity and Experror.
        '''
    #Initialization variables to store
        Q = np.array( [] )
        I = np.array( [] )
        E = np.array( [] )
    
        for i in range(2, len(lines)):      #looping through all lines
            aux =  lines[i].split(' ')          #variable containing the line splitted
            #NOw we read element by element that file
            aux_2 = np.array([])            #varibale storing the 3 numbers
        
            for element in aux:
                if element != '':       #if element different from blank (''), do this
                    aux_2 = np.append(aux_2, float(element) )       #storing the values
        
            #Now I can store the values easily
            Q = np.append( Q, aux_2[0] )
            I = np.append( I, aux_2[1])
            E = np.append( E, aux_2[2] )
        
        #I can compute now 2theta:
        Lambda = 1.54184            #[Ang] wavelength of Kalpha of Cu
        Theta = np.arcsin(Lambda * Q / (4* np.pi )) * 180 / np.pi
        
        #Now lets create a dictionary to create a dataframe:
        dic = {'2Theta' : 2 * Theta, 'I' : I, 'Delta_I' : E}
        df = pd.DataFrame.from_dict( dic)
        #
        #Finally return the df:
        return df
