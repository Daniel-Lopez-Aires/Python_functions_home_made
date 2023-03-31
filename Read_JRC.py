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

def Read_ICPMS_excel (exc_name,D_f_data, sheet_name = 'Df_cps' ):
    '''
    Function that will read the excel file from ICPMS and will return df with the relevant
    information, for easier handling /plotting. Note the excel should be a bit preprocessed:
            1) Including sample preparation (neccesary to get the Dilution factor). In this excel,
               the sample names (1st colum where the Df is) must be the same as the names in the top
               of the ICPMS data (do it manually, lazy spaniard, less siesta and more work!)
            2) Clean sheet where only the intensity data is, to load it. THe name must be:
                "To_read"
            3) Prepare a similar sheet for the %rsd, called "%rsd_to_read"
    
    Maybe that could be automatized? note that requires computing stuff from different sheets,
    and the D_f position could differ from file to file, so maybe more challenging that simply
    wworking with the excels a bit (Eww)
    
    *Inputs:
        .exc_name: string with the name of the excel, without the .xlsx
        .sheet_name: string with the name of the sheet with the data to read (with counts, in the
            future maybe also concentration values?). Default value: 'Df_cps' (from acid vs no acid test)
        .D_f_data: array with the column number and the row interval in which that data is found.
            Df_data = [1, 3, 5] means from row 1 to 3 (included) and column 5 (E in letters)
        
    *Outputs:
        .several df with the Df_cps, %rsd, cps. Note that if I return X outputs, if I want
        to obtain a variable per output, in the script I should call X variables, like:
            a, b, .. = Read_ICPMS_excel(name)
            
            
    ######## TO DO ########
        .Include way to substract blank if desired (indicating which one is blank, etc) ? Note now
        you just read the sheet you prepared to read, maybe could be more optimized?
        .Include plotting, also sorting out what happens with the automatization of indexes (-5, -4,
                                                                                             etc)
        '''
    
    
    ########### 1) Raw load ###########
    '''
    Can be done easily with pandas. Since th excel sheet containing the cps and the excel file only
    differs in the .xlsx we can define the excel sheet name with the name given as input:
    '''
    excel_name = exc_name + '.xlsx'
    
    #Load
    Dat_cps = pd.read_excel(excel_name, sheet_name, header = [1])
        #header 1 means take row 1 to give names to the columns
        #That contains the cps and cps*dil factor
        
    Dat_rsd = pd.read_excel(excel_name, '%rsd_to_read', header = [1])
        #This contains the rsd values (measured 3 times, automatically computed average)
    Dat_sa_prep = pd.read_excel(excel_name, 'Sampl_prep', header = None)
                #Sample prep sheet. Ensure it has that name!!!!    
    
        
    '''
    From the sampl prep sheet I should get the dilutions factors, useful for correcting
    for it in both the RSD and in the cps
    '''
    D_f = Dat_sa_prep.iloc[D_f_data[0]-1 : D_f_data[1], D_f_data[2]-1 ]   #Dilution factor (pandas Series)
            #The -1 is because python start in 0 while excel in 1 for counting rows
    #I need to put the correct index names for the operations, that can be done like:
    D_f.index = Dat_sa_prep.iloc[D_f_data[0]-1 : D_f_data[1], 0 ]     #proper index name (to operate)
    
    
    ############### 2) Clean df, 1 ############
    
    '''
    Note there, the 1st row, the isotopes row, have no name, since stefaan do the excel in the 
    way he do it, so we need to change the name manually like:
    
    '''
    Dat_cps.rename(columns = {'Unnamed: 0' : 'Isotopes'}, inplace = True)   #changing column
                            #name from Unnamed to Isotopes
    Dat_rsd.rename(columns = {'Unnamed: 0' : 'Isotopes'}, inplace = True)   

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
    df_rsd = Dat_rsd.drop([0,1,2,3])    #Removing rows 01,2,3 (their index)
    
    '''
    A further step could be the blank removal, a bit trickier, but possibly could be though. 
    To do it in the future when I need it...
    '''
    
    '''
    Note the rsd file have the samples, then blank, then std, and then wash. 2 more columns in acid
    vs no acid than teh cps file. I should delete them in order to operate 
    (correct for the Df). For that, I could just create a new rsd df with the columns in common,
    since I know the name
    
    '''
    
    
    ################## 3) Derived calcs ###################################
    '''
    Note the sig is RSD = relative standard deviation. I should compute sigma (std dev),
    which could be plotted in the temporal plots for ex as error bar. Thats easy:
            RSD = sigma / <x> * 100 ==> sigma = RSD * <x> /100
    
    Note RSD I have in 1 df (for each measurement (column) I have lot of elements (rows) ),
    and 1 df is for RSD; the other is for the mean values, so I need to do operations between 
    them. Column 1 is the isotopes name so not to be sued, the rest can be used! 
    
    I can get all the columns but the first by doing:
        df_cps.iloc[:,1:]   gives everything but 1st column!
        
    Still, note the cps are multiplied by the dilution factor Df and applied some corrections
    (blanks, IS). i will forget about the 2nd things, more or less minor, and will consider the
    Df correction, which is essentially multiplying by it. So, I should multiply the RSD by the Df,
    and then apply that
    '''
    df_std= pd.DataFrame(df_cps.iloc[:,1:] * df_rsd.iloc[:,1:].values / 100, 
                  columns=df_cps.columns, index=df_cps.index)   #std df
    
    '''
    Now, I can apply the Dilution factor to both the std and the cps, should be straightforward, same
    fashion than above. Well, not as simple, since for the multiplication the indexes should be the same
    so, I redefined (above) the Df indexes so they matched the Df ones, and then that calc is straightforward
    '''
    
    df_stdDf = pd.DataFrame(df_std.iloc[:,1:] * D_f, 
                  columns=df_std.columns, index=df_std.index)           # std*Df df
    df_cpsDf = pd.DataFrame(df_cps.iloc[:,1:] * D_f, 
                  columns=df_cps.columns, index=df_std.index)           # std*Df df    
    
    '''
    TO end that, note the isotopes column non is NaN since the multiplications, so lets redefine it again:
    '''
    df_stdDf['Isotopes'] = Dat_cps['Isotopes']              #redefining the isotopes column
    df_cpsDf['Isotopes'] = Dat_cps['Isotopes']
    
    
    ############## 4) Ouptuts ######################
    
    debug_df = {'df_cps': df_cps,'df_std': df_std,
                'raw_cps' : Dat_cps, 'raw_rsd' : Dat_rsd }   #Dataframes for debug                
    
    return df_cpsDf, df_stdDf, D_f, debug_df                #return of values
    
    
    '''
    Sucessfully debbugged, enjoy bro!
    '''
    
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
