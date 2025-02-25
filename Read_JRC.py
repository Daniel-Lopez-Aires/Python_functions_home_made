# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 10:29:58 2023

@author: lopedan

This script will contain all the functions I create to read files 
from experimental measurements, say TGA, XRD, etc (the ones that needed ofc)

I will include in all the plots, minor grid lines, for better analysis from the
plot!

old: plt.grid(True)             #Only major grid

new: plt.minorticks_on()             #enabling minor grid lines
plt.grid(which = 'minor', linestyle=':', linewidth=0.5) 
                    #which both to plot major and minor grid lines
plt.grid(which = 'major')

For printing a break of line, \n goes at the end of the print, not at the beginning!
"""

#%%######### 0) General packages ###########
######################################


import matplotlib.pyplot as plt  #for simplicity, to not write matplotlib.pyplot
        #everytime we want to plot something
#from scipy.stats import norm               ##norm.fit() fit to gaussian
import numpy as np
    #np contain linspaces as np.linspace(a,b,N)
import pandas as pd
import os, sys                   #to import functions from other folders!!
sys.path.insert(0, 
 '//net1.cec.eu.int/jrc-services/KRU-Users/lopedan/Desktop/PhD_Residuos_nucleares/Python/Functions')   
                                    #path where I have the functions
sys.path.insert(0, '/home/dla/Python/Functions_homemade')
import Fits, Peak_analyis_spectra
import time as tr                                #to measure the running time
from scipy.optimize import curve_fit             #Fit tool
import warnings
warnings.filterwarnings('ignore')       #To ignore and not print warning
warnings.simplefilter("ignore")         #to ignore all warnings

#############################################################

#Useful stuff
Bent_color = {'Sard' : (.68,.24,.31), 'Tur' :  '#F6BE00', 'BK' : 'grey'} 
    #'Tur' :  '#EEE8AA' is perfect, but not for real visulaiztion xD
# Isot_rel = ['Si28', 'Al27', 'Mg24', 'Mn55', 'Fe56', 'Ca44', 'Na23', #bentonite elements
#             'Sr88', 'Cs133', 'Eu151', 'La139', 'U238']     #CL eleements
#             #Reserve: 'Ti46', 'Ti47', 'Ti48', 'Ti49', 'Ti50',
#             #Eu151 less abundant as Eu153, but Eu153 sufffer interferences from
#             #Ba oxides,so for low Eu concentrations, Eu151 better!! [Stef]        

Isot_rel = ['Si28(MR)', 'Al27(MR)', 'Mg24(MR)', 'Mn55(MR)', 'Fe56(MR)', 
            'Ca44(MR)', 'Na23(MR)', #bentonite elements
            'Sr88(LR)', 'Cs133(LR)', 'Eu151(LR)', 'Eu153(LR)',
            'La139(LR)', 'U238(LR)']     #CL eleements
            #Reserve: 'Ti46', 'Ti47', 'Ti48', 'Ti49', 'Ti50',
            #Eu151 less abundant as Eu153, but Eu153 sufffer interferences from
            #Ba oxides,so for low Eu concentrations, Eu151 better!! [Stef]   

"""
Isot rele Cs are from the Cs sep
"""
Font = 18               #Fontsize, for the plots (labels, ticks, legends, etc)           
            
#############################################################            
#%%## ## 1.1) ICPMS excel reader #############
#####################################

def Read_ICPMS_excel (excel_name, cps_sheet_name = 'To_read', 
                      return_debug = False):
    '''
    Function that will read the excel file from ICPMS and will return a df with 
    the relevant information, for easier handling /plotting. Note the excel 
    should be a bit preprocessed:
            
        1) Clean sheet where only the relevant data (cps, ppb, whatever) is, 
        to load it. You can remove ICPMS blanks (std 0, etc). THe Isotopes column
        in Column A in excel THe 1st isotope, Co59(LR) in row 7 in excel). 
        Sample names in row 2 in excel
            
    You could use this function to get the raw data (output from ICPMS) or to correct 
    them for the ICPMS dilution factor.     

    Note sometimes some random NaN data from excel can be added. Easy solution, 
    go to the excel sheet, and delete those rows/columns. No clue why this 
    happens, but that solves it (:   
                                                                                 
    Recommended to delete the "wash" sample, will only bring problems in analyis xD

    
    *Inputs:
        .excel_name: string with the name of the excel, with the .xlsx. note if 
        you select the file on the FIle viewer and do copy paste, the name will 
        be there. But dont forget the '', it must be an string! Eg: 'Excel.xlsx'
        .cps_sheet_name: string with the name of the sheet with the data to read 
            future maybe also concentration values?). Default value: 'To_read' 
        (from acid vs no acid test)
        .return_debug: if you want to get some extra df for debug 
        (raw data, without cleaning, so like the excel). Default value = False

        
    *Outputs:
        .several df with the cps/%rsd or whatver it is reading. Depending whether
        you want the debug you may  have 1 or 2 outputs. THe isotopes are the 
        index of the df, so the columns are the sample data! the column names
        are the sample names
            
    
    Note that if you have N outputs, if you want to obtain a variable per output, 
    in the script I should call X variables, like:
            a, b, ..n = Read_ICPMS_excel(name)
    If you write less, say 1, 2, some variable will contain more data, in a dictionary
            
            
    ######## TO DO ########
        . Arbitrary inputs (type *arb_arguments, with other stuff like, variable1,
            *aribtrary_arg) so you choose if you want rsd or no, to save some time?
            
    #######################
        '''
    
    
    ########### 1) Raw load ###########
    '''
    Can be done easily with pandas. Since th excel sheet containing the cps and
    the excel file only differs in the .xlsx we can define the excel sheet name 
    with the name given as input:
    '''    
    #Load
    Dat = pd.read_excel(excel_name, cps_sheet_name, header = [1], index_col=0)
        #header 1 means take row 1 to give names to the columns
        #That contains the cps and cps*dil factor
        #index col = 0 to use first column as index!
    '''
    Note once I suffered that the dimesions of those were not similar, and in 
    one sheet they were loadingn NaN values. I just erase those empty stuff in
    excel (selecting and delete) and after it worked!
    '''
    
    ############### 2) Clean df ############
    
    '''
    After the raw load, we can clean that a bit, creating a handful df, not the
    preovious, which are literally the excel in a df. That is, only collecting 
    the relevant columns and putting them in a df, for further analysis 
    (plotting, etc).
    
    Nevertheless, here could be relevant the fact of removing the blank or not, 
    and for that you should say if you have a blank and if you want to delete it.
    
    The 1st cleaning is removing the 1st 4 rows, which contain bullshit, so we 
    could do it with the  .drop method. Note the index are no longer used, simply 
    erased.
    '''
    #df_cps = Dat.drop(index = [Dat.index[0], Dat.index[1], Dat.index[2],Dat.index[3] ], axis = 0)   
         #Removing rows 0, 1,2,3 (their index)
         #that does not work if I load directly the Blk corr, so lets do it simpler:
    df_cps = Dat.iloc[4:,:] 
    
    #Another cleaning will be putting the df in numeric format. It is in object 
    #format, which gives problems
    #(YeroDivisionerror) with divisions, while for numbers there is no problem.

    df_cps = df_cps.apply(pd.to_numeric)    
                
    
  
    
    ############## 3) Ouptuts ######################
        
    if return_debug == True:    #return the debug
        #
        return df_cps, Dat              #return of values
    
    else: #no return debug, but yes corrections
            return df_cps

    '''
    Sucessfully debbugged, enjoy bro! Note I return the last the Df corrected because 
    I may not be interested in that, for ex if I want to plot the raw data.
    '''

    '''
    Manu read all the sheet of an excel like that:
    def readAllSheets(filename):
    if not os.path.isfile(filename):
        return None
    
    xls = pd.ExcelFile(filename)
    sheets = xls.sheet_names
    results = {}
    for sheet in sheets:
        results[sheet] = xls.parse(sheet)
        
    xls.close()
    
    return results, sheets
    '''



#%% ############## 1.1.5) Convert data from ppb to M ###########
#############################################################

def ICPMS_ppb_to_M(df_ppb, df_ppb_std, m_s = 1000, V_s =1, 
                   Delta_m_s = 0.0001, Delta_V_s = 0.001):
    '''
    This function will convert the ppb data (ng/g) from the ICPMS to the M (mol/L)
    data, the common format that papers use. The atomic weights values are
    taken from Stefaans excel, which come from an IUPAC publication!
    
    The data is from ppb, ng/gtot, we can convert it to M easily:
        10*-9 g/gtot * 1mol /At g * gtot/Vtot ==>
        ppb * 10**-9 * rho/At = M
    
    
    *Inputs
    .df_ppb,df_ppb_std: df with the ppb and std value
    .m_s, V_s: array with the mass and volumes of the samples. Default values:
        m_S = 1000g, V_s = 1L, to have desntiy 1kg/L
    .Delta_m_s, Delta_V_s: values of the uncertainties. Defaults:
            Delta_m_s = 0.0001g, Delta_V_s = 0.001L
            
    *Outputs
    .df with the concentration in M
    .df with the uncertainties in M, ASSUMING no uncertainty in the atomic weights
    '''
    
    
    ################ 3) Conversion from ppb to M #########3
    '''
    If I am reading a ppb data, I could convert it to Mol, I would only need
        1) Atomic weights
        2) Density of the samples, which I could have, since I measured both
    
    The operation is easy:
        ng/gtot *gtot/Vtot *mol/g = M ==> ppb * rho * At = M
    '''

    At_we = pd.read_excel('/home/dla/Python/at_wt_natural_elements_SVW.xlsx',
                          'To_read_atom_weight', index_col=0)   #read of the
                    #excel with the atomic weights
                    #
    
    rho = m_s/V_s               #density of the samples
    Delta_rho = rho*np.sqrt((Delta_m_s/m_s)**2 + (Delta_V_s/V_s)**2)
                #uncertainty of the density
                
    #Creation of empty df to store the data
    df_M = pd.DataFrame( index =df_ppb.index, columns = df_ppb.columns )
                    #empty df, but with defined columns and rows
    
    #Now to apply it I would need a loop
    for i in range(df_ppb.shape[0] ): 
        #loop through all isotopes, loop and stop when u see IS conc
        df_M.loc[df_ppb.index[i], :]= df_ppb.loc[df_ppb.index[i],
                df_ppb.columns[:]] * 10**-9 / At_we["atomic mass (u)"][
                At_we.index==df_ppb.index[i][:-4]].values * rho
        #print('Iter ' + str(i) + ' performed successfully')
    
    #It is not numeric, so lets convert it:
    df_M = df_M.apply(pd.to_numeric)  
    
    #The uncertainty calc can be outside the loop
    df_M_std = df_M * np.sqrt((Delta_rho/rho)**2 + (df_ppb_std/df_ppb)**2)
            #ASSUMING no error in the atomic weights!!!!!!!      


    ############# 2) Output ##########

    return df_M, df_M_std              
                    
                    
#%%########## 1.2) Future std computer!!!!!!!!!!!!!!!!!!!!!! #############
#####################################

def ICPMS_std_calculator (df_cps, df_rsd):
    '''
    Function that will compute the std from the %rsd data and the cps data. The data
    have already been read with the reader function.
    
    *Inputs:
        .df_cps: df containing the cps/ppb data. Isotopes are indexes, each 
        row a measurement, 1st 1st repl, then 2nd replicate, etc
	.df_rsd: df containing the rsd data, same format as the other

        
    *Outputs:
        .several df with the cps, %rsd, std, and a df series with the DIlution factor
            for the ICPMS sample prep. Depending whether you want the debug you may 
            have 2 or 3 outputs (always the raw and Df returned)
            
    Note that if you have N outputs, if you want to obtain a variable per output, 
    in the script I should call X variables, like:
            a, b, ..n = Read_ICPMS_excel(name)
    If you write less, say 1, 2, some variable will contain more data, in a dictionary
    
    ############ To Do ############
        Not tested, but was copied literally from a function that worked 
        (the reader was bigger), so it should work, right? not even changed 
        the names xD
    
    ##################
            
        '''
    
    
    ########### 1) Calc ###########
    '''
    Note the sig is RSD = relative standard deviation. I should compute sigma 
    (std dev), which could be plotted in the temporal plots for ex as error bar.
    Thats easy:
            RSD = sigma / <x> * 100 ==> sigma = RSD * <x> /100
    
    Note RSD I have in 1 df (for each measurement (column) I have lot of elements 
    (rows) ), and 1 df is for RSD; the other is for the mean values, so I need 
    to do operations between them. 
        
    Still, note the cps are multiplied by the dilution factor Df and applied 
    some corrections (blanks, IS). i will forget about the 2nd things, more or
    less minor, and will consider the Df correction, which is essentially
    multiplying by it. So, I should multiply the RSD by the Df, and then apply that
    
    '''
    df_std= pd.DataFrame(df_cps * df_rsd / 100, 
                  columns=df_cps.columns, index=df_cps.index)   #std df    
    
    '''
    TO apply the Df corrections or not we set a variable to choose, true or false
    '''    
        
    ############## 2) Ouptuts ######################

    return df_std

    


#%%######## 1.3) ICPMS Dilution factor finder #############
#####################################

def ICPMS_Df_finder (excel_name, D_f_data, samp_prep_sheet_name = 'Sample_prep'):
    '''
    Function that will find the ICPMS dilution factor from the excel containing
    the ICPMS sample preparation info.
    Note There is 2 dilution factors involved:
            1) Dlution actor for the ICPMS sample preparation (simeq 50). In this case you add 
                                    .8.8mL HNO3 1M
                                    .1mL IS (2IS; 0.5mL each)
                                    .0.2mL sample
            2) Dilution factor for the sample you use for the ICPMS sample prep (simeq 1). 
                    In this case you add some HNO3 conc (65% w/w)to the sample, to stabilize it.
                                                                            
    Note the excel should be a bit preprocessed:
            
        1) Including sample preparation (neccesary to get the Dilution factor). In this excel,
               the sample names must be the same as the names in the top
               of the ICPMS data (do it manually, lazy spaniard, less siesta and more work!). This sheet
               must be called 'Sample_prep'. Note the structure is like ICPMS results, column for sample.
               This is needed for further operations (corrections and so)
               
    *Inputs:
        .excel_name: string with the name of the excel, with the .xlsx. note if you select the file
        on the FIle viewer and do copy paste, the name will be there. But dont forget the '', it must
        be an string! Eg: 'Excel.xlsx'
        .D_f_data: array with the rows where the sample names and the Dilution factor 
            (ICPMS sample prep) is found.
            Note that info is also important for getting the index labels!
            df_samples = [1, 3] means names in row 1, Df in row 3 (in the excel)
        .samp_prep_sheet_name: string with the name of the sheet containing
            the ICPMS sample prep data, for the dilution factor. Default: 'Sample_prep'

    *Outputs:
        .ICPMS sample prep dilution factor (circa 50), in a pandas series format
        '''
    
    
    ########### 1) Read excel #################
    '''
    1st excel is read the excel
    '''
    Dat_sa_prep = pd.read_excel(excel_name, samp_prep_sheet_name, header = None)
                #Sample prep sheet. 
                
    ######## 2) Df extraction #######
    '''
    The Df extraction from the excel is simlpe:
    '''

    D_f = Dat_sa_prep.iloc[D_f_data[1]-1, 1: len(Dat_sa_prep)]   #Dilution factor (pandas Series)
            #The -1 is because python start in 0 while excel in 1 for counting rows and columns
    
    
    #I need to put the correct index names for the operations, that can be done like:
    D_f.index = Dat_sa_prep.iloc[D_f_data[0]-1, 1: len(Dat_sa_prep)]     #proper index name (to operate)
    
    
    ######## 3) cleaning ###############
    '''
    Since its object type, I will make it numeric, since everything will be easier
    with it (and strictly its true)
    '''
    D_f = D_f.apply(pd.to_numeric)                    #conversion to numeric data type
    
    ########### 4) Return #############
    '''
    Finally we return Df, which is a pandas serie
    '''
    
    return D_f            #return


#%% ##### 1.4) ICPMS Dilution factor corrector #############
#####################################

def ICPMS_Df_corrector (df_data, Df):
    '''
    Function that will apply the correction for the dilution factor to the ICPMS results.
    Note There is 2 dilution factors involved:
            1) Dlution factor for the ICPMS sample preparation (simeq 50). In this case you add 
                                    .8.8mL HNO3 1M
                                    .1mL IS (2IS; 0.5mL each)
                                    .0.2mL sample
            2) Dilution factor for the sample you use for the ICPMS sample prep (simeq 1). 
                    In this case you add some HNO3 conc (65% w/w)to the sample, to stabilize it.
                                                                            
                                                                            
    The correction is essentially scalatin for that factor, so the results takes into account
    that only a portion was measuring. So:
            df_data * Df (Df >=1)
    
    Its fundamental the labelling is appropiate! The data_exp sheet AND Sample_prep sheet must have
    same labels as the ICPMS output!!!!! Once the ICPMS resutls are there, change name to both, otherwise
    will return NaN!!!!

    *Inputs:
        .df_data: dataframe containing the cleaned data, to which the correction should be applied. 
		the isotopes are the index, so all columns are data!
        .D_f: pandas series containing the dilution factor to apply. Note the labelling is crutial,
            both inputs should have same labels (remember that you change the name in the exp sheet to
                                                 match the names that Stefaan used)

    *Outputs:
        .df with the correction factor (Df) applied

	CAn the [:,:] be removed ?? I would say yes, but check!!!
        '''
    
    
    ########### 1) Calcs ###########
    '''Now, I can apply the Dilution factor to both the std and the cps, should be straightforward, same
    fashion than above. Well, not as simple, since for the multiplication the indexes should be the same
    so, I redefined (above) the Df indexes so they matched the Df ones, and then that calc is straightforward
    '''
    
    df_corrected = pd.DataFrame(df_data * Df,
                                columns = df_data.columns, 
                                index = df_data.index)  #computing the correction
    
    #Just in case that df is not numeric (depends on Df mostly), I will convert
    #it to numeric
    
    df_corrected = df_corrected.apply(pd.to_numeric)
    
    ########### 2) Return #############
    return df_corrected             #return



#%% ######## 1.5) ICPMS Sample blank substraction #############
#####################################

def ICPMS_Sample_Blk_corrector (df_data, Nrepl = 2):
    '''
    Function that will apply the sample blank correction to the other samples in the ICPMS results df.
    This version is the 2 replicates version. Remember we already applied in the excel the IS correction
    and the ICPMS blanks corrections. Now its time for the blanks you did in your experiment.

    The correction is essentially substracting the blank (number 1) to the rest of the samples, but
    involved some operations since we have dataframes. So those will be here. Necessary that the data
    contain no Div0, ensure in the excel by deleting those!
    that only a portion was measuring. 
    
    Only for 2 and 3 replicates, but could be generalized for N replicates easily.


    *Inputs:
        .df_data: dataframe containing the data, the full data, with the 2 replicates. This is the 
        output fromt he reader function. You could apply this before or after Df corrections. Format:
            isotopes as index, columns the samples, 1st 1st replicate, then 2nd replicate. 2 replicates assume
            this function!!!!
        .Nrepl: number of replicates. Default value: 2 (2 replicates). Can also be 3

    *Outputs:
        .df with the correction factor (Df) applied
        '''
    
    
    ########### 1) Calcs ###########
    '''
    I must treat the 2 experiments are different, I should substract the blank 1 to the 1st emasurements
    and the 2 to the others. Since I ordered it in the right way (1st replicacte 1, then replicate 2, 
    I could) split it easily :D
            df.shape gives the shape of the df, n_rows, n_columns
    
    Note the df have number of samples * 2 replicates columns.
    
    Then, I will create a new dataframe substracting that data. To do so, I need to get rid
    of the isotopes column, since is text, and then add it again. Watch, the substraction is 
    easy with a pandas mehotd.

    I shuold then remove those columns
    from there, and replace negatives values for 0, for a good plot
    '''

    if Nrepl ==2:               #2 replicates, standard case
        df_1 = df_data.iloc[ :, 0: round( ( df_data.shape[1] ) / 2 ) ]      #1st replicate
        df_2 = df_data.iloc[ :, round( ( df_data.shape[1] ) / 2 ) :  ]       #replicate 2
    
        #The next step is blank substraction
        df_1_blk = df_1.subtract(df_1.iloc[:,0], axis = 0 ) 
                #0 since 1st columns is the blank  
        df_1_blk.drop( [df_1.iloc[:,0].name], axis = 1, inplace = True)     #removing the value I use to substract
    
    #Finally, for plotting purposes, we will replace negative values with 0:
        df_1_blk[df_1_blk < 0] = 0          #substituyin negative values with 0! 
                                    #This needs that no Div0 in excel!    
    #And now we do the same for replicate 2   
    
        df_2_blk = df_2.subtract(df_2.iloc[:,0],  axis = 0 )        #substraction
        df_2_blk.drop( [df_2.iloc[:,0].name], axis = 1, inplace = True)     #removing the value I use to substract
        df_2_blk[df_2_blk < 0] = 0      #replcaing negative values with 0
            #Those 3 lines would be needed for a loop, so sohuld be easy, if needed
    
    #Finally, lets store it

        df_blk = pd.concat( [df_1_blk, df_2_blk], axis = 1)         #mergind the 2 little df ina  huge one
    
    elif Nrepl == 3:         #3 replicates case
        #Copy paste the before code, now for 3 replicates
        df_1 = df_data.iloc[:, : round(df_data.shape[1] / 3)]
        df_2 = df_data.iloc[:, round(df_data.shape[1] / 3): 2*round(df_data.shape[1] / 3)]
        df_3 = df_data.iloc[:, 2*round(df_data.shape[1] / 3) :]

        df_1_blk = df_1.subtract(df_1.iloc[:,0], axis = 0 ) 
        df_2_blk = df_2.subtract(df_2.iloc[:,0], axis = 0 ) 
        df_3_blk = df_3.subtract(df_3.iloc[:,0], axis = 0 ) 

        df_1_blk.drop( [df_1.iloc[:,0].name], axis = 1, inplace = True)
        df_2_blk.drop( [df_2.iloc[:,0].name], axis = 1, inplace = True)
        df_3_blk.drop( [df_3.iloc[:,0].name], axis = 1, inplace = True)
        
        df_1_blk[df_1_blk < 0] = 0                      #replcaing negative values with 0
        df_2_blk[df_2_blk < 0] = 0     
        df_3_blk[df_3_blk < 0] = 0      
        
        df_blk = pd.concat( [df_1_blk, df_2_blk, df_3_blk], axis = 1)         #mergind the 2 little df ina  huge one
    
    else:                   #Otherwise, error
        print('Wrong number of replicates Nrepl!, nothing has been done!')
        df_blk = 0                      #Error case, not doing anything
    
    
    ########### 2) Return #############
    return df_blk             #return



#%%######## 1.6) ICPMS: Get Mass number #############
#####################################

def Get_A_Resol(isotope_name):
    '''
    Function to get the mass number A and the resolution type from the isotope name of the ICPMS excel. 
    The name is in the format:
            .BB111(MR)
    where BB can be 1 or 2 leters indicating the chemical symbol, and 111 can be 2 numbers also, the
    mass number. (MR) or (LR) is low or high resolution. There are 2 excepctions, Ar40Ar40(LR/MR)
    and U238O16(LR/MR).
    
    * Input
        .isotope_name: string with the ICPMS excel format, eg 'U238(LR)', 'Co59(MR)', etc
        
    #Output
        .Mass number A of the isotope. For U238(LR) 238, for Co59(MR) 59, etc  
        .Resolution format: LR or MR
    '''
    
    '''
    1st thing to do is split the (LR) and the other thing, can be done easily:
    '''
    
    aux = isotope_name.split('(')     #This gives 2 strings in an array, [ BB111,'MR)']
    resol = aux[1][:-1]
    ''' 
    Now in the first eleemnt of that I need to check:
        1) Name with 1 or 2 letters? with isalpha() you get True if letter. isdigit() for numbers)
    I need to introduce the 2 expcetions as different cases also. All in a loop, since if i get the 
    exceptions, I do not want to continue
    '''
    
    if aux[0] == 'U238O16':         #exception
        mass = 238
    
    elif aux[0] == 'Ar40Ar40':     #exception 2
        mass = 80
    
    
    elif aux[0][1].isalpha():     #If 2nd letter is an alpha ==> name with 2 letters
        mass = aux[0][2:]
    else:           # 2nd letter is a number, so name with 1 letter (eg, U)
        mass = aux[0][1:]
        
    return mass, resol             #mass is an string!




#%%######## 1.7) ICPMS IS sens calculation/plotter #############
#####################################

def IS_sens_calculator_plotter(df_cps_ppb, df_std,
        IS_meas = ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)'],
        name_IS_sens_LR_plot = 'IS_sensLR_plot', 
        name_IS_sens_MR_plot = 'IS_sensMR_plot'):
    '''
    Function part of the ICPMS Data processing! Function that will compute the IS sens (cps/ppb) 
    for a df containing the cps and ppb, and another df with its std. The format is like
    Stefaans excel, his first sheet. It IS needed the ppb table below the cps data. How below? Do not care, I
    find the data using locate functions ;)
    Note std is computed using the squared error propagation methods!  

    
    !!!!!!ASSUMPTION!!!! : Delta(ppb IS)/ppb IS = 1%, completely arbitrary, to simplify everything, as a 1st approx!
    it may be changed after! Those ppb come from perkin Elmer multielementa solutions, and derived calcs 
        (logbooks, calculations of concentrations, etc)
        !!!!!!!!!!!!!!!!!!
    
    Note the plots appear in the plots pannel in Spyder, but are also saved ;)

    * Input
        .df_cps_ppb: df containing the cps and ppb data. THe fashion is like Stefaans raw sheet. 
        Isotopes are index. Take care of the names of the columns (like in IS conc table or so, the values you find),
        if they dont exist will give error! For the ppb data, the names are just Co-59, Ho-165, In-115, Th-232.
        This is extremely important. Otherwise it will not work!
        . df_std: df containng the std of the cps. Default value: None, so I dont need to give it, since in the
            past I was not giving it, not to obtian error with the old scripts
        .IS_meas: array containing in a list the measured Internal Standards, containing its resolution, like 
            how they appear in the isotopes column. Default value:
                ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)']
            That mean those isotopes were measured. If Ho165(MR) also measured, just included it, and fine ;)
        .name_IS_sens_LR_plot: name for the plot of the IS sens for LR case. Similar for MR. Default values:
            'IS_sensLR_plot' and 'IS_sensMR_plot'. To that is added the .png to create the file.
              
    #Output
        .df with the IS sens data
        .df with the std of the IS sens.
        .Plot with the IS sens plot, 2 plots, 1 for LR and other for MR, in .png file
    
    Note that if this only one argument for output written, then the output will be a tuple with 2 df        
        
    ###### TO Do
        .Function to compute the ppb data, from Sum isobars? May be complex, possibly too time consuming ?
        . Think a way so that std can be optionally given, so my old scripts are still working? The alternative
        is modifying them, inclduding the std, will not be fatal though xD
    '''
    
    # if df_std == 'No':      #if no std passed
    #     df_std = df_cps_ppb * 0 #Creating a ceros matrix if no error provided, so we avoid doing an if loop
    
    
    ############ Data finder ################
    '''
    1st thing to find the ppb data. Should be somehow below the cps data, which 
    has Stefaans format.
    To do the find in a loop way I define arrays which the elements to find. 
    Given as an input now, easier, KISS!

    Note they are associated, you take i element of both arrays, they go in pairs. 
    I ahve seen a way, is to create in the loop the lines, and add them separately. 
    The loop find the elements and perform the division cps/ppb
    
    !
    The dataframes will be out of the for loop, for coding efficiency!
    https://stackoverflow.com/questions/36489576/why-does-concatenation-of-dataframes-get-exponentially-slower
    '''

    df_IS_sens = pd.DataFrame()         #Empty df to store the values
    df_IS_sens_std = pd.DataFrame()        
    list_aux = np.array([])
    list_aux2 = np.array([])
    
    for i in range(0, len(IS_meas)):
        value_to_find = IS_meas[i]
        
        ##### Getting the ppb values #########
        '''
        This will be made from the IS_meas variable. Since each element has a 
        different letter, I can just read the 1st letter and say:
            C ==> Co-59 in ppb
            I ==> In115 in ppb
            etc
        '''
        if value_to_find[0] == 'C':  #Co59 case
            value_to_find_ppb = 'Co-59'
        elif value_to_find[0] == 'I':  #In115 case
            value_to_find_ppb = 'In-115'
        elif value_to_find[0] == 'H':  #Ho165 case
            value_to_find_ppb = 'Ho-165'
        elif value_to_find[0] == 'T':  #Th232 case
            value_to_find_ppb = 'Th-232'
            
        cps_IS = df_cps_ppb.loc[df_cps_ppb.index == value_to_find]      #cps of the IS (give full row)
        std_IS = df_std.loc[df_std.index == value_to_find]        #%rstd of the IS 
        ppb_IS = df_cps_ppb.loc[df_cps_ppb.index == value_to_find_ppb]  #ppb of the IS

        #Now the operations, element wise, in array mode to avoid the index problems
        aux = cps_IS.iloc[:,:].values / ppb_IS.iloc[:,:].values
            #eleement wise operation, like this you dont care about having diferent 
            #indexes, since you are multiplying arrays. I erased the 1st
            #value (Co name), so I needit to add it again
        
        #To compute the error of the sens, I assume that Delta(ppb) = 1% ppb ==> 
                    #Delta(ppb)/ppb = 1/100!!
        aux2 = aux * np.sqrt((std_IS/cps_IS)**2 + (1/100)**2)     #std values
        
        #TO store temporarily those values I create an auxiliary list
        #     #no df, not efficient!
        # list_aux = np.append(list_aux, aux) #append it (in rows not possible)
        # list_aux2 = np.append(list_aux2, aux2)
        
        df_aux = pd.DataFrame(data = aux)
        df_aux2 =pd.DataFrame(data = aux2)
        
        #And I add that to the storing df
        df_IS_sens = pd.concat([df_IS_sens,df_aux], ignore_index= True)
        df_IS_sens_std = pd.concat([df_IS_sens_std, df_aux2], ignore_index= True)
    
    '''
    Now, to add the isotopes as index, we first add them as a column, and then we set
    the column to index:we need to insert the isotopes column, giving as values the isotopes names or so:
    '''
    #df_IS_sens['Isotopes'] = ['Co LR', 'In LR', 'Ho LR', 'Th LR', 'Co MR', 'In MR']
    df_IS_sens['Isotopes'] = IS_meas        #setting isotopes name from the loop variable, better
    df_IS_sens.set_index('Isotopes', inplace = True)
        
    df_IS_sens_std['Isotopes'] = IS_meas
    df_IS_sens_std.set_index('Isotopes', inplace = True)

    '''
    As a final stylish, I can put the same column names as the df used to create the IS df:
    '''
    df_IS_sens.columns = df_cps_ppb.columns
    df_IS_sens_std.columns = df_cps_ppb.columns
    
    
    'I can print the %rstd values = std/mean * 100 of the IS sens, useful to spot fluctuations!'
    rstd = df_IS_sens.std(axis =1) / df_IS_sens.mean(axis =1) * 100                   #%rstd of the IS sens
                                    #ofc agrees with excel!
    print('##################################################################### \n')
    print('%rstd of the IS sens:')
    print(rstd)
    print('Values >= 5/6% start to be suspicious, something happened! (Th oxidation for ex?) \n')
    print('##################################################################### \n')
    
    
    #################################################
    ############# IS sens plotter ####################
    ########
    '''
    I need to do 2 plots, 1 for the LR and other for the MR. I should sort them somehow. I found how xD, explicit plotting,
    recommended by python, better than implicit(what you normally do xD)
    '''
    
    pltL = plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default  LR plot!
    axL = pltL.subplots()
    axL.set_title("IS sens LR along measuring sequence", fontsize=22, wrap=True)           #title
   
    pltM = plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default   MR plot!
    axM = pltM.subplots()
    axM.set_title("IS sens MR along measuring sequence", fontsize=22, wrap=True)           #title    
    
    for i in range(0, df_IS_sens.shape[0]):   #looping through rows of the df
        if df_IS_sens.index[i][-4:] == '(LR)':      #low resolution
            axL.plot(list(range(0, df_IS_sens.shape[1])), df_IS_sens.iloc[i,:],'-o' ,label = df_IS_sens.index[i]) 
            
        else:   #MR
            axM.plot(list(range(0, df_IS_sens.shape[1])), df_IS_sens.iloc[i,:],'-o' ,label = df_IS_sens.index[i])  
    
    #Final styling of the plot
    axL.legend(fontsize = Font)
    axM.legend(fontsize = Font)
    axL.grid(True)
    axM.grid(True)
    axL.set_xlabel('Sample number', size = Font)
    axM.set_xlabel('Sample number', size = Font)
    axL.set_ylabel("cps/ppb", size = Font)
    axM.set_ylabel("cps/ppb", size = Font)   
    axL.tick_params(axis='both', labelsize= Font)              #size of axis
    axM.tick_params(axis='both', labelsize= Font)              #size of axis
    pltL.savefig(name_IS_sens_LR_plot + '.png', format='png', bbox_inches='tight')    #note I call plt, not ax!
    pltM.savefig(name_IS_sens_MR_plot + '.png', format='png', bbox_inches='tight')
    
         
    ############## Return of values ##########################
    return df_IS_sens, df_IS_sens_std            #mass is an string!


#%% ########## 1.8) ICPMS IS sens correction #############
#####################################

def IS_sens_correction(df_raw, df_raw_std, df_IS_sens, df_IS_sens_std,
                       IS_meas = ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)']):
    '''
    PART OF THE ICPMS Data processing!
    
    Function that will apply the IS sens correction to the cps data, from the raw ICPMS data.
    This is aprt of the ICPMS data analysis. You need to correct for IS sens, then apply ICPMS
    blanks, and then calibrate with 2+4 and 3+5 IS.
    
    'IS conc [ppb]' should be the box that indicates where the ppb chart start. That is, in A column in excel,
    that should be written, since that is used to find the ppb data!!!
    
    Note the plots appear in the plots pannel in Spyder, and that the IS measured decide the correction, I defined
    them with if statements, have not yet come up with a better idea...

    
    * Input
        .df_raw: df containing the raw cps and ppb data. THe fashion is like Stefaans raw sheet
        .df_raw_std: df containing the std from the cps data. 
        .df_IS_sens: df containing the IS sens data (cps/ppb). Ensure the size andindexing are
        appropiates!
        .df_IS_sens_std: df containing the IS sens data std (cps/ppb)
        .IS_meas: array containing in a list the measured Internal Standards, containing its resolution, like 
            how they appear in the isotopes column. Default value:
                ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)']
        
    #Output
        .df with the IS corrected data, containing the cps and ppb data
        .df with the std of the IS corrected data
        
        
    ###### TO Do #############
        .Switch cases if other corrections needed (avoid one IS because the measuring was wrong?)
    #########################
    '''
    
    ############## 1) Calcs #########################
    '''
    This is simple. If the 4 IS are fine, I do the correction in the following fashion:
            .From 59 to 80 apply Co59
            .From 81 to 138 apply In115
            .From 139 to 209 Ho165
            .From 210 to 248 Th232

    If there are less IS, then I will make less divisions. Ex, if for MR there is only Co and In,
            .Co from 59 to 80
            .Rest In115
            
    I shold be able to detech the name of th eisotope (row), get number, and apply those limits.
    The correction is data / IS sens * <IS sens>, being sens = cps/ppb

    The name format is:
        .A11(LR)
        .AB11(LR)
        .A111(LR)
        .AB111(MR)
        .And some random elements, like Ar40Ar40(L/MR), U238O16(L/MR)

    So, avoiding the weird stuff, the chemical name can have 1 or 2 letters, and the number can be 2 or 3 ciphers.
    I created a function to get the mass number. Note the operations should be
    line by line. The structure would be loop like, at the beginning I check the mass number, and apply. Lets do a
    single one
    
    Now I can try to do a loop. I could give where the IS sens start from the df inspection (240), and 
    I know the order:
        .Co LR
        .In LR
        .Ho R
        . Th LR
        .CO MR
        . In MR
        
    Next step is to create a function to get the mass and apply one or other correction. We already have the 
    function, so now we need to create a loop (after a function) that get the mass, and
    apply one or other correction
    '''   
    
    df_IS_co = df_raw.copy()       #Thats the proper way to create it, copy so it is not asigned to it!
    df_IS_co_std = df_raw.copy() 
    
    	#The loop should go until the alst isotope, which I can find by finding IS conc ppb, the first thing 
        #for the ppb chart! so, THIS data is mandatory that exist like that!!
    '''
    To generalize, I will do if statement for the different IS measured cases. By far only 2, the sequence 
    done in the past, the 4 in LR and 2 in MR; and now (8/23) 4 in MR also, the most detailed case. 11/23,
    only LR elements!
    '''
    
    if IS_meas == ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)']:   
                                                                    ########case 1, the old scenario (exp bef 8/23)
     
        for i in range(np.where(df_raw.index == 'IS conc [ppb]')[0][0]): #loop through all isotopes
                #df_raw.loc[229,'Isotopes'] = df_raw.iloc[226,0]#relation iloc, loc
        
            mass, resol = Get_A_Resol(df_raw.index[i])
        #No we neeed to see the resolution, since for LR we have 4, and for MR only 2
            if resol == 'LR':           #low resolution
            #Now, which mass?
                if int(mass) in range(59,81): #mass from 59 to 80 included both ==> Use Co59 (element 0 in IS sens)
                    df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],
                            df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[0], df_IS_sens.columns[:] 
                            ] * df_IS_sens.loc[df_IS_sens.index[0], df_IS_sens.columns[:]].mean() 
                    #Take a look a that, I can split the line makins advantage of the parenthesis existing
                    df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[0], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[0], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[0],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[0], df_IS_sens.columns[:]].mean() )**2
                                )                   #quadratic error prop!                   
                elif int(mass) in range (81, 139):
                    df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],
                        df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[1], df_IS_sens.columns[:] 
                        ] * df_IS_sens.loc[df_IS_sens.index[1], df_IS_sens.columns[:]].mean() 
                    df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[1], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[1], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[1],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[1], df_IS_sens.columns[:]].mean() )**2
                                )
                elif int(mass) in range (139, 210):
                    df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],
                            df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[2],df_IS_sens.columns[:] 
                            ] * df_IS_sens.loc[df_IS_sens.index[2],df_IS_sens.columns[:]].mean() 
                    df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[2], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[2], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[2],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[2], df_IS_sens.columns[:]].mean() )**2
                                )
                else:   #if mass in range(210, 249):    rest of mass range
                    df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i], 
                        df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[3], df_IS_sens.columns[:] 
                        ] * df_IS_sens.loc[df_IS_sens.index[3],df_IS_sens.columns[:]].mean() 
                    df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[3], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[3], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[3],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[3], df_IS_sens.columns[:]].mean() )**2
                                )
            else:       #Medium resolution,  MR
                if int(mass) in range(59,81): #mass from 59 to 80 included both
                    df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]
                        ] / df_IS_sens.loc[df_IS_sens.index[4],df_IS_sens.columns[:] ] * df_IS_sens.loc[
                            df_IS_sens.index[4],df_IS_sens.columns[:]].mean() 
                    df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[4], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[4], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[4],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[4], df_IS_sens.columns[:]].mean() )**2
                                )   
                else:       #rest of mass ranges, only 2 IS here xD
                    df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]
                        ] / df_IS_sens.loc[df_IS_sens.index[5],df_IS_sens.columns[:] ] * df_IS_sens.loc[
                            df_IS_sens.index[5],df_IS_sens.columns[:]].mean() 
                    df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[5], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[5], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[5],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[5], df_IS_sens.columns[:]].mean() )**2
                                )               
    elif  IS_meas == ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 
               'Co59(MR)', 'In115(MR)', 'Ho165(MR)', 'Th232(MR)']:   #######case 2, all IS measured in LR and MR! 

        for i in range(np.where(df_raw.index == 'IS conc [ppb]')[0][0]): #loop through all isotopes
                #df_raw.loc[229,'Isotopes'] = df_raw.iloc[226,0]#relation iloc, loc
        
            mass, resol = Get_A_Resol(df_raw.index[i])
        #No we neeed to see the resolution, since for LR we have 4, and for MR only 2
            if resol == 'LR':           #low resolution
            #Now, which mass?
                if int(mass) in range(59,81): #mass from 59 to 80 included both
                    df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]
                        ] / df_IS_sens.loc[df_IS_sens.index[0],df_IS_sens.columns[:] ] * df_IS_sens.loc[
                            df_IS_sens.index[0],df_IS_sens.columns[:]].mean() 
                    df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[0], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[0], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[0],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[0], df_IS_sens.columns[:]].mean() )**2
                                )
                elif int(mass) in range (81, 139):
                    df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]
                        ] / df_IS_sens.loc[df_IS_sens.index[1],df_IS_sens.columns[:] ] * df_IS_sens.loc[
                            df_IS_sens.index[1],df_IS_sens.columns[:]].mean() 
                    df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[1], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[1], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[1],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[1], df_IS_sens.columns[:]].mean() )**2
                                )                                                     
                elif int(mass) in range (139, 210):
                    df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]
                        ] / df_IS_sens.loc[df_IS_sens.index[2],df_IS_sens.columns[:] ] * df_IS_sens.loc[
                            df_IS_sens.index[2],df_IS_sens.columns[:]].mean() 
                    df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[2], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[2], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[2],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[2], df_IS_sens.columns[:]].mean() )**2
                                )                
                else:   #if mass in range(210, 249):    rest of mass range
                    df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]
                        ] / df_IS_sens.loc[df_IS_sens.index[3],df_IS_sens.columns[:] ] * df_IS_sens.loc[
                            df_IS_sens.index[3],df_IS_sens.columns[:]].mean() 
                    df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[3], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[3], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[3],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[3], df_IS_sens.columns[:]].mean() )**2
                                ) 
            else:       #Medium resolution,  MR
                if int(mass) in range(59,81): #mass from 59 to 80 included both
                    df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]
                        ] / df_IS_sens.loc[df_IS_sens.index[4],df_IS_sens.columns[:] ] * df_IS_sens.loc[
                            df_IS_sens.index[4],df_IS_sens.columns[:]].mean() 
                    df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[4], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[4], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[4],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[4], df_IS_sens.columns[:]].mean() )**2
                                )             
                elif int(mass) in range (81, 139):
                    df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]
                        ] / df_IS_sens.loc[df_IS_sens.index[5],df_IS_sens.columns[:] ] * df_IS_sens.loc[
                            df_IS_sens.index[5],df_IS_sens.columns[:]].mean() 
                    df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[5], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[5], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[5],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[5], df_IS_sens.columns[:]].mean() )**2
                                )            
                elif int(mass) in range (139, 210):
                    df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]
                        ] / df_IS_sens.loc[df_IS_sens.index[6],df_IS_sens.columns[:] ] * df_IS_sens.loc[
                            df_IS_sens.index[6],df_IS_sens.columns[:]].mean() 
                    df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[6], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[6], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[6],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[6], df_IS_sens.columns[:]].mean() )**2
                                )            
                else:   #if mass in range(210, 249):    rest of mass range
                    df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]
                        ] / df_IS_sens.loc[df_IS_sens.index[7],df_IS_sens.columns[:] ] * df_IS_sens.loc[
                            df_IS_sens.index[7],df_IS_sens.columns[:]].mean() 
                    df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[7], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[7], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[7],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[7], df_IS_sens.columns[:]].mean() )**2
                                )     
    elif  IS_meas == ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)']:   #######case 3, only LR IS!
        for i in range(np.where(df_raw.index == 'IS conc [ppb]')[0][0]): #loop through all isotopes
                #df_raw.loc[229,'Isotopes'] = df_raw.iloc[226,0]#relation iloc, loc
        
            mass, resol = Get_A_Resol(df_raw.index[i])
            #Now, which mass?
            if int(mass) in range(59,81): #mass from 59 to 80 included both
                df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]
                    ] / df_IS_sens.loc[df_IS_sens.index[0],df_IS_sens.columns[:] ] * df_IS_sens.loc[
                        df_IS_sens.index[0],df_IS_sens.columns[:]].mean() 
                df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[0], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[0], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[0],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[0], df_IS_sens.columns[:]].mean() )**2
                                )             
            elif int(mass) in range (81, 139):
                df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]
                    ] / df_IS_sens.loc[df_IS_sens.index[1],df_IS_sens.columns[:] ] * df_IS_sens.loc[
                        df_IS_sens.index[1],df_IS_sens.columns[:]].mean() 
                df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[1], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[1], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[1],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[1], df_IS_sens.columns[:]].mean() )**2
                                )                           
            elif int(mass) in range (139, 210):
                df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]
                    ] / df_IS_sens.loc[df_IS_sens.index[2],df_IS_sens.columns[:] ] * df_IS_sens.loc[
                        df_IS_sens.index[2],df_IS_sens.columns[:]].mean() 
                df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[2], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[2], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[2],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[2], df_IS_sens.columns[:]].mean() )**2
                                )                   
            else:   #if mass in range(210, 249):    rest of mass range
                df_IS_co.loc[df_raw.index[i],df_IS_co.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]
                    ] / df_IS_sens.loc[df_IS_sens.index[3],df_IS_sens.columns[:] ] * df_IS_sens.loc[
                        df_IS_sens.index[3],df_IS_sens.columns[:]].mean() 
                df_IS_co_std.loc[df_raw.index[i],df_IS_co_std.columns[:]] = df_IS_co.loc[df_raw.index[i],
                            df_IS_co.columns[:]] * np.sqrt(
                                (df_raw_std.loc[df_raw_std.index[i], df_raw_std.columns[:]] / df_raw
                                 .loc[df_raw.index[i], df_raw.columns[:]] )**2 + (df_IS_sens_std
                                .loc[df_IS_sens_std.index[3], df_IS_sens_std.columns[:] ] / 
                                df_IS_sens.loc[df_IS_sens.index[3], df_IS_sens.columns[:] ]
                                )**2 + ( df_IS_sens.loc[df_IS_sens.index[3],
                                        df_IS_sens.columns[:]].std() / df_IS_sens
                                        .loc[df_IS_sens.index[3], df_IS_sens.columns[:]].mean() )**2
                                )      
    else:           #no contempled case
        print('\n---------------------------\n')
        print('The Internal standards used are not defined as a case, do it! (if statements)\n')
        print("Das Lebe ist kein Ponyhof / Embeses la bida no e kmo keremo / C'st la vie, mon ami\n")
    
    ################## Return ###############
    '''
    After the loop ends, we can return it
    '''
    return df_IS_co, df_IS_co_std



#%%########### 1.9) ICPMS Blank correction #############
#####################################
def ICPMS_ICPMSBlanks_corrector(df_IS_co, df_IS_co_std, columns_blks):
    
    '''
    PART OF THE ICPMS Data processing!
    
    Function to be run after the IS sensitivity correction, to apply the next step in the ICPMS data
    analysis process, the ICPMS blanks correction (std 0ppt and blanks). 
    
    Note that in the excel, int he part of the ppb table, you should delete the names, that stefaan write twice,
    for this to work!!!! 
    Also needed that the ppb data table is introduced by:  "IS conc [ppb]". beware, sometimes Stef writes (ppb), not [ppb]!!
    
    
    *Inputs:
        .columns_blks: np.array([]) indicating the number of the columns contaning blanks (std 0ppt and blank std).
            Ex: columns_blanks = np.array([1, 2]). Those columns number are from excel!
        .df_IS_co: df containing the cps (also ppb, not needed but there it is), output from the IS sens correction
            funciton. The formatting is like the excel. 
        .df_IS_co_std: df containing the std of the cps corrected from the IS
    *Outputs:
        .df containinng the cps data corrected for the ICPMS blanks. Note it also contains the ppb data, 
            but now modified so they are random numbers. 
        .df containing the std of the cps corrected for the ICPMS blanks. Quadratic error propagation used!
            
            
    ##### To DO ###########
        *Improve style and remove ppb data (would need modifications for the IS sens function)
    
    '''
    
    ###################################################################
    '''
    So, for the Blank correction, I need:
        1) To put ina  df all of the blanks. I could indicate column numbers, KISS
        2) Compute average from the cps of the blank
        3) Substract either the average or the ICPMS blanks, or the minimum cps value of the samples (from IS corrected, no 
                 ICPMS blanks), to ensure no negative values. 

    I would really need the IS corrected data without the sens also bro, but I could get it easily I think
    '''
    
    columns_blks = columns_blks - 1 #to adapt to the system in the df, isotopes are index, not columns!
    
    #1) is trivial:
    df_blanks = df_IS_co.iloc[:, columns_blks ]          #df with the ICPMS blanks
                    #no isotope label contained!!
    df_blanks_std = df_IS_co_std.iloc[:, columns_blks ]
    
    ######### Debug: printing the labels of df ICPMS blanks  ##########
    
    print('\n ######################################')
    print('\n The columns used as ICPMS blanks are (should be std 0ppt and blank, ' + 
          'if not, correct the columns that you give as an input!): \n')
    print(df_blanks.columns)            #printing columns
    print('\n ######################################')
    
    #2) could be done with df
    df_blanks['<cps>'] = df_blanks.mean(axis = 1)            #I append the mean values. Gives a warning, 
                                                            #but nothing is wrong!
    df_blanks_std ['<cps>'] = df_blanks.std(axis = 1) 

    #3) The substraction
    '''
    THe first step is getting the minimun values of the cps data
    Those values are for the columns of the samples except the blanks!
    '''
    df_Is_co_no_blk = df_IS_co.drop(df_IS_co.columns[columns_blks ], axis = 1)
                        #df initial without the ICPMS blanks!! 

    #I could do the min to that, but since also contain the cps data there are some NaN, which make things not 
    #work. If I do fill nan with the mean values, that could fix it. Lets try bro! 

    df_Is_co_no_blk.fillna(999999.999, inplace = True)       #filling NaN values with 99999 (easily recognizable)

    #Now the min values are, storing them in the same df:
    df_Is_co_no_blk['min cps'] = df_Is_co_no_blk.min(axis = 1) 
                #axis = 1 for columns!
    '''
    For this I would need a loop to check for each row which value to substract. I would say we always substract
    the mean, and Stefaan also think that. adn I probe that matematically, so it is like that xD

    So, now the problem is perform that operation here on python, since the df are different.

    I needed to delete some column nsames in teh ppb data!!!!!!!!!!
    
    The loop is in the same fashion as the one for IS sens correction :D
    '''

    df_IS_blk_co = df_IS_co.copy(  )        #initilization df for the loop
    df_IS_blk_co_std = df_IS_co.copy(  )

    for i in range(np.where(df_IS_co.index == 'IS conc [ppb]')[0][0]): #loop through all isotopes
                    #[0]s we get the desired value!
        if df_Is_co_no_blk['min cps'][i] <= df_blanks['<cps>'][i]:     #min low, so substracting min of
                        #cps data (no blanks)
            
            df_IS_blk_co.loc[df_IS_co.index[i],df_IS_blk_co.columns[:]]= df_IS_co.loc[df_IS_co.index[i],
                                df_IS_co.columns[:]] - df_Is_co_no_blk['min cps'][i]
                            #substraction of the min value        
            df_IS_blk_co_std.loc[df_IS_co_std.index[i],df_IS_blk_co_std.columns[:]] = np.sqrt(2) * df_IS_co_std.loc[
                df_IS_co_std.index[i], df_IS_co_std.columns[:]]                #error
        else:           #mean lower, so substracting the mean of ICPMS blanks
            df_IS_blk_co.loc[df_IS_co.index[i],df_IS_blk_co.columns[:]]= df_IS_co.loc[df_IS_co.index[i],
                                df_IS_co.columns[:]] - df_blanks['<cps>'][i]
            df_IS_blk_co_std.loc[df_IS_co_std.index[i],df_IS_blk_co_std.columns[:]] = np.sqrt(
                (df_IS_co_std.loc[df_IS_co_std.index[i], df_IS_co_std.columns[:]]
                 )**2 + (df_blanks_std['<cps>'][i])**2
                )
    
    #FInally we delete the blanks before returning the data: 
        
    df_IS_blk_co_final = df_IS_blk_co.drop(df_IS_blk_co.columns[columns_blks], axis = 1)
                #for some reason I can not do the drop and inplace = True, dont work
    df_IS_blk_co_std_final = df_IS_blk_co_std.drop(df_IS_blk_co_std.columns[columns_blks], axis = 1)


    ########### Return ###########
    return df_IS_blk_co_final, df_IS_blk_co_std_final         #return of the data



#%%######## 1.10) ICPMS data processing automatized #############
#####################################

def ICPMS_data_process(df_cps, df_rsd, ICPblk_columns, 
                       name_plot_LR_bef = 'IS_sensLR_plotBEF', name_plot_MR_bef = 'IS_sensMR_plotBEF',
                       name_plot_LR_aft = 'IS_sensLR_plot', name_plot_MR_aft = 'IS_sensMR_plot',
                       IS_meas = ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)'],
                       excel_name = 'df_IS_Blks_corr.xlsx'):
    '''
    SUITE of ICPMS Data Processing!
    
    Function that will apply all the steps for the automatization of the ICPMS data processing:
        1) Read cps amd ppb data and %rsd 
        2) Compute sens = cps/ppb and plot it
        3) Apply the sensitivity correction
        4) plot the new sensitivtiy plots (should be straight lines)
        5) Substract ICPMS blanks
        6) Save that data (to calibrate after, by hand, for the moment..)
    This function computes as well the std using quadratic uncertainty propagation. Note the excel sheets containing
    the std also contain info on the ppb, but tha tinfo is just nonsense, I didnt erase it not to have problem with
    dimensions!
    
    Important notes:
        . To do 2), its needed that the ppb data table begins with IS conc [ppb]! Befpre there must be only
                                                the cps data, nothing else, no text nor anything!!
        . To do 3), you define the IS cases (which IS are to be used), so beware, 
                maybe your case its not (yet) defined!!  
        
    *Inputs:
        .df_cps: df containing the cps data and also the ppb data, in the classic format. Must not contain the wash, 
            will give errors (divide by zero). Take care of the names (like for the cps table), they are crutial 
            for the sens calc and  correction! Also about the format, not rows with 0s etc. 
                    Take a lot of care!!!
        
        .IS_meas: array containing in a list the IS and its resolution, like how they appear in the isotopes column. 
            Default value:
                ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)']
            That mean those isotopes were measured. If Ho165(MR) also measured, just included it, and fine ;)
        .ICPblk_columns: np array containing the columns numbers where the ICPMS blanks are (blank and std 0ppt). 
            Numbers from the excel (A = 0, B = 1, etc)
        .name_plot_L(M)R_bef(aft): name of the plots of the IS sensitivity before and after doing the IS correction. 
        .excel_name: string that will be the name of the file containing the output. Default: 'df_IS_Blks_corr.xlsx'
        
    *Output:
        .df containing the data applying the ICPMS blanks and sensitivity corrections :)
    
    
    To Do:
            .Automatize more stuff? such as the sens corrections, not by case something better, more general??
            .Optimize excel saving!
            :Delete the ppb data when not needed, so in excels it doesnt exist neither?
    '''
    
    ##### 0 ) Std calc ###########
    'From the %rsd and the cps the std is trivial, %rsd = std/cps * 100 ==> std = cps*%rstd/100'
    
    df_std = ICPMS_std_calculator(df_cps, df_rsd)           #Std of the cps
    
    ###### 1) IS sens calc ######
    df_IS_sens, df_IS_sens_std = IS_sens_calculator_plotter(df_cps, df_std, IS_meas, 
                    name_IS_sens_LR_plot = name_plot_LR_bef, name_IS_sens_MR_plot = name_plot_MR_bef)         
                #calculation the IS correction
                #I define names of the plot, so the other ones have the default name
    
    print('\n ########################')
    print('Step 1. done, computed IS sens :) ')
    print('#############################')
        #printing that tha step went good, would be good for debuggind!
        
        
    ###### 2) IS sens correction and new sens calc ##########
    df_IS_co, df_IS_co_std = IS_sens_correction(df_cps, df_std, df_IS_sens, df_IS_sens_std, IS_meas)
                               #applying the IS correction to the cps data
    
    print('\n ########################')
    print('Step 2.1 done, applyed the IS correction :)) ')
    print('#############################')    
    
    df_IS_sens_co, df_IS_sens_co_std = IS_sens_calculator_plotter(df_IS_co, df_IS_co_std, IS_meas, 
                name_IS_sens_LR_plot= name_plot_LR_aft, name_IS_sens_MR_plot= name_plot_MR_aft)    
                                #getting and plotting new IS sens
    
    print('\n ########################')
    print('Step 2.2. donce, copmuted the new IS sens :)) ')
    print('#############################')
    
    
    ##### 3)ICPMS Blk correction #########
    df_Blks_co, df_Blks_co_std = ICPMS_ICPMSBlanks_corrector(df_IS_co, df_IS_co_std, ICPblk_columns)    
                                        #correcting for ICPMS blanks
    print('\n ########################')
    print('Step 3 (final). done, applyed the Blk correction :))) ')
    print('#############################')

    ##### 4) Saving and Output #########
    '''
    Here I want to save the df after IS correction, and after the Blk correction. Both steps would be nice, 
    for debugging!
    
    '''
    writer = pd.ExcelWriter(excel_name, engine = 'xlsxwriter')      #excel writer

    #########Blk correction data
    df_Blks_co.to_excel(writer, sheet_name = 'Blk_correction', startrow = 6, header = False, freeze_panes = (6, 1))            
                #saving to excel. I freeze row 6 and column 1, so I can see all the data in a good way :)
                        #Chatgpt helped me to get the format I want, the one that Stefaan uses :)

    excel_sheet = writer.sheets['Blk_correction']
    bold_format = writer.book.add_format({'bold': True})      
    excel_sheet.write_row('B2', df_Blks_co.columns, bold_format)      #Write from cell B2 with the numer of columns
                    #Note B2 is column B, row 2
    excel_sheet.write_row('B1', range(1, len(df_Blks_co.columns) + 1), bold_format)  #2nd row with columns names
    excel_sheet.write_row('B3', [None] * len(df_Blks_co.columns))            #row 3 empty
    excel_sheet.write_row('B4', ['Net <Int>'] * len(df_Blks_co.columns))         #row 4 with a value repeated
    excel_sheet.write_row('B5', ['[cps]'] * len(df_Blks_co.columns))     #row 5 with a value repeated

    #########std Blk correction data
    df_Blks_co_std.to_excel(writer, sheet_name = 'Blk_correction_std', startrow = 6, header = False, freeze_panes = (6, 1))            

    excel_sheet = writer.sheets['Blk_correction_std']
    bold_format = writer.book.add_format({'bold': True})      
    excel_sheet.write_row('B2', df_Blks_co_std.columns, bold_format)
    excel_sheet.write_row('B1', range(1, len(df_Blks_co_std.columns) + 1), bold_format)  #2nd row with columns names
    excel_sheet.write_row('B3', [None] * len(df_Blks_co_std.columns))            #row 3 empty
    excel_sheet.write_row('B4', ['Net Int std'] * len(df_Blks_co_std.columns))         #row 4 with a value repeated
    excel_sheet.write_row('B5', ['[cps]'] * len(df_Blks_co_std.columns))     #row 5 with a value repeated
    
    
    ##########Is correction data
    df_IS_co.to_excel(writer, sheet_name = 'IS_correction', startrow = 6, freeze_panes = (6, 1),
                             header = False)        #saving to excel in another sheet
            #THe start trow make that Co59 is on row 7, as it should be!
    df_IS_sens_co.to_excel(writer, sheet_name = 'IS_correction', startrow = 6 + df_IS_co.shape[0] + 2,
                           header = False)
                        #putting the new IS sensitivity below the IS corrected data!!
    
    excel_sheet2 = writer.sheets['IS_correction']
    excel_sheet2.write_row('B2', df_IS_co.columns, bold_format)      #1st row with numbers of columns
    excel_sheet2.write_row('B1', range(1, len(df_IS_co.columns) + 1), bold_format)  #2nd row with columns names
    excel_sheet2.write_row('B3', [None] * len(df_IS_co.columns))            #row 3 empty
    excel_sheet2.write_row('B4', ['<Int>'] * len(df_IS_co.columns))         #row 4 with a value repeated
    excel_sheet2.write_row('B5', ['[cps]'] * len(df_IS_co.columns))     #row 5 with a value repeated
    
    ############ IS correction std data
    df_IS_co_std.to_excel(writer, sheet_name = 'IS_correction_std', startrow = 6, freeze_panes = (6, 1),
                             header = False)        #saving to excel in another sheet
    df_IS_sens_co_std.to_excel(writer, sheet_name = 'IS_correction_std', startrow = 6 + df_IS_co_std.shape[0] + 2,
                           header = False)
                        #putting the new IS sensitivity below the IS corrected data!!
    
    excel_sheet2 = writer.sheets['IS_correction_std']
    excel_sheet2.write_row('B2', df_IS_co_std.columns, bold_format)      #1st row with numbers of columns
    excel_sheet2.write_row('B1', range(1, len(df_IS_co_std.columns) + 1), bold_format)  #2nd row with columns names
    excel_sheet2.write_row('B3', [None] * len(df_IS_co_std.columns))            #row 3 empty
    excel_sheet2.write_row('B4', ['std'] * len(df_IS_co_std.columns))         #row 4 with a value repeated
    excel_sheet2.write_row('B5', ['[cps]'] * len(df_IS_co_std.columns))     #row 5 with a value repeated    
    
    
    #writer.save()                                           #critical step, save the excel xD 
    writer.close()      #save was deprecated bro xD

    ################ Return ####################
    '''
    I return the df with the Blk corrected data, although I only need the excel, which is what I
    will use for the calib
    '''
    
    return df_Blks_co



#%%######## 1.11) ICPMS Isotope selector #############
#####################################
def ICPMS_Isotope_selector(df_cps, Isotopes):
    '''
    Function to choose from all the elements measured in the DF some selected elements. This funciton will be used mainly
    to plot all those elements in a single plot, since for that having one df witht he info is the best.
    
    
    *Inputs:
        .df_cps: df with the cps/ppb from the ICPMS measurement. Same format as "always" (jsut impleemnted xD), index are the isotopes
        then samples, 1streplicate, then 2nd
        .Isotopes: np.array containing the isotopes that I want, including the reoslution: np.array[('Sr84(LR)', 'Si28(MR)')]
        
        columns_blks: np.array([]) indicating the number of the columns contaning blanks (std 0ppt and blank std).
            Ex: columns_blanks = np.array([1, 2])
        .df_IS_co: df containing the cps (also ppb, not needed but there it is), output from the IS sens correction funciton.
            The formatting is like the excel. 
    *Outputs:
        .df like the input one, but containing only the results for the desired isotopes
            

    '''
    #
    '''
    TO get the columns, given the relevant elements, I could do what I do
    in the cps/ppb function:
        I have a np aray with the elemnts to find (IS_meas), and I find them like that:
            
    value_to_find = IS_meas[i]
    matching_rows = df_cps_ppb_dat.loc[df_cps_ppb_dat.iloc[:,0] == value_to_find]  #give the full row!
    '''
    
    df_elem_rel = pd.DataFrame()         #Empty df to store the values

    for i in range(0, len(Isotopes)):     #Loop through all the relevant elements

       matching_rows = df_cps.loc[df_cps.index == Isotopes[i]]  #give the full row!
                   #rows whose elements are the relevant elements
           
           #TO store temporarily those values I create an auxiliary df
       df_aux = pd.DataFrame(data = matching_rows)

           #And I add that to the storing df
       df_elem_rel = df_elem_rel.append(df_aux, ignore_index= False)
           #ignore index false to store the index, key!

    ########### Return ###########
    return df_elem_rel        #return of the data





#%%############ 1.12) Kd calculaor #############
#####################################

def ICPMS_KdQe_calc (df_data, df_VoM_disol, df_m_be, Nrepl = 2, 
                     ret_Co__Ceq = False):
    '''
    Function that will compute the distribution constant Kd and the adsorption 
    quantity Q_e from the ppb data obtained with ICPMS. Note that data must be 
    corrected for the dilutions factors. 
    
    Based on the blk corrector function.
    
    THe distribution constant Kd is:
        K_d = (C_0 - C_eq)/C_eq * V/m = Q_e / C_eq;
    being            
        C_0 = initial concentration [M]
        C_eq = concentration at equilibrium [M]
        m = mass of dry bentonite [g]
        V = volume of the solution [L]
        Q_e = absorbed quantity in equil [g soluto / g bent]
        
    In our case, that we measure C in ppb = ng soluto /g disol, we need to mutiply
    by g tot / m bent, to achieve the same units in Q_e!! And we do not have
    equilibirum since its kinetic, so Qe, Ke ==> Q(t), K(t)

    Necessary that the data contain no Div0, ensure in the excel by using the 
    iferror(operation, 0) function!
    
    Note this requires a df series with the volume, that you were not measuring
    in the first exp
    (up to 8/23). Note ten that units are involved!. If measuring mass ing and 
    volumes in L, Q


    *Inputs:
        .df_data: dataframe containing the ICPMS data, the full data, with the 2 or 3 
            replicates. Should be Dfs corrected. Format: isotopes as index, 
            columns the samples, 1st 1st replicate, then 2nd replicate. Note 
            the 1st sample in each replicate must be the sample blank! This
            samples will be removed from the Qe and Kd df!
        .df_VoM_disol: pd series containing the volume [mL] added to the bottle
            of the solution, BIC, or whatever. normally 50ml or the total 
            mass of the solution [g]. If df_data in ppb, this must be the total
            mass so that Q_e is in g/g !
        .df_m_bent: pd series contaning the mass of bentonite [g] in the bottle
        (normally 250mg)
        .Nrepl: number of replicates. Default value = 2. 3 also accepted
        ret_Co__Ceq: if True, returns a df with C_0 - C_eq = False
    
    *Outputs (in that order):
        .df with the Kd data
        .df with q_e data
        
        '''
    
    
    ########## 0) Precalcs ##########
    '''
    To avoid the div0 error, which occurs when I have 0 as a values in the df, 
    which I have for all the elements that were not found by ICPMS, I can just 
    put NaN instead, since that will not give the Div0 error when computing Kd
    
    13/2/25, To do this, I was doing:
        df_data.replace(0, np.nan, inplace=True)  #replace 0 values with NaN, 
                                                  #to avoid Div0 error!
        
    But that introduced a problem, NaN in the ppb values, since I do inplace =
    True. I want to do that mod only for this calcs, so I should create a new
    variable.
    '''
    
    df_data_aux = df_data.copy()
    df_data_aux.replace(0, np.nan, inplace=True)  #replace 0 values with NaN, 
                                              #to avoid Div0 error!
    
    
    ########### 1) Calcs ###########
    '''
    The operations to perform are:
        1) C_0 - C_eq (>0)
        2) 1) * V/m = q_e
        3) 2)	1/C_eq = Kd
    
    I must treat the 2 experiments are different, I should substract the blank 1 
    to the 1st emasurements and the 2 to the others. Since I ordered it in the 
    right way (1st replicacte 1, then replicate 2, I could) split it easily :D
            df.shape gives the shape of the df, n_rows, n_columns
    
    Note the df have number of samples * 2 replicates columns.
    
    Then, I will create a new dataframe substracting that data. To do so, I 
    need to get rid of the isotopes column, since is text, and then add it again.
    Watch, the substraction is easy with a pandas mehotd.

    I shuold then remove those columns
    from there, and replace negatives values for 0, for a good plot
    '''
    
    #So, lets split into the 2 replicates!
    '''
    For 2 replicates its easy, for 3 it could be more tricky. Beware! TO create
    a function you should say the number of replicates and so!
    '''
    
    if Nrepl == 2:          #Standard case, 2 replicates
        df_1 = df_data_aux.iloc[ :, 0: round( ( df_data_aux.shape[1] ) / 2 ) ] #1st replicate
        df_2 = df_data_aux.iloc[ :, round( ( df_data_aux.shape[1] ) / 2 ) :  ] #replicate 2
    
        df_VoM_1 = df_VoM_disol.iloc[ 0: round( ( df_VoM_disol.shape[0] ) / 2 ) ]      #1st replicate
        df_VoM_2 = df_VoM_disol.iloc[ round( ( df_VoM_disol.shape[0] ) / 2 ) :  ]       #replicate 2
            #Achtung! In shape I put 0, because they are series, so 1D!!!!
            #VoM = Volume or Mass!
        df_m_1 = df_m_be.iloc[ 0: round( ( df_m_be.shape[0] ) / 2 ) ] #1st replicat
        df_m_2 = df_m_be.iloc[ round( ( df_m_be.shape[0] ) / 2 ) :  ] #replicate 2    
    
    #######Future note: here you see the automatization to N-replicates, doing 
        #this with a function. Then the operations you can done them, grouping 
        #the df in an array, and for element in array, perform them!
    
    ###### 1) C_0 - C_eq = - (C_eq - C0)
    #I will do C_eq - C0, and then invert that, since its easier. C0 is the blank data, 
    #thats why is easier, so I can copy paste the blank substraction
    
        dfCeq__C0_1 = df_1.subtract(df_1.iloc[:,0], axis = 0 ) #doing the substraction
        dfCeq__C0_1.drop( [df_1.iloc[:,0].name], axis = 1, inplace = True) #drop blank column
        #
        dfCeq__C0_2 = df_2.subtract(df_2.iloc[:,0], axis = 0 ) #Replicate 2
        dfCeq__C0_2.drop( [df_2.iloc[:,0].name], axis = 1, inplace = True)
    
        #Now lets invert the sign:
        dfC0__Ceq_1 = - dfCeq__C0_1
        dfC0__Ceq_2 = - dfCeq__C0_2

    ######## 2) Apply the V/ m giving q_e (from Df_exp)
    #For this I ned to remove the blank columns to both m and V, since from 
    #C0-Ceq they are removed!

        df_m_1 = df_m_1[1:]  #fast way to delete 1st elemen (blank) in a series
                #new_series = data.drop(data.index[0]) also work, from Chatgpt
        df_m_2 = df_m_2[1:]
        df_VoM_1 = df_VoM_1[1:]
        df_VoM_2 = df_VoM_2[1:]
    
    #And now I can operate:
        
        df_Qe_1 = dfC0__Ceq_1 * df_VoM_1 / df_m_1
        df_Qe_2 = dfC0__Ceq_2 * df_VoM_2 / df_m_2 
    
    ######## 3) Apply 1/C_eq = Kd
        df_Kd_1 = df_Qe_1 / df_1.drop( [df_1.iloc[:,0].name], axis = 1)   
         #Not df_1 contains blk (1st column), so I remove it for the operation!    
         #This also works: df_Qe_1.div(df_1.drop( [df_1.iloc[:,0].name], axis = 1) )
        df_Kd_2 = df_Qe_2 / df_2.drop( [df_2.iloc[:,0].name], axis = 1)
    
    #Now lets add them together
        
        df_C0__Ceq = pd.concat( [dfC0__Ceq_1, dfC0__Ceq_2], axis = 1) 
                        #merging the 2 df (replicates) in a hugeone 
        df_Kd = pd.concat( [df_Kd_1, df_Kd_2], axis = 1)         
        df_Qe = pd.concat( [df_Qe_1, df_Qe_2 ] , axis = 1)

    elif Nrepl == 3:            #3 replicates
    #Gathering the replicates sepparately    
        df_1 = df_data_aux.iloc[:, : round(df_data_aux.shape[1] / 3)]
        df_2 = df_data_aux.iloc[:, round(df_data_aux.shape[1] / 3): 2*round(df_data_aux.shape[1] / 3)]
        df_3 = df_data_aux.iloc[:, 2*round(df_data_aux.shape[1] / 3) :]
        
        df_VoM_1 = df_VoM_disol.iloc[ 0: round( ( df_VoM_disol.shape[0] ) / 3 ) ]      
                                                #1st replicate
        df_VoM_2 = df_VoM_disol.iloc[ round( ( df_VoM_disol.shape[0] ) / 3 ): 2* round( ( df_VoM_disol.shape[0] ) / 3 ) ]      
        df_VoM_3 = df_VoM_disol.iloc[ 2 *round( ( df_VoM_disol.shape[0] ) / 3 ) : ]      
                                            #3rd replicate
        
        df_m_1 = df_m_be.iloc[ : round( ( df_m_be.shape[0] ) / 3 ) ] 
        df_m_2 = df_m_be.iloc[ round( ( df_m_be.shape[0] ) / 3 ) : 2* round( ( df_m_be.shape[0] ) / 3 )]
        df_m_3 = df_m_be.iloc[ 2 *round( ( df_m_be.shape[0] ) / 3 ) : ]  
    
        #1) 
        dfCeq__C0_1 = df_1.subtract(df_1.iloc[:,0], axis = 0 )#doing the substraction
        dfCeq__C0_2 = df_2.subtract(df_2.iloc[:,0], axis = 0 )           
        dfCeq__C0_3 = df_3.subtract(df_3.iloc[:,0], axis = 0 )       
        
        dfCeq__C0_1.drop( [df_1.iloc[:,0].name], axis = 1, inplace = True)   
                                    #drop blank column
        dfCeq__C0_2.drop( [df_2.iloc[:,0].name], axis = 1, inplace = True)
        dfCeq__C0_3.drop( [df_3.iloc[:,0].name], axis = 1, inplace = True)     
    
        dfC0__Ceq_1 = - dfCeq__C0_1
        dfC0__Ceq_2 = - dfCeq__C0_2
        dfC0__Ceq_3 = - dfCeq__C0_3
    
    ######## 2) Apply the V/ m giving q_e (from Df_exp)
    #For this I ned to remove the blank columns to both m and V, since from 
    #C0-Ceq they are removed!

        df_m_1 = df_m_1[1:]  #fast way to delete 1st elemen (blank) in a series
        df_m_2 = df_m_2[1:]
        df_m_3 = df_m_3[1:]
        df_VoM_1 = df_VoM_1[1:]
        df_VoM_2 = df_VoM_2[1:]
        df_VoM_3 = df_VoM_3[1:]
        
    #And now I can operate:
        
        df_Qe_1 = dfC0__Ceq_1 * df_VoM_1 / df_m_1
        df_Qe_2 = dfC0__Ceq_2 * df_VoM_2 / df_m_2 
        df_Qe_3 = dfC0__Ceq_3 * df_VoM_3 / df_m_3
        
        
    ######## 3) Apply 1/C_eq = Kd
        df_Kd_1 = df_Qe_1 / df_1.drop( [df_1.iloc[:,0].name], axis = 1)   
        df_Kd_2 = df_Qe_2 / df_2.drop( [df_2.iloc[:,0].name], axis = 1)
        df_Kd_3 = df_Qe_3 / df_3.drop( [df_3.iloc[:,0].name], axis = 1)
    
    #Now lets add them together
        df_Kd = pd.concat( [df_Kd_1, df_Kd_2, df_Kd_3], axis = 1)         
        df_Qe = pd.concat( [df_Qe_1, df_Qe_2, df_Qe_3 ] , axis = 1)  
        df_C0__Ceq = pd.concat( [dfC0__Ceq_1, dfC0__Ceq_2, dfC0__Ceq_3], axis = 1)
    
    else:               #Error case
        print('Error case, wrong Nrepl introduced!')    
        df_Kd = 0
        df_Qe =0
                
      
   ####Checked that the Qe calc is correct, and then Kd must be also :)     
      
    ########### 2) Return #############
    #Here the if for returning or not C_0 - C(t) applies
    #Before returning them, lets convert them into numeric:
        
    df_Kd = df_Kd.apply(pd.to_numeric)       #converting it to numeric 
    df_Qe = df_Qe.apply(pd.to_numeric) 
    df_C0__Ceq =df_C0__Ceq.apply(pd.to_numeric)        
    
    if ret_Co__Ceq == False:            #do not return it
        return df_Kd, df_Qe             #return
    
    else:                   #return C0-Ceq
        return df_Kd, df_Qe, df_C0__Ceq


#%% ########## 1.13) Kd calculaor, Adsorption version #############
#####################################
def ICPMS_KdQe_calc_Ad (df_mother_sol, df_samples, df_VoM_disol, df_m_be, 
                        df_m_liq = False, ret_Co__Ceq = False, 
                        removed_sample_1 = True):
    '''
    Function that will compute the distribution constant Kd and the adsorption 
    quantity Q_e from the ppb data obtained with ICPMS. Note that data must be
    corrected for the dilutions factors. 
    
    Based on the blk corrector function. The sorbed quantity Q_e is:
        Q_e = (C_0 - C_e) * V/m
    
    THe distribution constant Kd is:
        K_d = (C_0 - C_eq)/C_eq * V/m = Q_e / C_eq;
    being            
        C_0 = initial concentration [M]
        C_eq = concentration at equilibrium [M]
        m = mass of dry bentonite [g]
        V = volume of the solution [L]
        Q_e = absorbed quantity in equil [g soluto / g bent]
        
    In our case, that we measure C in ppb = ng soluto /g disol, we need to 
    mutiply by g tot / m bent, to achieve the same units in Q_e!! Then our 
    equations are:
        Q_e = (C_0 - C_e) * m_tot /m
        K_d = Q_e / C_e
        
    Improvement made thanks to Espe (3/2024). Note that after the solid liquid 
    sepparation,
    you retrieve a bit less of liquid, since the bentonite adsorbs some. This
    means that the concentration you measure with the ICP-MS is the concentration
    of this liquid, with mass m_liq_e <= m_liq_0. Then the Q_e should be 
    calculated like this, using that final mass!.
    
    
    Note that in the case of no equilibrium (eg Kinetics exp), Q_e, Kd ==> 
    Q_e(t), K(t). 
    
    Here we have different moher solutions, whihc is the main different
    from the other function, where all the samples had a common mother solution 
    (C_0). Here we have several C_0!

    Necessary that the data contain no Div0, ensure in the excel by using the 
    iferror(operation, 0) function!
    
    Note this requires a df series with the volume, that you were not measuring
    in the first exp (up to 8/23). Note ten that units are involved!. If 
    measuring mass ing and volumes in L, Q.
    
    Note as well that the 1st sample/replicate is the procedural blank, 
    not used here!!
    beware if first sample is not procedural blank!!!!!!!


    *Inputs:
        .df_samples: dataframe containing the data, the full data, with the 3 
        replicates. Should be Dfs corrected
            Format: isotopes as index, columns the samples, 1st 1st replicate, 
            then 2nd replicate.
        .df_mother_sol: df containng the mother solution data, C_0
        .df_VoM_disol: pd series containing the volume [mL] added to the bottle 
        of the solution, BIC, or whatever. normally 50ml OR the total mass of 
        the solution [g]. If df_samples in ppb, this must be the total mass
        so that Q_e is in g/g !
        .df_m_bent: pd series contaning the mass of bentonite [g] in the bottle 
        (normally 250mg)
        .Nrepl: number of replicates. Default value = 2. 3 also accepted
        ret_Co__Ceq: if True, returns a df with C_0 - C_eq = False
        .removed_sample_1: string to say if you want to remove the 1st sample 
        (number 1) or not (True means remove),
            that it is the procedural blank. Default value: True (like I was 
                                                    doing before 7/24!)
    
    *Outputs (in that order):
        .df with the Kd data
        .df with q_e data
        
        '''
    
    
    ########## 0) Precalcs ##########
    #1st, in case of the data is not numeric, lets convert it to numeric
    df_VoM_disol =  df_VoM_disol.apply(pd.to_numeric) 
    df_m_be = df_m_be.apply(pd.to_numeric) 

    '''
    To avoid the div0 error, which occurs when I have 0 as a values in the df, 
    which I have for all the elements that were not found by ICPMS, I can just
    put NaN instead, since that will not give the Div0 error when computing Kd
    '''
    
    df_samples_aux = df_samples.copy()
    df_samples_aux.replace(0, np.nan, inplace=True)  #replace 0 values with NaN, 
                                              #to avoid Div0 error!    
    
    if df_m_liq== False:    #if no data for the volume of the liquid provided, 
                            #then set the mass of the liquid initially
        df_m_liq = df_VoM_disol - df_m_be 
    
    
    ########### 1) Calcs ###########
    '''
    The operations to perform are:
        1) C_0 - C_eq (>0)
        2) 1) * V/m = q_e
        3) 2)	1/C_eq = Kd
    
    I must treat the 2 experiments are different, I should substract the blank 1
    to the 1st emasurements and the 2 to the others. Since I ordered it in the 
    right way (1st replicacte 1, then replicate 2, I could) split it easily :D
            df.shape gives the shape of the df, n_rows, n_columns
    
    Note the df have number of samples * 2 replicates columns.
    
    Then, I will create a new dataframe substracting that data. To do so, I 
    need to get rid of the isotopes column, since is text, and then add it again.
    Watch, the substraction is easy with a pandas mehotd.

    I shuold then remove those columns
    from there, and replace negatives values for 0, for a good plot
    
    
    WELL, I MUST MODIFY EVERYTHING. I should do now
    1) C_0 * V/m
    2) C_eq * V_eq /m
    3) 1) - 2)
    4) 3) / C_e
                                                                                                      
                                                            
    
    '''
    
    #So, lets split into the 2 replicates!
    '''
    For 2 replicates its easy, for 3 it could be more tricky. Beware! TO create 
    a function you should say the number of replicates and so!
    '''
    
    #Gathering the replicates sepparately    
    df_1 = df_samples_aux.iloc[:, : round(df_samples_aux.shape[1] / 3)]
            #1st replicate
    df_2 = df_samples_aux.iloc[:, round(df_samples_aux.shape[1] / 3): 2*round(df_samples_aux.shape[1] / 3)]
    df_3 = df_samples_aux.iloc[:, 2*round(df_samples_aux.shape[1] / 3) :]
        
    df_VoM_1 = df_VoM_disol.iloc[ 0: round( ( df_VoM_disol.shape[0] ) / 3 ) ]     
    df_VoM_2 = df_VoM_disol.iloc[ round( ( df_VoM_disol.shape[0] ) / 3 ): 2* round( 
                    ( df_VoM_disol.shape[0] ) / 3 ) ]      
    df_VoM_3 = df_VoM_disol.iloc[ 2 *round( ( df_VoM_disol.shape[0] ) / 3 ) : ]      
        
    df_m_1 = df_m_be.iloc[ : round( ( df_m_be.shape[0] ) / 3 ) ]     
    df_m_2 = df_m_be.iloc[ round( ( df_m_be.shape[0] ) / 3 ) : 2* round( ( df_m_be.shape[0] ) / 3 )] 
    df_m_3 = df_m_be.iloc[ 2 *round( ( df_m_be.shape[0] ) / 3 ) : ]     
    
    #1) 
    '''
    Note that I have N different mother solutions, which also means N different
    samples. I could do that with a for loop, but I found a better version. 
    I can ubstrcat df ignoring their indexes by doing df.values!
    '''
    N = df_mother_sol.shape[1]       #number of samples in the mother solution
    
    dfCeq__C0_1 = pd.DataFrame(df_1.iloc[:,:N].values - df_mother_sol.values,
                               index = df_mother_sol.index, columns = df_1.columns) 
                                            ##doing the substraction checked!
    dfCeq__C0_2 = pd.DataFrame(df_2.iloc[:,:N].values - df_mother_sol.values,
                               index = df_mother_sol.index, columns = df_2.columns) 
    dfCeq__C0_3 = pd.DataFrame(df_3.iloc[:,:N].values - df_mother_sol.values,
                               index = df_mother_sol.index, columns = df_3.columns)
    
    if removed_sample_1 == True:
        '''
        If this is true, I do not want the 1st sample in each serie, the procedurla
        blank, so let´s remove it! how? by removing it in:
            .m/V data
            .Ceq-C0 data
            .Ceq data (Fro Kd calc only used)
        '''
    #Since 1st sample for each replicate is the procedural bllank, that 
    #IN PRINCIPLE is not needed, I remove it!!     
        dfCeq__C0_1.drop( [df_1.iloc[:,0].name], axis = 1, inplace = True)   
                #drop blank column
        dfCeq__C0_2.drop( [df_2.iloc[:,0].name], axis = 1, inplace = True)   
        dfCeq__C0_3.drop( [df_3.iloc[:,0].name], axis = 1, inplace = True)   
        #        
        df_m_1 = df_m_1[1:]  #fast way to delete 1st elemen (blank) in a series
        df_m_2 = df_m_2[1:]
        df_m_3 = df_m_3[1:]
        #
        df_VoM_1 = df_VoM_1[1:]
        df_VoM_2 = df_VoM_2[1:]
        df_VoM_3 = df_VoM_3[1:]
        #
        df_1.drop( [df_1.iloc[:,0].name], axis = 1, inplace = True)     
                            #removing sample 1 in repl 1
        df_2.drop( [df_2.iloc[:,0].name], axis = 1, inplace = True)
        df_3.drop( [df_3.iloc[:,0].name], axis = 1, inplace = True)
        print('####################################')
        print('Sample 1 in each replicate removed!!!!')
        print('####################################')
        #
    
    dfC0__Ceq_1 = - dfCeq__C0_1
    dfC0__Ceq_2 = - dfCeq__C0_2
    dfC0__Ceq_3 = - dfCeq__C0_3
    
    ######## 2) Apply the V/ m giving q_e (from Df_exp)
    #For this I ned to remove the blank columns to both m and V, since from 
    #C0-Ceq they are removed!
        
    #And now I can operate:
        
    df_Qe_1 = dfC0__Ceq_1 * df_VoM_1 / df_m_1
    df_Qe_2 = dfC0__Ceq_2 * df_VoM_2 / df_m_2 
    df_Qe_3 = dfC0__Ceq_3 * df_VoM_3 / df_m_3
        
        
    ######## 3) Apply 1/C_eq = Kd
    df_Kd_1 = df_Qe_1 / df_1  
    df_Kd_2 = df_Qe_2 / df_2
    df_Kd_3 = df_Qe_3 / df_3
    
    #Now lets add them together
    df_Kd = pd.concat( [df_Kd_1, df_Kd_2, df_Kd_3], axis = 1)         
    df_Qe = pd.concat( [df_Qe_1, df_Qe_2, df_Qe_3 ] , axis = 1)  
    df_C0__Ceq = pd.concat( [dfC0__Ceq_1, dfC0__Ceq_2, dfC0__Ceq_3], axis = 1)
    
    
   #Lets finally convert those variables into numeric ones, better
    df_Kd = df_Kd.apply(pd.to_numeric)       #converting it to numeric 
    df_Qe = df_Qe.apply(pd.to_numeric) 
    df_C0__Ceq =df_C0__Ceq.apply(pd.to_numeric)        
         
   ####Checked that the Qe calc is correct, and then Kd must be also :)     
      
    ########### 2) Return #############
    #Here the if for returning or not C_0 - C(t) applies
    
    if ret_Co__Ceq == False:            #do not return it
        return df_Kd, df_Qe             #return
    
    else:                   #return C0-Ceq
        return df_Kd, df_Qe, df_C0__Ceq



#%%####### 1.14) Mean/std of replicates calculator #############
#####################################


def ICPMS_MeanStd_calculator (df_data, Nrepl = 2):
    '''
    Function that will compute the mean and std of:
        .the measuring sequence (df). Note we measure 2 or 3 replicates,
        OR
        .a df series, like the shaking time
    so this function will ocmpute the average and sigma, ready to plot them :)

    If we have N measurements: x_1, x_2,..., x_N, the mean value is simply
        <x> = sum (x_i) /N , i=1,..N
    The variance (with respect to the mean value) is
        var = sigma**2 = sum(x_i-<x>)**2, i =1, ..N
    The std is 
        std = sqrt(var)/N-1 
    

    *Inputs:
        .df_data: dataframe containing the data, the full data, with the 2 
        or 3 replicates. Should be Dfs corrected
            Format: isotopes as index, columns the samples, 1st 1st replicate, 
            then 2nd replicate.
        .Nrepl: number of replicates. Default value = 2. 3 also accepted
    
    *Outputs (in that order):
        .df with the < >
        .df with the std
        Both df have same column names: Data + number. Data and not sample, 
        not to mix it, since the Data 1 may be from replicates 2 for ex, in 
        case of Qe, kd. But could be from replicates 1 for variables like 
        Conc for ex
        
    **TO Do:
        .Simplify, removing unncesary variables! There are many!!!!
        '''
    
    #df_data.replace(0, np.nan, inplace=True)     
                #replace 0 values with NaN, to avoid Div0 error!
    
    
    ########### 1) Calcs ###########
    '''
    The operations to perform are:
        1) < > of each sample number (1,2,..) 
        2) std if each sample number (Divideb by N-1, as excel do when suing
            STDEV() function! )
    
    Those cals are really easy since they are implemented in pandas. I just 
    need to sort the proper way of sepparating the df and getting the results, 
    see below. I need to discriminate whether the data is a df or a Series.
    
    '''
    
    if type(df_data) == pd.core.frame.DataFrame :         #If True data is a df

        df_mean = pd.DataFrame()         #Empty df to store the values
        df_std = pd.DataFrame()        #Empty df to store the std
    
        if Nrepl == 2:          #Standard case, 2 replicates
             df_1 = df_data.iloc[ :, 0: round( ( df_data.shape[1] ) / 2 ) ]      
                         #1st replicate
             df_2 = df_data.iloc[ :, round( ( df_data.shape[1] ) / 2 ) :  ]       
                         #replicate 2
        #
             for i in range(df_1.shape[1]):   
                 #loop thorugh all elements, but with index to work with 2 df
                 df_temp = df_data.iloc[:, [i, i+ df_1.shape[1] ] ]        
                         #df containing the 2 replicates of the number i
            
                #TO store temporarily those values I create an auxiliary df
                 df_aux1 = pd.DataFrame(data = df_temp.mean(1), columns = [i] )
                 df_aux2 = pd.DataFrame(data = df_temp.std(1), columns = [i] )

                #And I add that to the storing df (add as columns)
                 df_mean['Data ' + str(i +1) ] = df_aux1
                 df_std['Data ' + str(i +1) ] = df_aux2
        

        elif Nrepl == 3:            #3 replicates
        #Gathering the replicates sepparately    
            df_1 = df_data.iloc[:, : round(df_data.shape[1] / 3)]
            df_2 = df_data.iloc[:, round(df_data.shape[1] / 3): 2*round(df_data.shape[1] / 3)]
            df_3 = df_data.iloc[:, 2*round(df_data.shape[1] / 3) :]
        
            for i in range(df_1.shape[1]):         #loop thorugh all elements, but with index to work with 2 df
                df_temp = df_data.iloc[:, [i, i+ df_1.shape[1], i+ 2 * df_1.shape[1] ] ]        
                                                #df containing the 3 replicates of the number i
            
                #TO store temporarily those values I create an auxiliary df
                df_aux1 = pd.DataFrame(data = df_temp.mean(1), columns = [i] )      #< >. 1 indicates compute by columns
                df_aux2 = pd.DataFrame(data = df_temp.std(1), columns = [i] )       # Std 

                #And I add that to the storing df (add as columns)
                df_mean['Data ' + str(i +1) ] = df_aux1
                df_std['Data ' + str(i +1) ] = df_aux2

    elif type(df_data) == pd.core.frame.Series :      #if data is a serie
        df_mean = pd.Series()         #Empty df to store the values
        df_std = pd.Series()        #Empty df to store the std
        
        if Nrepl == 2:          #Standard case, 2 replicates
             df_1 = df_data.iloc[ 0: round( ( df_data.shape[0] ) / 2 ) ]      #1st replicate
             df_2 = df_data.iloc[ round( ( df_data.shape[0] ) / 2 ) :  ]       #replicate 2
            #
             for i in range(df_1.shape[0]):         #loop thorugh all elements, but with index to work with 2 df
                df_temp = df_data.iloc[[i, i+ df_1.shape[0] ] ]        #df containing the 2 replicates of the number i
                
                #TO store temporarily those values I create an auxiliary df
                #df_temp.mean and .std gives mean and std
                
                df_mean['Data ' + str(i +1 ) ] = df_temp.mean()      #<>
                df_std['Data ' + str(i +1 ) ] = df_temp.std()         #std

        
        elif Nrepl == 3:            #3 replicates
        #Gathering the replicates sepparately    
            df_1 = df_data.iloc[: round(df_data.shape[0] / 3)]
            df_2 = df_data.iloc[round(df_data.shape[0] / 3): 2*round(df_data.shape[0] / 3)]
            df_3 = df_data.iloc[ 2*round(df_data.shape[0] / 3) :]
            
            for i in range(df_1.shape[0]):         #loop thorugh all elements, but with index to work with 2 df
                df_temp = df_data.iloc[ [i, i+ df_1.shape[0], i+ 2 * df_1.shape[0] ] ]        
                                                    #df containing the 3 replicates of the number i

                #And I add that to the storing df (add as columns)
                df_mean['Data ' + str(i +1) ] = df_temp.mean()
                df_std['Data ' + str(i +1) ] = df_temp.std()    
    
    
    else:       #weird data
        print('############################')
        print('\n Bro WTF? Wrong data type! pd.Series or pd.DataFrame only!')
        print('############################')
        
        df_mean = False
        df_std = False
    
########### 2) Return #############
    return df_mean, df_std             #return



#%% ######## 1.15 Cs sep, cumulative conc computer ##############
#################################################################

def ICPMS_Cumulative_conc_calc(df_ppb, df_ppb_std, V, delta_V = 5,
                               rho = 1, delta_rho = .001, MS_here = True, 
                         ret_Nmass = False ):
    '''
    Function that will compute the cumulative concentration of a given nuclide.
    This assume mass not measured, so needs the volume and density to compute
    the mass of the sample, to compute the cumulative conc. 
    
    This is, imagine I have several samples, with a given [Cs133]. This code 
    will compute the total [Cs133], the concentration if all the samples are 
    merged into a single one. IN a succesive way, so:
        .sample 2 cumulative = sample 1 + smaple 2
        .sample 3 cumulative = sample 1 + sample 2 + sample 3
    Until the last one, which will contain the total cumulative concentration.
    The formula is:
    
    [X]_tot = m_X_tot/m_Tot = sum (m_X)/sum(m_tot) = sum([X]*m_tot) / sum (m_tot)
    = sum([X]*rho*V) / sum(V*Rho),
    
    being [X] the concentration of X, the nuclide of interest, m_tot the mass 
    of each sample, and m_Tot the total mass of all the solution, the sum of m_tot.
    There it is assumed we do not know the masses, rather volumes, so the calcs
    is based on the volumes, and densities. Volumes are given as input, as well
    as densities.
    
    *Input
        .df_ppb: df with the ppb values of the samples
        .df_ppb_std: df with the std of the ppb values
        .MS_here: boolean to indicate wether the MS is in the df_ppb
                or not. Default: True
        .rho: density, np array or single value. Default: 1g/mL
        .delta_rho: np array with rho uncertainty. Default: .001g/mL
        .V: np.array containing the volumes of the samples, or single value.
                eg: V+ np.array([100,100,100,100,100]) for 5 samples
        .delta_V: np array or single value, for the uncertainty of V. 
                Default: 5mL
        
        .ret_Nmass: boolean to indicate weather the df with the nuclide mass
            should ber etunred or not. Default: False, not return it
    *Output
        . df with the cumulative conc, for each sample
        .df with the std of the cumulative conc
        .df with the nuclide mass, if desired
        .df with the std of the nuclide mass, if desired
        
    '''
    
    ########### 0) Pre work ###############
    """
    The uncertainty of the volume was a number, but should be an array, but
    we can do it easily
    
    !!!!!!!
    In the future, I might need to do the same I do here, but for densities....
    !!!!!!!!!
    """
    delta_V = delta_V * np.ones(len(V))         #array creation
    
    ########### 1) Mass of nuclide (Nmass) calc ####################
    '''
    1st step, to calc m_X = [X] * m_tot = [X] * rho * V, for each sample
    '''
    if MS_here:     #exclude 1st column of the ppb df, the MS
        df_Nmass = df_ppb.iloc[:,1:]*V*rho
        df_Nmass_std = df_Nmass * np.sqrt((delta_V/V)**2 + (delta_rho/rho)**2 +
                            (df_ppb_std.iloc[:,1:]/df_ppb.iloc[:,1:])**2)
                #Uncertainty of Nmass
    else:       #MS not here
        df_Nmass = df_ppb*V*rho
        df_Nmass_std = df_Nmass * np.sqrt((delta_V/V)**2 + (delta_rho/rho)**2 +
                            (df_ppb_std/df_ppb)**2) 
    
    ############## 2) Cumulative concentration calc ###################
    '''
    Now I need to do the sums, 
    [X]_tot = sum([X]*rho*V) / sum(V*Rho)
    
    To do it, first I will create the df I will fill.
    '''
    df_cum = pd.DataFrame(index = df_Nmass.index, columns = df_Nmass.columns)
                #empty df creation, with desired column and row names :)
    df_cum_std = pd.DataFrame(index = df_Nmass.index, columns = df_Nmass.columns)
            #empty df, for the uncertainties!
            
    '''
    Now, I should fill that df, with a foor loop, or 2, as needed. Essentially 
    I need to do:
    [U238] = sum of df_ppb*m_U238, suming in columns, and that for each row (
        nuclide). note 1st column does not need to be sum, so I can compute them
    now:
    '''
 
    df_cum.iloc[:,0] = df_Nmass.iloc[:,0]/(V[0]*rho)    #1st row, non cumulative
    df_cum_std.iloc[:,0] = df_cum.iloc[:,0] * np.sqrt(
                (delta_V[0]/V[0])**2 + (delta_rho/rho)**2 )
                    #quadratic error prop, as always
    
    df_Nmass_std2 = df_Nmass_std**2         #square, will be needed for cum error
    
    for i in range(1,df_cum.shape[1]): #loop thourhg all columns of the df
        df_cum.iloc[:,i] = df_Nmass.iloc[:,0:i+1].sum(axis = 1) /(sum(V[0:i+1])*rho)
            #sum in columns the mass of the nuclide
        df_cum_std.iloc[:,i] = df_cum.iloc[:,i] * np.sqrt(
         (df_Nmass_std2.iloc[:,0:i+1].sum(axis = 1) /
          df_Nmass.iloc[:,0:i+1].sum(axis = 1)**2) + 
         sum(rho*V[0:i+1]* ((delta_rho/rho)**2 + (delta_V[0:i+1]/V[0:i+1])**2) ) / 
             ((sum(V[0:i+1])*rho)**2)
             ) 
                #Note that sum(V*rho)  =sum(V)*rho, since rho is a scalar!
                #Bro, fuck, what a mess to write that xD
                
    df_cum = df_cum.apply(pd.to_numeric)        #make it numeric, since it is not
    df_cum_std = df_cum_std.apply(pd.to_numeric) 
    
    ############ 3) Output ########   
    
    if ret_Nmass:           #also returns the std!
        return df_cum, df_cum_std, df_Nmass, df_Nmass_std
    else:       #Not return it
        return df_cum, df_cum_std

    
#%% ######## 1.16 Substraction of bentonite leached elements ##############
#################################################################

def ICPMS_Removal_Bent_leach(df_ppb, df_ppb_std, df_MS, df_MS_std,
                             return_leached = False, Elem_rel = Isot_rel,
                             N_repl = 3):
    '''
    This function will remove the nuclides leached by bentonite from the ICPMS
    data. Having the ppb data and the MS (with their errors), ASSUMING the
    1st sample of each replicate (3 replicates ASSUMED!) is the procedural blank,
    only bentonite and BIC water, and that the 1st MS is the blank, BIC water,
    This function, for each replicate, will substract the contribution of the
    bentonite to the concentration of all the elements.
    
    We have 3 replicates: 1_1, 1_2,.... 2_1,2_2,.., 3_1, 3_2, ...
    with their mother solutions; 0_1,0_2,.... Then, this function will compute 
    the elements leached by the bentonite for each replicate:
                    1_1-0_1
                    2_1-0_1
                    3_1-0_1
    And it will substract it to the other samples:
        1_2- (1_1-0_1), 1_3 - (1_1-0_1),...
        2_2 - (2_1-0_1), 2_3 - (2_1-0_1)...
        ...
        
    *INPUTS
        .df_ppb: df with the ppb data. It should have the replicates in order:
            1_1, 1_2,1_3,...,2_1,2_2,2_3,...,3_1,3_2,...
        .df_ppb_std: df with the std of the ppb data. Same order
        .df_MS: df with the mother solutions. In the order:
            0_1, 0_2,...
        .df_MS_std: df with the std of the MS
        .return_leached: boolean, to indicate wether df with the ppb (and 
            their std) leached by bentonite should be returned or not. Default:
            False
        .Elem_rel: array containing the name of the relevant nuclides.
        Eg: np.array(['U238(LR)', 'Sr88(LR)'])
        .N_repl: number of replicates. Default: 3
        
    *OUTPUTS
        .df_ppb_br: df with the ICPMS data with the bentonite contribution 
            removed (br). Note the procedural blank will still be there. Though
            they may not be useful, at least for the Qe function they will
        .df_ppb_br_std: df of the std of the ppb_br
    '''

    if N_repl == 3:         #Case of 3 replicates 
        ############### 1) Data preparation #######################
        #We need to separate the replicates, in order to perform the substraction
        df_1 = df_ppb.iloc[:, : round(df_ppb.shape[1] / 3)]
        #1st replicate
        df_2 = df_ppb.iloc[:, round(df_ppb.shape[1] / 3): 2*round(df_ppb.shape[1] / 3)]
        df_3 = df_ppb.iloc[:, 2*round(df_ppb.shape[1] / 3) :]    
        #Also their std are needed:    
        df_1_std = df_ppb_std.iloc[:, : round(df_ppb_std.shape[1] / 3)]
        df_2_std = df_ppb_std.iloc[:, round(df_ppb_std.shape[1] / 3): 2*round(df_ppb_std.shape[1] / 3)]
        df_3_std = df_ppb_std.iloc[:, 2*round(df_ppb_std.shape[1] / 3) :]      
        #
        ############## 2) Calc of the contribution of bentonite
        #This will be a df Series!
        Ser_ppb_leach_1 = df_1.iloc[:,0] - df_MS.iloc[:,0]
            #df series, with 1 column, named 0, and all indexes (nuclei)
        Ser_ppb_leach_2 = df_2.iloc[:,0] - df_MS.iloc[:,0]               
        Ser_ppb_leach_3 = df_3.iloc[:,0] - df_MS.iloc[:,0]
        #And their std:
        Ser_ppb_std_leach_1 = np.sqrt(df_1_std.iloc[:,0]**2 + df_MS_std.iloc[:,0]**2)
        Ser_ppb_std_leach_2 = np.sqrt(df_2_std.iloc[:,0]**2 + df_MS_std.iloc[:,0]**2)
        Ser_ppb_std_leach_3 = np.sqrt(df_3_std.iloc[:,0]**2 + df_MS_std.iloc[:,0]**2)
        '''
        Okay, this data I could print, and even store it and give it as an output,
        if desired.
    
        First I will join them in a df, and then printing it and saving it
        '''
        df_leached = pd.concat([Ser_ppb_leach_1, Ser_ppb_leach_2, 
                            Ser_ppb_leach_3], axis = 1)
        df_leached.columns = ['Repl 1','Repl 2','Repl 3']
        df_leached_std = pd.concat([Ser_ppb_std_leach_1, Ser_ppb_std_leach_2, 
                            Ser_ppb_std_leach_3], axis = 1)
        df_leached_std.columns = ['Repl 1','Repl 2','Repl 3']
    #I could print these data, or at least for the relevant Nuclei
    ############# 3) Substraction ##################
        df_1_br = df_1.subtract(Ser_ppb_leach_1.values, axis = 0) #substr in repl 1
        df_2_br = df_2.subtract(Ser_ppb_leach_2.values, axis = 0) 
        df_3_br = df_3.subtract(Ser_ppb_leach_3.values, axis = 0) 
        '''
    The uncertainty will be more complicate, since I need to do the sqrt of the
    sum of the squares. The simples thing, 
    np.sqrt(df_1_std**2 + Ser_ppb_std_leach_1**2)
    
    ofc does not work, because it treated the series as a row vector, not
    a column one. Asking chatgpt, it gave me a solution that seemed to work
    
        '''
        df_1_br_std = np.sqrt(df_1_std**2 + Ser_ppb_std_leach_1.values[:, None]**2)
        df_2_br_std = np.sqrt(df_2_std**2 + Ser_ppb_std_leach_2.values[:, None]**2)
        df_3_br_std = np.sqrt(df_3_std**2 + Ser_ppb_std_leach_3.values[:, None]**2)
        ################ 4) Output #############
        #First we need to mergue them
        df_ppb_br = pd.concat([df_1_br, df_2_br, df_3_br], axis = 1)
        df_ppb_std_br = pd.concat([df_1_br_std, df_2_br_std, df_3_br_std], axis = 1)
        #
        #
    elif N_repl == 2:           #2 replc
        df_1 = df_ppb.iloc[:, : round(df_ppb.shape[1] / 2)]   #1st replicate
        df_2 = df_ppb.iloc[:, 2*round(df_ppb.shape[1] / 2) :]    
        #Also their std are needed:    
        df_1_std = df_ppb_std.iloc[:, : round(df_ppb_std.shape[1] / 2)]
        df_2_std = df_ppb_std.iloc[:, 2*round(df_ppb_std.shape[1] / 2) :]      
        #
        ############## 2) Calc of the contribution of bentonite
        #This will be a df Series!
        Ser_ppb_leach_1 = df_1.iloc[:,0] - df_MS.iloc[:,0]
            #df series, with 1 column, named 0, and all indexes (nuclei)
        Ser_ppb_leach_2 = df_2.iloc[:,0] - df_MS.iloc[:,0]               
        #And their std:
        Ser_ppb_std_leach_1 = np.sqrt(df_1_std.iloc[:,0]**2 + df_MS_std.iloc[:,0]**2)
        Ser_ppb_std_leach_2 = np.sqrt(df_2_std.iloc[:,0]**2 + df_MS_std.iloc[:,0]**2)
        #
        df_leached = pd.concat([Ser_ppb_leach_1, Ser_ppb_leach_2], axis = 1)
        df_leached.columns = ['Repl 1','Repl 2']
        df_leached_std = pd.concat([Ser_ppb_std_leach_1, Ser_ppb_std_leach_2],
                                   axis = 1)
        df_leached_std.columns = ['Repl 1','Repl 2']
    #I could print these data, or at least for the relevant Nuclei
    ############# 3) Substraction ##################
        df_1_br = df_1.subtract(Ser_ppb_leach_1.values, axis = 0) #substr in repl 1
        df_2_br = df_2.subtract(Ser_ppb_leach_2.values, axis = 0) 
        #
        df_1_br_std = np.sqrt(df_1_std**2 + Ser_ppb_std_leach_1.values[:, None]**2)
        df_2_br_std = np.sqrt(df_2_std**2 + Ser_ppb_std_leach_2.values[:, None]**2)   
        #
        #
        df_ppb_br = pd.concat([df_1_br, df_2_br], axis = 1)
        df_ppb_std_br = pd.concat([df_1_br_std, df_2_br_std], axis = 1)
    else:
        print('Number of replicates different from 2 or 3? nothing done bro')

    
    ###Note that those include the procedural blank, shall I remove it?
    #Nooo, because it is needed for the Qe calc, since I did the funciton in
    #the way I did. SO just keep it!
    
    #After the loop, I can print that
    print('####### Concentration in ppb of leached relevant nuclides from the bentonites to the BIC solution')
    print(df_leached.loc[Elem_rel])
    print('######################################\n')

    print('########## Uncertainty of those: ########')
    print(df_leached_std.loc[Elem_rel])
    print('######################################\n')
    
    ########### Output ###################
    
    if return_leached:      #True, so return it
        return df_ppb_br, df_ppb_std_br, df_leached, df_leached_std
    else:       #False, dont return it
        return df_ppb_br, df_ppb_std_br
    
    
#%% ######## 1.17 Ratio based correction Bentonite leached ##############
#################################################################    

def ICPMS_Removal_Bent_leach_ratio(df_ppb, df_ppb_std, df_MS, df_MS_std,
                             return_leached = False, Elem_rel = Isot_rel,
                             N_repl = 3):
    '''
    This function will remove the nuclides leached by bentonite from the ICPMS
    data, but in a ratio based way. Having the ppb data and the MS 
    (with their errors), ASSUMING the
    1st sample of each replicate (3 replicates ASSUMED!) is the procedural blank,
    only bentonite and BIC water, and that the 1st MS is the blank, BIC water,
    This function, for each replicate, will substract the contribution of the
    bentonite to the concentration of all the elements.
    
    We have 3 replicates: 1_1, 1_2,.... 2_1,2_2,.., 3_1, 3_2, ...
    with their mother solutions; 0_1,0_2,.... Then, this function will compute 
    the elements leached by the bentonite (C_leach) for each replicate:
                    1_1-0_1
                    2_1-0_1
                    3_1-0_1
    And it will perform the following correction to the other samples:
        C_f corr = C_f * C_0/(C_0+C_leach)
    which would be like:
        1_2 corr = 1_2 * 0_1/(0_1+C_leach 1), 1_3 corr = 1_3 * 0_1/(0_1+C_leach 1),
        ...
        2_2 corr = 2_2 * 0_2/(0_2+C_leach), 2_3 cor = 2_3 * 0_2/(0_2+C_leac)
        ...
        
    *INPUTS
        .df_ppb: df with the ppb data. It should have the replicates in order:
            1_1, 1_2,1_3,...,2_1,2_2,2_3,...,3_1,3_2,...
        .df_ppb_std: df with the std of the ppb data. Same order
        .df_MS: df with the mother solutions. In the order:
            0_1, 0_2,...
        .df_MS_std: df with the std of the MS
        .return_leached: boolean, to indicate wether df with the ppb (and 
            their std) leached by bentonite should be returned or not. Default:
            False
        .Elem_rel: array containing the name of the relevant nuclides.
        Eg: np.array(['U238(LR)', 'Sr88(LR)'])
        .N_repl: number of replicates. Default: 3
        
    *OUTPUTS
        .df_ppb_br: df with the ICPMS data with the bentonite contribution 
            removed (br). Note the procedural blank will still be there. Though
            they may not be useful, at least for the Qe function they will
        .df_ppb_br_std: df of the std of the ppb_br
    '''

    if N_repl == 3:         #Case of 3 replicates 
        ############### 1) Data preparation #######################
        #We need to separate the replicates, in order to perform the substraction
        df_1 = df_ppb.iloc[:, : round(df_ppb.shape[1] / 3)]
        #1st replicate
        df_2 = df_ppb.iloc[:, round(df_ppb.shape[1] / 3): 2*round(df_ppb.shape[1] / 3)]
        df_3 = df_ppb.iloc[:, 2*round(df_ppb.shape[1] / 3) :]    
        #Also their std are needed:    
        df_1_std = df_ppb_std.iloc[:, : round(df_ppb_std.shape[1] / 3)]
        df_2_std = df_ppb_std.iloc[:, round(df_ppb_std.shape[1] / 3): 2*round(df_ppb_std.shape[1] / 3)]
        df_3_std = df_ppb_std.iloc[:, 2*round(df_ppb_std.shape[1] / 3) :]      
        #
        ############## 2) Calc of the contribution of bentonite
        #This will be a df Series!
        Ser_ppb_leach_1 = df_1.iloc[:,0] - df_MS.iloc[:,0]
            #df series, with 1 column, named 0, and all indexes (nuclei)
        Ser_ppb_leach_2 = df_2.iloc[:,0] - df_MS.iloc[:,0]               
        Ser_ppb_leach_3 = df_3.iloc[:,0] - df_MS.iloc[:,0]
        #And their std:
        Ser_ppb_std_leach_1 = np.sqrt(df_1_std.iloc[:,0]**2 + df_MS_std.iloc[:,0]**2)
        Ser_ppb_std_leach_2 = np.sqrt(df_2_std.iloc[:,0]**2 + df_MS_std.iloc[:,0]**2)
        Ser_ppb_std_leach_3 = np.sqrt(df_3_std.iloc[:,0]**2 + df_MS_std.iloc[:,0]**2)
        '''
        Okay, this data I could print, and even store it and give it as an output,
        if desired.
    
        First I will join them in a df, and then printing it and saving it
        '''
        df_leached = pd.concat([Ser_ppb_leach_1, Ser_ppb_leach_2, 
                            Ser_ppb_leach_3], axis = 1)
        df_leached.columns = ['Repl 1','Repl 2','Repl 3']
        df_leached_std = pd.concat([Ser_ppb_std_leach_1, Ser_ppb_std_leach_2, 
                            Ser_ppb_std_leach_3], axis = 1)
        df_leached_std.columns = ['Repl 1','Repl 2','Repl 3']
        #This data will be printed at the end of the function..
    ############# 3) Correction ##################
        #I need fist to compute the correction factor:
        corr_1 = df_MS/(df_MS + Ser_ppb_leach_1.values[:, None]) 
                #correction for replicate 1
        corr_2 = df_MS/(df_MS + Ser_ppb_leach_2.values[:, None])
        corr_3 = df_MS/(df_MS + Ser_ppb_leach_3.values[:, None])        
        Delta_corr_1 = corr_1*np.sqrt(
            (df_MS_std/df_MS)**2 + (df_MS_std**2 + Ser_ppb_std_leach_1.values[:, None]**2)/(
            df_MS**2+Ser_ppb_leach_1.values[:, None]**2) ) #uncertainty of the corr
        Delta_corr_2 = corr_2*np.sqrt(
            (df_MS_std/df_MS)**2 + (df_MS_std**2 + Ser_ppb_std_leach_2.values[:, None]**2)/(
            df_MS**2+Ser_ppb_leach_2.values[:, None]**2) )
        Delta_corr_3 = corr_3*np.sqrt(
            (df_MS_std/df_MS)**2 + (df_MS_std**2 + Ser_ppb_std_leach_3.values[:, None]**2)/(
            df_MS**2+Ser_ppb_leach_3.values[:, None]**2) )
        #
        #Applying them:
        df_1_br = df_1*corr_1.values
        df_2_br = df_2*corr_2.values
        df_3_br = df_3*corr_3.values
        df_1_br_std = df_1_br*np.sqrt( (df_1_std/df_1)**2 + (Delta_corr_1/corr_1).values**2)
        df_2_br_std = df_2_br*np.sqrt( (df_2_std/df_2)**2 + (Delta_corr_2/corr_2).values**2)
        df_3_br_std = df_1_br*np.sqrt( (df_3_std/df_3)**2 + (Delta_corr_3/corr_3).values**2)
        '''
    The uncertainty will be more complicate, since I need to do the sqrt of the
    sum of the squares. The simples thing, 
    np.sqrt(df_1_std**2 + Ser_ppb_std_leach_1**2)
    
    ofc does not work, because it treated the series as a row vector, not
    a column one. Asking chatgpt, it gave me a solution that seemed to work
    
        '''
        ################ 4) Output #############
        #First we need to mergue them
        df_ppb_br = pd.concat([df_1_br, df_2_br, df_3_br], axis = 1)
        df_ppb_std_br = pd.concat([df_1_br_std, df_2_br_std, df_3_br_std], axis = 1)
        #
        #
    elif N_repl == 2:           #2 replc
        df_1 = df_ppb.iloc[:, : round(df_ppb.shape[1] / 2)]   #1st replicate
        df_2 = df_ppb.iloc[:, 2*round(df_ppb.shape[1] / 2) :]    
        #Also their std are needed:    
        df_1_std = df_ppb_std.iloc[:, : round(df_ppb_std.shape[1] / 2)]
        df_2_std = df_ppb_std.iloc[:, 2*round(df_ppb_std.shape[1] / 2) :]      
        #
        ############## 2) Calc of the contribution of bentonite
        #This will be a df Series!
        Ser_ppb_leach_1 = df_1.iloc[:,0] - df_MS.iloc[:,0]
            #df series, with 1 column, named 0, and all indexes (nuclei)
        Ser_ppb_leach_2 = df_2.iloc[:,0] - df_MS.iloc[:,0]               
        #And their std:
        Ser_ppb_std_leach_1 = np.sqrt(df_1_std.iloc[:,0]**2 + df_MS_std.iloc[:,0]**2)
        Ser_ppb_std_leach_2 = np.sqrt(df_2_std.iloc[:,0]**2 + df_MS_std.iloc[:,0]**2)
        #
        df_leached = pd.concat([Ser_ppb_leach_1, Ser_ppb_leach_2], axis = 1)
        df_leached.columns = ['Repl 1','Repl 2']
        df_leached_std = pd.concat([Ser_ppb_std_leach_1, Ser_ppb_std_leach_2],
                                   axis = 1)
        df_leached_std.columns = ['Repl 1','Repl 2']
    #I could print these data, or at least for the relevant Nuclei
    ############# 3) Correction ##################
        #I need fist to compute the correction factor:
        corr_1 = df_MS/(df_MS + Ser_ppb_leach_1.values[:, None]) 
                #correction for replicate 1
        corr_2 = df_MS/(df_MS + Ser_ppb_leach_2.values[:, None])     
        Delta_corr_1 = corr_1*np.sqrt(
            (df_MS_std/df_MS)**2 + (df_MS_std**2 + Ser_ppb_std_leach_1.values[:, None]**2)/(
            df_MS**2+Ser_ppb_leach_1.values[:, None]**2) ) #uncertainty of the corr
        Delta_corr_2 = corr_2*np.sqrt(
            (df_MS_std/df_MS)**2 + (df_MS_std**2 + Ser_ppb_std_leach_2.values[:, None]**2)/(
            df_MS**2+Ser_ppb_leach_2.values[:, None]**2) )
        #
        #Applying them:
        df_1_br = df_1*corr_1.values
        df_2_br = df_2*corr_2.values
        df_1_br_std = df_1_br*np.sqrt( (df_1_std/df_1)**2 + (Delta_corr_1/corr_1).values**2)
        df_2_br_std = df_2_br*np.sqrt( (df_2_std/df_2)**2 + (Delta_corr_2/corr_2).values**2)
    else:
        print('Number of replicates different from 2 or 3? nothing done bro!')

    
    ###Note that those include the procedural blank, shall I remove it?
    #Nooo, because it is needed for the Qe calc, since I did the funciton in
    #the way I did. SO just keep it!
    
    #After the loop, I can print that
    print('####### Concentration in ppb of leached relevant nuclides from the bentonites to the BIC solution')
    print(df_leached.loc[Elem_rel])
    print('######################################\n')

    print('########## Uncertainty of those: ########')
    print(df_leached_std.loc[Elem_rel])
    print('######################################\n')
    
    ########### Output ###################
    
    if return_leached:      #True, so return it
        return df_ppb_br, df_ppb_std_br, df_leached, df_leached_std
    else:       #False, dont return it
        return df_ppb_br, df_ppb_std_br





# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
####################### PLOTTERS ####################################
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#%% ########## 1.16) ICPMS Single Bar plotter #############
#####################################

def ICPMS_1Barplotter (df_1, df_2, ylabel_1 = 'I [cps]' , folder_name = 'Bar_plots',
                      pre_title_plt = "Concentration of ", 
                      pre_save_name = 'Conc_rsd', Elem_rel = Isot_rel, 
                      plot_everything = False, Logs = False ):
    '''
    Function that will do a single bar plots of ICPMS data, using
    the 2nd data as the errorbars. By raw I mean withoutany corrections/
    calibrations (neither the dilution factor correction). This is a preliminary 
    plot, to check if  everything is allright, with the rstd (rstd should be < 20%).
    
    *Inputs:
        .df_1, df_2: dataframes to plot. Initially containing the cps and 
        the std. 
        df_1 is a DataFrame, df_2 could be a DataFrame or a df.Series!
        Note the Isotopes column are the index for the df
        .folder_name: folder name. Default value: 'Bar_plots'
        .ylabel_1,2. Labels to be used for the legend and the y axes. Default: 
            ylabel_1 = 'I [cps]' , ylabel_2 = "$\sigma_{rel}$ [%]"
        .pre_title_plt : title of the graph, part that appears before the name 
        of the elements (thats why pre title).
            Detault value: "Concentration of " (note the space after of, so 
            the element is not together with that!)
        . pre_save_name: name of the graph files to save. Default: 'Conc', 
        giving Conc_Mg24.png for ex           
        .Elem_rel: array containing the name of the relevant elemtns, which 
        are the elements that will be saved in a specific folder. D
        efault value: (see above in the script)   
        .Plot_everything: boolean stating if we plot everything or only the 
        relevant elements (in Elem_rel). 
            Default: False
        .logs: value defining if applying log scale to y axis. Default: False
        
    *Outputs:
        .Plots (saving them) of the raw data
    
    ######### To DO #########
    .TInclude option to plot all the elements or only the relevants, since 
    for all (250) take 140s!
    
    '''
    
    
    ############# 1) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder a 
    subfolder with the relevant elements, to be given, will be created
    '''
    
    path_bar_pl = os.getcwd() + '/' + folder_name + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)

    #Subfolder with relevant plots:
    path_bar_pl_rel = os.getcwd() + '/' + folder_name + '/' + 'Relevants' + '/' 
        #folder path for the relevant plots
    if not os.path.exists(path_bar_pl_rel):
        os.makedirs(path_bar_pl_rel)   
    
    
    ######### 2) plotting ###############
    '''
    This is a loop plot, so beware, will take long (2-3mins!).
    '''
    t_start = tr.time()       #[s] start time of the plot execution
    
    ###Plot
    #Some parameters for the plot
    X_axis = np.arange( len(df_1.axes[1]))                 #To do the 2 bar plot
            #choosing the number of columns. [0] for rows

    for i in list( range(df_1.shape[0] ) ):     #Loop thorugh all rows
                    #df_1.index give the index values, low and high
                    #
        if df_1.index[i] in Elem_rel or plot_everything == True: 
            #if the element is relevant you plot it!
        ########### Bar plot ###########################
            plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
            plt.title(pre_title_plt + df_1.index[i][:-4], fontsize=22, wrap=True)           #title
            if isinstance(df_2, pd.Series):                         #if df_2 is a pd.Series
                aa =  plt.bar( X_axis, df_1.loc[df_1.index[i]], yerr = df_2, edgecolor="black", 
                               align='center') 
            else:               #df_2 is a DataFrame
                aa = plt.bar( X_axis, df_1.loc[df_1.index[i]], yerr = df_2.loc[df_1.index[i]], edgecolor="black", 
                             align='center')  ##df version!! works!
                #
            if Logs == True:                  #put scale
                plt.yscale('log')  
            #
            plt.ylabel(ylabel_1, fontsize= Font)              #ylabel
            plt.xlabel('Sample', fontsize = Font )
            plt.tick_params(axis='both', labelsize= Font)              #size of axis
            plt.minorticks_on()             #enabling minor grid lines
            plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
            plt.grid(which = 'major')
            plt.xticks(X_axis, [ df_1.columns[:][j] for j in range(len(df_1.axes[1]) ) ]
                , rotation = 90)
        
            #Saving in the folder
            if df_1.index[i] in Elem_rel:  #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
                plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_' + df_1.index[i][:-4] + '.png', format='png', bbox_inches='tight')
            #
            else:        #if the element is not relevant
                plt.savefig(folder_name +'/' +  
                        pre_save_name + '_' + df_1.index[i][:-4] +'.png', format='png', bbox_inches='tight')
                    #To save plot in folder
            plt.close()             #to clsoe the plot not to consume too much resources
    
    
    ######### 3) Running time displaying ###############
    '''
    The last thing will be to see and display the time needed
    '''
    
    t_run = tr.time() - t_start     #Running time

    print('###############################################')
    print('Plotting running time: ' + str(t_run) + 's')
    print('###############################################')


#%% ########## 1.17) ICPMS Single Bar plotter #############
#####################################

def ICPMS_1Bar_1line_plotter (df_1, df_2, df_3, ylabel_1 = 'I [cps]' , 
                              folder_name = 'Bar_plots',
                      pre_title_plt = "Concentration of ", 
                      pre_save_name = 'Conc_rsd', Elem_rel = Isot_rel, 
                      plot_everything = False, Logs = False ):
    '''
    Function that will do a single bar plots of ICPMS data, using
    the 2nd data as the errorbars. This function inlcudes a hline plot, with the
    data from another df, df_3. It will plot the 1st column of that df.
    
    Function created for the Cs sep plotting
    
    *Inputs:
        .df_1, df_2: dataframes to plot. Initially containing the cps and 
        the std. 
        .df_3: df from the hline
        df_1 is a DataFrame, df_2 could be a DataFrame or a df.Series!
        Note the Isotopes column are the index for the df
        .folder_name: folder name. Default value: 'Bar_plots'
        .ylabel_1,2. Labels to be used for the legend and the y axes. Default: 
            ylabel_1 = 'I [cps]' , ylabel_2 = "$\sigma_{rel}$ [%]"
        .pre_title_plt : title of the graph, part that appears before the name 
        of the elements (thats why pre title).
            Detault value: "Concentration of " (note the space after of, so 
            the element is not together with that!)
        . pre_save_name: name of the graph files to save. Default: 'Conc', 
        giving Conc_Mg24.png for ex           
        .Elem_rel: array containing the name of the relevant elemtns, which 
        are the elements that will be saved in a specific folder. D
        efault value: (see above in the script)   
        .Plot_everything: boolean stating if we plot everything or only the 
        relevant elements (in Elem_rel). 
            Default: False
        .logs: value defining if applying log scale to y axis. Default: False
        
    *Outputs:
        .Plots (saving them) of the raw data
    
    ######### To DO #########
    .Generalize it, so it can plot a df series for hline, not only df!!!!!
    
    '''
    
    
    ############# 1) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder a 
    subfolder with the relevant elements, to be given, will be created
    '''
    
    path_bar_pl = os.getcwd() + '/' + folder_name + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)

    #Subfolder with relevant plots:
    path_bar_pl_rel = os.getcwd() + '/' + folder_name + '/' + 'Relevants' + '/' 
        #folder path for the relevant plots
    if not os.path.exists(path_bar_pl_rel):
        os.makedirs(path_bar_pl_rel)   
    
    
    ######### 2) plotting ###############
    '''
    This is a loop plot, so beware, will take long (2-3mins!).
    '''
    t_start = tr.time()       #[s] start time of the plot execution
    
    ###Plot
    #Some parameters for the plot
    X_axis = np.arange( len(df_1.axes[1]))                 #To do the 2 bar plot
            #choosing the number of columns. [0] for rows

    for i in list( range(df_1.shape[0] ) ):     #Loop thorugh all rows
                    #df_1.index give the index values, low and high
                    #
        if df_1.index[i] in Elem_rel or plot_everything == True: 
            #if the element is relevant you plot it!
        ########### Bar plot ###########################
            plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
            plt.title(pre_title_plt + df_1.index[i][:-4], fontsize=22, wrap=True)           #title
            if isinstance(df_2, pd.Series):                         #if df_2 is a pd.Series
                aa =  plt.bar( X_axis, df_1.loc[df_1.index[i]], yerr = df_2, edgecolor="black", 
                               align='center', label = 'Samples') 
            else:               #df_2 is a DataFrame
                aa = plt.bar( X_axis, df_1.loc[df_1.index[i]], yerr = df_2.loc[df_1.index[i]], edgecolor="black", 
                             align='center', label = 'Samples')  ##df version!! works!
                #
            plt.hlines(df_3.loc[df_3.index[i]][0], min(X_axis), max(X_axis),
                       label = "MS", color = 'r')   #PLot the hline!
            if Logs == True:                  #put scale
                plt.yscale('log')  
            #
            plt.ylabel(ylabel_1, fontsize= Font)              #ylabel
            plt.xlabel('Sample', fontsize = Font )
            plt.tick_params(axis='both', labelsize= Font)              #size of axis
            plt.minorticks_on()             #enabling minor grid lines
            plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
            plt.grid(which = 'major')
            plt.xticks(X_axis, [ df_1.columns[:][j] for j in range(len(df_1.axes[1]) ) ]
                , rotation = 90)
            plt.legend(fontsize = Font)
        
            #Saving in the folder
            if df_1.index[i] in Elem_rel:  #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
                plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_' + df_1.index[i][:-4] + '.png', format='png', bbox_inches='tight')
            #
            else:        #if the element is not relevant
                plt.savefig(folder_name +'/' +  
                        pre_save_name + '_' + df_1.index[i][:-4] +'.png', format='png', bbox_inches='tight')
                    #To save plot in folder
            plt.close()             #to clsoe the plot not to consume too much resources
    
    
    ######### 3) Running time displaying ###############
    '''
    The last thing will be to see and display the time needed
    '''
    
    t_run = tr.time() - t_start     #Running time

    print('###############################################')
    print('Plotting running time: ' + str(t_run) + 's')
    print('###############################################')


#%% ########## 1.15) ICPMS 2 Bar plotter #############
#####################################

def ICPMS_Barplotter (df_1, df_2, ylabel_1 = 'I [cps]' , 
                      ylabel_2 = "$\sigma_{rel}$ [%]", folder_name = 'Bar_plots',
                      pre_title_plt = "Concentration of ", 
                      pre_save_name = 'Conc_rsd', Elem_rel = Isot_rel, 
                      plot_everything = False, Logs = 0 ):
    '''
    Function that will do bar plots of the raw data from the ICPMS, the cps and
    the rstd. By raw I mean withoutany corrections/calibrations (neither the 
    dilution factor correction). This is a preliminary plot, to check if 
    everything is allright, with the rstd (rstd should be < 20%).
    
    *Inputs:
        .df_1, df_2: dataframes to plot. Initially containing the cps and 
        the relative standard deviation. 
        df_1 is a DataFrame, df_2 could be a DataFrame or a df.Series!
        Note the Isotopes column are the index for the df
        .folder_name: folder name. Default value: 'Bar_plots'
        .ylabel_1,2. Labels to be used for the legend and the y axes. Default: 
            ylabel_1 = 'I [cps]' , ylabel_2 = "$\sigma_{rel}$ [%]"
        .pre_title_plt : title of the graph, part that appears before the name 
        of the elements (thats why pre title).
            Detault value: "Concentration of " (note the space after of, so 
            the element is not together with that!)
        . pre_save_name: name of the graph files to save. Default: 'Conc', 
        giving Conc_Mg24.png for ex           
        .Elem_rel: array containing the name of the relevant elemtns, which 
        are the elements that will be saved in a specific folder. D
        efault value: (see above in the script)   
        .Plot_everything: boolean stating if we plot everything or only the 
        relevant elements (in Elem_rel). 
            Default: False
        .logs: value defining if applying log scale to 1, 2, or both. 
            0: no apply
            1: apply to df 1. 
            2: apply to df 2. 
            3: apply to both
        
    *Outputs:
        .Plots (saving them) of the raw data
    
    ######### To DO #########
    .TInclude option to plot all the elements or only the relevants, since 
    for all (250) take 140s!
    
    '''
    
    
    ############# 1) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder a 
    subfolder with the relevant elements, to be given, will be created
    '''
    
    path_bar_pl = os.getcwd() + '/' + folder_name + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)

    #Subfolder with relevant plots:
    path_bar_pl_rel = os.getcwd() + '/' + folder_name + '/' + 'Relevants' + '/' 
        #folder path for the relevant plots
    if not os.path.exists(path_bar_pl_rel):
        os.makedirs(path_bar_pl_rel)   
    
    
    ######### 2) plotting ###############
    '''
    This is a loop plot, so beware, will take long (2-3mins!).
    
    
    ##### Mutiple bar plot; number of bars = 2 ############

The relations that must be satisfied are, if the width of the bar is w and the 
blank space between values (value is x value) is b<1, is:

        centroid = value - w/2
        centroid (right of value) = value + w/2
        2w + b = 1 ==> w = (1-b)/2
        
Setting b gives w. In fact the general equations for 2n bars per
 X tick (n = 1,2,..) are:
    w = (1-b)/2n
    value - w/2*n <= centroid <= value + w/2*n, in steps of w (of course)
    
    ##############################
    '''
    t_start = tr.time()       #[s] start time of the plot execution
    
    ###Plot
    #Some parameters for the plot
    X_axis = np.arange( len(df_1.axes[1]))                 #To do the 2 bar plot
            #choosing the number of columns. [0] for rows
    b = .4                              #[au] blank space between values <1
    w = (1-b)/2          #bar width

    for i in list( range(df_1.shape[0] ) ):     #Loop thorugh all rows
                    #df_1.index give the index values, low and high
                    #
        if df_1.index[i] in Elem_rel or plot_everything == True: 
            #if the element is relevant you plot it!
        ########### Bar plot ###########################
            plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
            plt.title(pre_title_plt + df_1.index[i][:-4], fontsize=22, wrap=True)           #title
            a = plt.bar(X_axis - w/2, df_1.loc[df_1.index[i]], width = w, edgecolor="black", 
                        label = ylabel_1, align='center') 
            #-2 not to plot the blank!! Remove it to plot it!
            plt.ylabel(ylabel_1, fontsize= Font)              #ylabel
            plt.xlabel('Sample', fontsize = Font )
            plt.tick_params(axis='both', labelsize= Font)              #size of axis
            if Logs ==1 or Logs == 3:           #put log scale
                plt.yscale('log') 
            plt.minorticks_on()             #enabling minor grid lines
            plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
            plt.grid(which = 'major')
            plt.xticks(X_axis, [ df_1.columns[:][j] for j in range(len(df_1.axes[1]) ) ]
                , rotation = 90)
            plt.twinx()             #For setting 2 axes
            if isinstance(df_2, pd.Series):                         #if df_2 is a pd.Series
                aa =  plt.bar( X_axis + w/2, df_2, width = w, edgecolor="black", 
                              label = ylabel_2, align='center', color = 'red') 
            else:               #df_2 is a DataFrame
                aa = plt.bar( X_axis + w/2, df_2.loc[df_1.index[i]], width = w, edgecolor="black", 
                             label = ylabel_2, align='center', color = 'red')  ##df version!! works!
                #
            if Logs == 2 or Logs == 3:                  #put scale
                plt.yscale('log')  
            plt.ylabel(ylabel_2, fontsize= Font)              #ylabel
            #
            aaa = [a, aa]
            plt.legend(aaa, [p_.get_label() for p_ in aaa])
        
            #Saving in the folder
            if df_1.index[i] in Elem_rel:  #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
                plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_' + df_1.index[i][:-4] + '.png', format='png', bbox_inches='tight')
            #
            else:        #if the element is not relevant
                plt.savefig(folder_name +'/' +  
                        pre_save_name + '_' + df_1.index[i][:-4] +'.png', format='png', bbox_inches='tight')
                    #To save plot in folder
            plt.close()             #to clsoe the plot not to consume too much resources
    
    
    ######### 3) Running time displaying ###############
    '''
    The last thing will be to see and display the time needed
    '''
    
    t_run = tr.time() - t_start     #Running time

    print('###############################################')
    print('Plotting running time: ' + str(t_run) + 's')
    print('###############################################')
    
    

#%%######## 1.16) ICPMS plotter #############
#####################################

def ICPMS_Plotter (x, df_cps, x_label, y_label, folder_name = 'Plots', 
                   pre_title_plt = "Concentration of ", pre_save_name = 'Conc',
                   Elem_rel = Isot_rel, 
                   plot_everything = False, Nrepl = 2 ):
    '''
    Function that will plots of the data from the ICPMS (cps) vs another variable,
    initially time, the cps and the rstd. This assume we have 2 replicates, 
    1 series after the other. Stimated running time around 80s.
    
    *Inputs:
        .x: x axis variable in the plot. This should be a df series
        .df_cps: dataframes containing the cps. Those are
        outputs for the Read_ICPMS_excel function. Note the 1st column must be 
        the one with the isotopes (like in the excel)
        .x_label: string that will be the x label for the plot (for math stuff, 
                                    use $$. eg: '$\Delta t[h]$')
        .y_label: string that will be the y label for the plot
        .folder_name: string defining the name of the folder to create to store
        the plots. default value: 'Plots'
        . plot_everything: string defining if you want to plot all the elements 
        or only the relevant ones. Default value: False (only plot relevants)
        .pre_title_plt : title of the graph, part that appears before the name 
        of the elements (thats why pre title).
            Detault value: "Concentration of " (note the space after of, so 
                the element is not together with that!)
        . pre_save_name: name of the graph files to save. Default: 'Conc', 
        giving Conc_Mg24.png for ex    
        . Nrepl: number of replicates. Default value: 2. 3 also accepted
        .Elem_rel: array containing the name of the relevant elemtns, which 
            are the elements that will be saved in a specific folder. 
            Default value: (see above in the script)     
                                    
    *Outputs:
        .Plots (saving them) of the x and df_cps data, cps vs x!
    
    
    ### TO DO: ####
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	.Implement error plotting (in an errorbar pyplot)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    '''
    
    
    ############# 1) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder a
    subfolder
    with the relevant elements, to be given, will be created
    '''
    
    path_bar_pl = os.getcwd() + '/' + folder_name + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)

    #Subfolder with relevant plots:
    path_bar_pl_rel = os.getcwd() + '/' + folder_name + '/' + 'Relevants' + '/' 
        #folder path for the relevant plots
    
    if not os.path.exists(path_bar_pl_rel):
        os.makedirs(path_bar_pl_rel)   
    
    
    ######### 2) plotting ###############
    '''
    This is a loop plot, so beware, will take long if you plot all the elements 
    (280) (2-3mins!). THe best way is loop with python style
    '''
    t_start = tr.time()       #[s] start time of the plot execution
    
    '''
    Before plotting, I need to see whether I have 2 or 3 replicates, to divide 
    the df in one or another way. That will be with an if loop, as usual
    '''
    
    if Nrepl == 2:          #2 replicates, standard case
    
    ###Plot
        for index, row in df_cps.iterrows():   #Loop throught all the rows, the isotopes
                    #df_cps.index give the index values, low and high
		   # 4 because of the way the df is created (and hence the excel tabelle)

            #Saving in the folder
            if index in Elem_rel:  #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
                plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                plt.title(pre_title_plt + index[:-4], fontsize=22, wrap=True)           #title
                plt.plot(x[:int(len(x)/2)], row[:int(len(x)/2)], 'bo--', 
                     markersize = 5, label = 'Repl_1') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x[int(len(x)/2):], row[int(len(x)/2) :], 'ro--', 
                     markersize = 5, label = 'Repl_2') 
                plt.ylabel(y_label, fontsize= Font)              #ylabel
                plt.xlabel(x_label, fontsize = Font)
                plt.tick_params(axis='both', labelsize= Font)              #size of axis
                #plt.yscale('log') 
                plt.minorticks_on()             #enabling minor grid lines
                plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
                plt.grid(which = 'major')
                plt.legend(fontsize = Font)
                plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_' + index[:-4] + '.png', format='png', bbox_inches='tight')
            #
            else:        #if the element is not relevant
                if plot_everything == True :        #if you want to plot all the elements (may be desired?)
                #    
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + index[:-4], fontsize=22, wrap=True)     #title
                    plt.plot(x[:int(len(x)/2)], row[:int(len(x)/2)], 'bo--', 
                         markersize = 5, label = 'Repl_1') 
                    plt.plot(x[int(len(x)/2):], row[int(len(x)/2):], 'ro--', 
                         markersize = 5, label = 'Repl_2') 
                    plt.ylabel(y_label, fontsize= Font)              #ylabel
                    plt.xlabel(x_label, fontsize = Font)
                    plt.tick_params(axis='both', labelsize= Font)              #size of axis
                    #plt.yscale('log') 
                    plt.minorticks_on()             #enabling minor grid lines
                    plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
                    plt.grid(which = 'major')
                    plt.legend(fontsize = Font)            
                    plt.savefig(folder_name +'/' +  
                        pre_save_name + '_' + index[:-4] +'.png', format='png', bbox_inches='tight')
                    #To save plot in folder
        
        
            plt.close()             #to clsoe the plot not to consume too much resources
    
    elif Nrepl == 3:             #3 replicates case!
   
    ###Plot
        for index, row in df_cps.iterrows():   #Loop throught all the rows, the isotopes
                    #df_cps.index give the index values, low and high
		   # 4 because of the way the df is created (and hence the excel tabelle)

            #Saving in the folder
            if index in Elem_rel:  #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
                plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                plt.title(pre_title_plt + index[:-4], fontsize=22, wrap=True)           #title
                plt.plot(x[:int(len(x)/3)], row[:int(len(x)/3)], 'bo--', 
                     markersize = 5, label = 'Repl_1') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x[int(len(x)/3): 2* int(len(x)/3)], row[int(len(x)/3) :2* int(len(x)/3)], 'ro--', 
                     markersize = 5, label = 'Repl_2') 
                plt.plot(x[2* int(len(x)/3):], row[2* int(len(x)/3) :], 'go--', 
                     markersize = 5, label = 'Repl_3') 
                plt.ylabel(y_label, fontsize= Font)              #ylabel
                plt.xlabel(x_label, fontsize = Font)
                plt.tick_params(axis='both', labelsize= Font)              #size of axis
                #plt.yscale('log') 
                plt.minorticks_on()             #enabling minor grid lines
                plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
                plt.grid(which = 'major')
                plt.legend(fontsize = Font)
                plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_' + index[:-4] + '.png', format='png', bbox_inches='tight')
            #
            else:        #if the element is not relevant
                if plot_everything == True :        #if you want to plot all the elements (may be desired?)
                #    
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + index[:-4], fontsize=22, wrap=True)     #title
                    plt.plot(x[:int(len(x)/3)], row[:int(len(x)/3)], 'bo--', 
                     markersize = 5, label = 'Repl_1') 
                    #+1 needed since the df contain a row with the column names!
                    plt.plot(x[int(len(x)/3): 2* int(len(x)/3)], row[int(len(x)/3) :2* int(len(x)/3)], 'ro--', 
                     markersize = 5, label = 'Repl_2') 
                    plt.plot(x[2* int(len(x)/3):], row[2* int(len(x)/3) :], 'go--', 
                     markersize = 5, label = 'Repl_3') 
                    plt.ylabel(y_label, fontsize= Font)              #ylabel
                    plt.xlabel(x_label, fontsize = Font)
                    plt.tick_params(axis='both', labelsize= Font)              #size of axis
                    #plt.yscale('log') 
                    plt.minorticks_on()             #enabling minor grid lines
                    plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
                    plt.grid(which = 'major')
                    plt.legend(fontsize = Font)            
                    plt.savefig(folder_name +'/' +  
                        pre_save_name + '_' + index[:-4] +'.png', format='png', bbox_inches='tight')
                    #To save plot in folder

            plt.close()             #to clsoe the plot not to consume too much resources   
            #
    else:            #Wrong number of replicates 
            print('!!!!! Wrong number of replicates Nrepl, nothing done!!!! \n')    
        
        
    ######### 3) Running time displaying ###############
    '''
    The last thing will be to see and display the time needed
    '''
    
    t_run = tr.time() - t_start     #Running time

    print('###############################################')
    print('Plotting running time: ' + str(t_run) + 's')
    print('###############################################')
    
    
    '''
    THe time when the loop was with iloc and moving a variable i, was 2 times faster! 10s vs 20!
    '''
    

 
#%% ########## 1.17) ICPMS plotter 3 bentonites #############
#####################################

def ICPMS_Plotter3 (x, df_cps, x_label, y_label, folder_name = 'Plots', 
                    plot_everything = False, pre_title_plt = "Concentration of ", 
                    pre_save_name = 'Conc', Elem_rel = Isot_rel ):
    '''
    Function that will plots of the data from the ICPMS (cps) vs another variable,
    initially time, for the 3 bentonites. This assume we have 2 replicates, 
    1 series after the other. Stimated running time around 80s.
    
    *Inputs:
        .x: x axis variable in the plot. dictionary, contaning the 3 df series 
        with the keys: Sard, Tur, BK
        .df_cps: dataframes containing the cps. dictionary containing the 3 df
        series, same order as x. Those are
        outputs for the Read_ICPMS_excel function. Note the 1st column must be
        the one
        with the isotopes (like in the excel)
        .x_label: string that will be the x label for the plot (for math stuff, 
                                    use $$. eg: '$\Delta t[h]$')
        .y_label: string that will be the y label for the plot
        .folder_name: string defining the name of the folder to create to store
        the plots
            default value: 'Plots'
        . plot_everything: string defining if you want to plot all the elements
        or only the relevant ones. Default value: False (only plot relevants)
        .pre_title_plt : title of the graph, part that appears before the name 
        of the elements (thats why pre title).
            Detault value: "Concentration of " (note the space after of, so the
                    element is not together with that!)
        . pre_save_name: name of the graph files to save. Default: 'Conc', 
        giving Conc_Mg24.png for ex    
        .Elem_rel: array containing the name of the relevant elemtns, which
        are the elements that will be saved in a specific folder. 
        Default value: (see above in the script)       
                                    
    *Outputs:
        .Plots (saving them) of the x and df_cps data, cps vs x!
    
    
    ### TO DO: ####
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	.Implement error plotting (in an errorbar pyplot)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    .Improve it, not soo god, averaging replicates and plotting average and 
    its error would be better!
    . TO include 3 replicates version!
    '''
    
    
    ############# 1) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder a 
    subfolder
    with the relevant elements, to be given, will be created
    '''
    
    path_bar_pl = os.getcwd() + '/' + folder_name + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)

    #Subfolder with relevant plots:
    path_bar_pl_rel = os.getcwd() + '/' + folder_name + '/' + 'Relevants' + '/' 
        #folder path for the relevant plots
    
    if not os.path.exists(path_bar_pl_rel):
        os.makedirs(path_bar_pl_rel)   
    
    
    ######### 2) plotting ###############
    '''
    This is a loop plot, so beware, will take long if you plot all the elements
    (280) (2-3mins!).
    
    '''
    Bent_color = {'Sard' : (.68,.24,.31), 'Tur' :  '#EEE8AA', 'BK' : 'grey'} 
                #Color for the ploting, color amtch the actual bentonite color
                
    t_start = tr.time()       #[s] start time of the plot execution
    
    ###Plot

    for i in list( range(df_cps['Sard'].shape[0] ) ):       #Loop thorugh all rows
		   # 4 because of the way the df is created (and hence the excel tabelle)
        #
        
        #Saving in the folder
        if df_cps['Sard'].index[i] in Elem_rel:        #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
            plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
            plt.title(pre_title_plt + df_cps['Sard'].index[i][:-4], fontsize=22, wrap=True)           #title
            #PLot bentonite 1, Sard
            plt.plot(x['Sard'][:int(len(x['Sard'])/2)], df_cps['Sard'].loc[df_cps['Sard'].index[i] ][:int(len(x['Sard'])/2)], 'o--', color = Bent_color['Sard'],
                     markersize = 5, label = 'Repl_1 S') 
                    #+1 needed since the df contain a row with the column names!
            plt.plot(x['Sard'][int(len(x['Sard'])/2):], df_cps['Sard'].loc[df_cps['Sard'].index[i] ][int(len(x['Sard'])/2):], 'o--', color = Bent_color['Sard'],
                     markersize = 5, label = 'Repl_2 S') 
            #PLot bentonite 2, T
            plt.plot(x['Tur'][:int(len(x['Tur'])/2)], df_cps['Tur'].loc[df_cps['Sard'].index[i] ][:int(len(x['Tur'])/2 )], 'o--', color = Bent_color['Tur'],
                     markersize = 5, label = 'Repl_1 T') 
                    #+1 needed since the df contain a row with the column names!
            plt.plot(x['Tur'][int(len(x['Tur'])/2):], df_cps['Tur'].loc[df_cps['Sard'].index[i] ][int(len(x['Tur'])/2 ):], 'o--', color = Bent_color['Tur'],
                     markersize = 5, label = 'Repl_2 T') 
            #PLot bentonite 3, BK
            plt.plot(x['BK'][:int(len(x['BK'])/2)], df_cps['BK'].loc[df_cps['Sard'].index[i] ][:int(len(x['BK'])/2 )], 'o--', color = Bent_color['BK'],
                     markersize = 5, label = 'Repl_1 BK') 
                    #+1 needed since the df contain a row with the column names!
            plt.plot(x['BK'][int(len(x['BK'])/2):], df_cps['BK'].loc[df_cps['Sard'].index[i] ][int(len(x['BK'])/2 ):], 'o--', color = Bent_color['BK'],
                     markersize = 5, label = 'Repl_2 BK') 
            plt.ylabel(y_label, fontsize= Font)              #ylabel
            plt.xlabel(x_label, fontsize = Font)
            plt.tick_params(axis='both', labelsize= Font)              #size of axis
            #plt.yscale('log') 
            plt.minorticks_on()             #enabling minor grid lines
            plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
            plt.grid(which = 'major')
            plt.legend(fontsize = Font)
            plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_'  + df_cps['Sard'].index[i][:-4]+ '.png', format='png', bbox_inches='tight')
            #
        else:        #if the element is not relevant
            if plot_everything == True :     #if you want to plot all the elements (may be desired?)
                #    
                plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
                plt.title(pre_title_plt + df_cps['Sard']['Isotopes'][i], fontsize=22, wrap=True)           #title
                #PLot bentonite 1, Sard
                plt.plot(x['Sard'][:int(len(x['Sard'])/2)], df_cps['Sard'].loc[df_cps['Sard'].index[i] ][:int(len(x['Sard'])/2)], 'o--', color = Bent_color['Sard'],
                     markersize = 5, label = 'Repl_1 S') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x['Sard'][int(len(x['Sard'])/2):], df_cps['Sard'].loc[df_cps['Sard'].index[i] ][int(len(x['Sard'])/2):], 'o--', color = Bent_color['Sard'],
                     markersize = 5, label = 'Repl_2 S') 
                    #PLot bentonite 2, T
                plt.plot(x['Tur'][:int(len(x['Tur'])/2)], df_cps['Tur'].loc[df_cps['Sard'].index[i] ][:int(len(x['Tur'])/2 )], 'o--', color = Bent_color['Tur'],
                     markersize = 5, label = 'Repl_1 T') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x['Tur'][int(len(x['Tur'])/2):], df_cps['Tur'].loc[df_cps['Sard'].index[i] ][int(len(x['Tur'])/2 ):], 'o--', color = Bent_color['Tur'],
                     markersize = 5, label = 'Repl_2 T') 
                    #PLot bentonite 3, BK
                plt.plot(x['BK'][:int(len(x['BK'])/2)], df_cps['BK'].loc[df_cps['Sard'].index[i] ][:int(len(x['BK'])/2 )], 'o--', color = Bent_color['BK'],
                     markersize = 5, label = 'Repl_1 BK') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x['BK'][int(len(x['BK'])/2):], df_cps['BK'].loc[df_cps['Sard'].index[i] ][int(len(x['BK'])/2 ):], 'o--', color = Bent_color['BK'],
                     markersize = 5, label = 'Repl_2 BK') 
                plt.ylabel(y_label, fontsize= Font)              #ylabel
                plt.xlabel(x_label, fontsize = Font)
                plt.tick_params(axis='both', labelsize= Font)              #size of axis
                #plt.yscale('log') 
                plt.minorticks_on()             #enabling minor grid lines
                plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
                plt.grid(which = 'major')
                plt.legend(fontsize = Font)
                plt.savefig(folder_name +'/' +  
                        pre_save_name + '_' + df_cps['Sard'].index[i][:-4] +'.png', format='png', bbox_inches='tight')   #To save plot in folder
        
        
        plt.close()             #to clsoe the plot not to consume too much resources
    
    
    ######### 3) Running time displaying ###############
    '''
    The last thing will be to see and display the time needed
    '''
    
    t_run = tr.time() - t_start     #Running time

    print('###############################################')
    print('Plotting running time: ' + str(t_run) + 's')
    print('###############################################')

    
    
#%%######### 1.18) ICPMS plotter blank appart #############
#####################################

def ICPMS_Plotter_blk (x, df_cps, x_label, y_label, folder_name = 'Plots', 
            plot_everything = False, pre_title_plt = "Concentration of ", 
            pre_save_name = 'Conc', Nrepl = 2, Blank_here = False, 
            Elem_rel = Isot_rel, Logs = 0 ):
    '''
    Function that will plots of the data from the ICPMS (cps) vs another variable,
    initially time, the cps and the rstd. This assume we have 2 replicates, 
    1 series after the other. Blank is plotted sepparately, so the data must 
    include a blank, which should be 1st, the number 1
    
    *Inputs:
        .x: x axis variable in the plot. This could be a pd.Series or 
        pd.DataFrame
        .df_cps: dataframes containing the cps. Those are
        outputs for the Read_ICPMS_excel function. Note the isotopes are the 
        index, so 1st column is 1_1!
        .x_label: string that will be the x label for the plot (for math stuff, 
                                    use $$. eg: '$\Delta t[h]$')
        .y_label: string that will be the y label for the plot
        .folder_name: string defining the name of the folder to create to 
        store the plots default value: 'Plots'
        . plot_everything: string defining if you want to plot all the elements
        or only the relevant ones. Default value: False (only plot relevants)
        .pre_title_plt : title of the graph, part that appears before the name 
        of the elements (thats why pre title). Detault value: 
            "Concentration of " (note the space after of, so the element is 
            not together with that!)
        . pre_save_name: name of the graph files to save. Default: 'Conc', 
        giving Conc_Mg24.png for ex    
        .Nrepl : number of replicates. Default value : 2. 3 value also accepted
        .Elem_rel: array containing the name of the relevant elemtns, which 
        are the elements that will be saved in a specific folder. 
        Default value: (see above in the script)      
         .Blank_here: True if the df contain the blank. Default: False
        .Logs: value defining if applying log scale to x axis, y axis or both:
            0: no log scale
            1: log scale on x axis
            2: log scale on y axis
            3: log scale on both axis
    *Outputs:
        .Plots (saving them) of the x and df_cps data, cps vs x!
    
    Note as since I included the option to have or not the blank, this 
    makes the old funciton ICPMS_Plotter unnecesary. But I
    will not remove it, since does not make any harm there :)
    
    ### TO DO: ####
	.Implement error plotting (in an errorbar pyplot)
    .Implement variable for modifying tittle (eg adding bentonite name, 
                                             as a input)
    '''
    
    
    ############# 1) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder 
    a subfolder with the relevant elements, to be given, will be created
    '''
    
    path_bar_pl = os.getcwd() + '/' + folder_name + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)

    #Subfolder with relevant plots:
    path_bar_pl_rel = os.getcwd() + '/' + folder_name + '/' + 'Relevants' + '/' 
        #folder path for the relevant plots
    
    if not os.path.exists(path_bar_pl_rel):
        os.makedirs(path_bar_pl_rel)   
    
    
    ############ 2) Data cleaning. pre-process #############
    '''
    Here I will clean everything a bit, getting the replicates and so. 
    This depends whether the x data is a df.Series or a df.DataFrame
    '''
    if isinstance(x, pd.Series) :             #x is a pd.Series. y (df_cps) is always pd.DataFrame
        if Nrepl ==2:           #2 replc    
            x_1 = x[:int(len(x)/2) ]
            x_2 = x[int(len(x)/2): ]
            y_1 = df_cps.iloc[ :, : round( df_cps.shape[1] / 2 ) ] 
            y_2 = df_cps.iloc[ :, round( df_cps.shape[1]  / 2 ): ] 
        elif Nrepl ==3 :        #3repl
            x_1 = x[:int(len(x)/3) ]
            x_2 = x[int(len(x)/3): int(2*len(x)/3) ]
            x_3 = x[int(2*len(x)/3):]
            y_1 = df_cps.iloc[ :, : round( df_cps.shape[1] / 3 ) ] 
            y_2 = df_cps.iloc[ :, round( df_cps.shape[1] / 3 ): round( 2*df_cps.shape[1] / 3 )  ] 
            y_3 = df_cps.iloc[ :, round( 2*df_cps.shape[1] / 3 ):  ] 
    elif isinstance(x, pd.DataFrame):
        if Nrepl == 2:
            x_1 = x.iloc[ :, : round( x.shape[1] / 2 ) ] 
            x_2 = x.iloc[ :, round( x.shape[1]  / 2 ): ] 
            y_1 = df_cps.iloc[ :, : round( df_cps.shape[1] / 2 ) ] 
            y_2 = df_cps.iloc[ :, round( df_cps.shape[1]  / 2 ): ] 
        elif Nrepl == 3:             
            x_1 = x.iloc[ :, : round( x.shape[1] / 3 ) ] 
            x_2 = x.iloc[ :, round( x.shape[1] / 3 ): round( 2*x.shape[1] / 3 )  ] 
            x_3 = x.iloc[ :, round( 2*x.shape[1] / 3 ):  ] 
            y_1 = df_cps.iloc[ :, : round( df_cps.shape[1] / 3 ) ] 
            y_2 = df_cps.iloc[ :, round( df_cps.shape[1] / 3 ): round( 2*df_cps.shape[1] / 3 )  ] 
            y_3 = df_cps.iloc[ :, round( 2*df_cps.shape[1] / 3 ):  ] 
            
    '''
    Note that I could plot that more easyly. Still I would need to iterate the df, with an index, to adequately plot it
    but the essence is in the way how I did it, but instead of :
        df_cps.loc[df_cps.index[i] ][int(len(x)/2)+1 : 2*int(len(x)/3)]
        y_2.index[i] should suffice, making everything muuch more clear.
        
    Do it puto!!!!!
    
    
    '''
    ######### 3) plotting ###############
    '''
    This is a loop plot, so beware, will take long if you plot all the elements
    (280) (2-3mins!).
    I inlcude in if statement the numer of replicates, currently only 2 and 3!
    
    I plot only the relevant elements (they are in a list, or all if plot all 
                                       included)
    
    Now it is generalized to be able to be used with pd series and pd df, so 
    another if statement needed!

    '''
    t_start = tr.time()       #[s] start time of the plot execution
      
    if isinstance(x, pd.Series):         #If x is a pd series (like time and so)  
        if Nrepl == 2:               #2 replicates, standard case (x is time for ex)    
            for i in list( range(df_cps.shape[0] ) ):       #Loop thorugh all rows (elements)
                if y_1.index[i] in Elem_rel or plot_everything == True:      #if the element is relevant
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + y_1.index[i][:-4], fontsize=22, wrap=True)     #title
                    if Blank_here:      #if Blank here ==> 1st colum is blk
                        plt.plot(x_1[1:], y_1.loc[y_1.index[i] ][1:], 'bo--', 
                                 markersize = 5, label = 'Repl_1')          #repl 1
                        plt.plot(x_2[1:], y_2.loc[y_2.index[i] ][1:], 'ro--', 
                                 markersize = 5, label = 'Repl_2')          #repl 2                                             
                        plt.hlines(y_1.loc[y_1.index[i] ][0], min(x_1), max(x_1) , label = 'Blk_1', color = 'b' )
                        plt.hlines(y_2.loc[y_2.index[i] ][0], min(x_2), max(x_2) , label = 'Blk_2' , color = 'r' )                                                                           
                    else:                   #No blank!
                        plt.plot(x_1, y_1.loc[y_1.index[i] ], 'bo--', 
                                 markersize = 5, label = 'Repl_1')          #repl 1
                        plt.plot(x_2, y_2.loc[y_2.index[i] ], 'ro--', 
                                 markersize = 5, label = 'Repl_2')          #repl 2
                    plt.ylabel(y_label, fontsize= Font)              #ylabel
                    plt.xlabel(x_label, fontsize = Font)
                    plt.tick_params(axis='both', labelsize= Font)              #size of axis
                    if Logs == 3:                               #x and y in log scale
                        plt.xscale('log')
                        plt.yscale('log')
                    elif Logs == 2:                             #yscale in log
                        plt.yscale('log') 
                    elif Logs == 1:                             #xscale in log
                        plt.xscale('log') 
                        #
                    plt.minorticks_on()             #enabling minor grid lines
                    plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
                    plt.grid(which = 'major')
                    plt.legend(fontsize = Font)
                    plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_'  + df_cps.index[i][:-4] + '.png', format='png', bbox_inches='tight')      
        elif Nrepl ==3:                     #3 replicates
            for i in list( range(df_cps.shape[0] ) ):       #Loop thorugh all rows (elements)
                if y_1.index[i] in Elem_rel or plot_everything == True:      #if the element is relevant
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + y_1.index[i][:-4], fontsize=22, wrap=True)     #title
                    if Blank_here:      #if Blank here ==> 1st colum is blk
                        plt.plot(x_1[1:], y_1.loc[y_1.index[i] ][1:], 'bo--', 
                                 markersize = 5, label = 'Repl_1')          #repl 1
                        plt.plot(x_2[1:], y_2.loc[y_2.index[i] ][1:], 'ro--', 
                                 markersize = 5, label = 'Repl_2')          #repl 2
                        plt.plot(x_3[1:], y_3.loc[y_2.index[i] ][1:], 'go--', 
                                 markersize = 5, label = 'Repl_3')          #repl 3                                               
                        plt.hlines(y_1.loc[y_1.index[i] ][0], min(x_1), max(x_1) , label = 'Blk_1', color = 'b' )
                        plt.hlines(y_2.loc[y_2.index[i] ][0], min(x_2), max(x_2) , label = 'Blk_2', color = 'r' )
                        plt.hlines(y_3.loc[y_3.index[i] ][0], min(x_3), max(x_3) , label = 'Blk_3', color = 'g' )                                                                             
                    else:                   #No blank!
                        plt.plot(x_1, y_1.loc[y_1.index[i] ], 'bo--', 
                                 markersize = 5, label = 'Repl_1')          #repl 1
                        plt.plot(x_2, y_2.loc[y_2.index[i] ], 'ro--', 
                                 markersize = 5, label = 'Repl_2')          #repl 2
                        plt.plot(x_3, y_3.loc[y_2.index[i] ], 'go--', 
                                 markersize = 5, label = 'Repl_3')          #repl 3 
                    plt.ylabel(y_label, fontsize= Font)              #ylabel
                    plt.xlabel(x_label, fontsize = Font)
                    plt.tick_params(axis='both', labelsize= Font)              #size of axis
                    if Logs == 3:                               #x and y in log scale
                        plt.xscale('log')
                        plt.yscale('log')
                    elif Logs == 2:                             #yscale in log
                        plt.yscale('log') 
                    elif Logs == 1:                             #xscale in log
                        plt.xscale('log') 
                        #
                    plt.minorticks_on()             #enabling minor grid lines
                    plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
                    plt.grid(which = 'major')
                    plt.legend(fontsize = Font)
                    plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_'  + df_cps.index[i][:-4] + '.png', format='png', bbox_inches='tight')          
    #
    elif isinstance(x, pd.DataFrame):           #if x is a DataFrame    
        if Nrepl == 2:               #2 replicates, standard case (x is time for ex)    
            for i in list( range(df_cps.shape[0] ) ):       #Loop thorugh all rows (elements)
                if y_1.index[i] in Elem_rel or plot_everything == True:      #if the element is relevant
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + y_1.index[i][:-4], fontsize=22, wrap=True)     #title
                    if Blank_here:      #if Blank here ==> 1st colum is blk
                        plt.plot(x_1.loc[x_1.index[i]][1:], y_1.loc[y_1.index[i] ][1:], 'bo--', 
                                 markersize = 5, label = 'Repl_1')          #repl 1
                        plt.plot(x_2.loc[x_2.index[i]][1:], y_2.loc[y_2.index[i] ][1:], 'ro--', 
                                 markersize = 5, label = 'Repl_2')          #repl 2                                           
                        plt.hlines(y_1.loc[y_1.index[i] ][0], min(x_1.loc[x_1.index[i]]), 
                                   max(x_1.loc[x_1.index[i]]) , label = 'Blk_1', color = 'b' ) 
                        plt.hlines(y_2.loc[y_2.index[i] ][0], min(x_2.loc[x_2.index[i]]), 
                                   max(x_2.loc[x_2.index[i]]) , label = 'Blk_2', color = 'r' )                                                                            
                    else:                   #No blank!
                        plt.plot(x_1.loc[x_1.index[i]], y_1.loc[y_1.index[i] ], 'bo--', 
                                 markersize = 5, label = 'Repl_1')          #repl 1
                        plt.plot(x_2.loc[x_2.index[i]], y_2.loc[y_2.index[i] ], 'ro--', 
                                 markersize = 5, label = 'Repl_2')          #repl 2
                    plt.ylabel(y_label, fontsize= Font)              #ylabel
                    plt.xlabel(x_label, fontsize = Font)
                    plt.tick_params(axis='both', labelsize= Font)              #size of axis
                    if Logs == 3:                               #x and y in log scale
                        plt.xscale('log')
                        plt.yscale('log')
                    elif Logs == 2:                             #yscale in log
                        plt.yscale('log') 
                    elif Logs == 1:                             #xscale in log
                        plt.xscale('log') 
                        #
                    plt.minorticks_on()             #enabling minor grid lines
                    plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
                    plt.grid(which = 'major')
                    plt.legend(fontsize = Font)
                    plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_'  + df_cps.index[i][:-4] + '.png', format='png', bbox_inches='tight')            
        if Nrepl == 3:               #2 replicates, standard case (x is time for ex)    
            for i in list( range(df_cps.shape[0] ) ):       #Loop thorugh all rows (elements)
                if y_1.index[i] in Elem_rel or plot_everything == True:      #if the element is relevant
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + y_1.index[i][:-4], fontsize=22, wrap=True)     #title
                    if Blank_here:      #if Blank here ==> 1st colum is blk
                        plt.plot(x_1.loc[x_1.index[i]][1:], y_1.loc[y_1.index[i] ][1:], 'bo--', 
                                 markersize = 5, label = 'Repl_1')          #repl 1
                        plt.plot(x_2.loc[x_2.index[i]][1:], y_2.loc[y_2.index[i] ][1:], 'ro--', 
                                 markersize = 5, label = 'Repl_2')          #repl 2                                           
                        plt.plot(x_3.loc[x_3.index[i]][1:], y_3.loc[y_3.index[i] ][1:], 'go--', 
                                 markersize = 5, label = 'Repl_3')          #repl 3       
                        plt.hlines(y_1.loc[y_1.index[i] ][0], min(x_1.loc[x_1.index[i]]), 
                                   max(x_1.loc[x_1.index[i]]) , label = 'Blk_1', color = 'b' )
                        plt.hlines(y_2.loc[y_2.index[i] ][0], min(x_2.loc[x_2.index[i]]), 
                                   max(x_2.loc[x_2.index[i]]) , label = 'Blk_2', color = 'r' )                                                                           
                        plt.hlines(y_3.loc[y_3.index[i] ][0], min(x_3.loc[x_3.index[i]]), 
                                   max(x_3.loc[x_3.index[i]]) , label = 'Blk_3', color = 'g' )                                                                          

                    else:                   #No blank!
                        plt.plot(x_1.loc[x_1.index[i]], y_1.loc[y_1.index[i] ], 'bo--', 
                                 markersize = 5, label = 'Repl_1')          #repl 1
                        plt.plot(x_2.loc[x_2.index[i]], y_2.loc[y_2.index[i] ], 'ro--', 
                                 markersize = 5, label = 'Repl_2')          #repl 2                                           
                        plt.plot(x_3.loc[x_3.index[i]], y_3.loc[y_3.index[i] ], 'go--', 
                                 markersize = 5, label = 'Repl_3')          #repl 3       
                    plt.ylabel(y_label, fontsize= Font)              #ylabel
                    plt.xlabel(x_label, fontsize = Font)
                    plt.tick_params(axis='both', labelsize= Font)              #size of axis
                    if Logs == 3:                               #x and y in log scale
                        plt.xscale('log')
                        plt.yscale('log')
                    elif Logs == 2:                             #yscale in log
                        plt.yscale('log') 
                    elif Logs == 1:                             #xscale in log
                        plt.xscale('log') 
                        #
                    plt.minorticks_on()             #enabling minor grid lines
                    plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
                    plt.grid(which = 'major')
                    #Mods to compare, comment!!
                    #plt.xlim(0,3500)
                    #plt.ylim(0, 700000)
                    #####
                    plt.legend(fontsize = Font)
                    plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_'  + df_cps.index[i][:-4] + '.png', format='png', bbox_inches='tight')           
        
        else:
            print('\n Wrong number of replicates given, nothing done xD')
    else:
        print('\n What is x? Not a pd.Series nor pd.DataFrame, so nothing done :) ')
    
    
    ######### 3) Running time displaying ###############
    '''
    The last thing will be to see and display the time needed
    '''
    
    t_run = tr.time() - t_start     #Running time

    print('###############################################')
    print('Plotting running time: ' + str(t_run) + 's')
    print('###############################################')
    
    
    
    
#%%########## 1.19) ICPMS plotter blank appart Average of replicates #############
#####################################

def ICPMS_Plotter_mean_blk (x, std_x, df_mean_cps, df_std_cps, 
                           x_label, y_label, folder_name = 'Plots', Blank_here = False, 
                           plot_everything = False,  LogScale = False, pre_title_plt = "Concentration of ", 
                           pre_save_name = 'Conc', Elem_rel = Isot_rel, Ben_type = 'BK' ):
    '''
    Function that will plots of the data from the ICPMS (cps) vs another variable, initially
    time, the cps and the rstd. This plots the avg value, with its std, so no replicates here.
    In those average values, blank could or not be there. If yes, the blk is plotted as an hor line
    

    *Inputs:
        .x: x axis variable in the plot (mean values). This could be a df series or a df
        .std_x: df Series or df with the std of the x variable
        .df_mean_cps: dataframes containing the cps, but average values (1,2,3,4, etc). From run of the mean and std calc
        .df_std_cps: df containing the std of the mean values of the cps.
        .x_label: string that will be the x label for the plot (for math stuff, 
                                    use $$. eg: '$\Delta t[h]$')
        .y_label: string that will be the y label for the plot
        .folder_name: string defining the name of the folder to create to store the plots
            default value: 'Plots'
        .Blank_here: True if the df contain the blank. Default: False
        . plot_everything: string defining if you want to plot all the elements or only the
            relevant ones. Default value: False (only plot relevants)
        .LogScale: if True, plot x and y variable (cps/ppb/Qe,etc) in logscale
        .pre_title_plt : title of the graph, part that appears before the name of the elements (thats why pre title).
                Detault value: "Concentration of " (note the space after of, so the element is not together with that!)
        . pre_save_name: name of the graph files to save. Default: 'Conc', giving Conc_Mg24.png for ex    
        .Elem_rel: array containing the name of the relevant elemtns, which are the elements that will be saved
            in a specific folder. Default value: (see above in the script)  
        .Ben_type: string with the bentonite type: 'S', 'T' or 'BK'. Default:'BK'. This defines color in plot
                                    
    *Outputs:
        .Plots (saving them) of the x and df_mean_cps data, cps vs x!
    
    Note that for df vs df plotting, Logscale with make both axis in log scale(Ad isoth)
    
    
    ### TO DO: ####
	.Plot 2 blk lines, <>+- std?
    '''
    
    ############# 1) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder a subfolder
    with the relevant elements, to be given, will be created
    '''
    
    path_bar_pl = os.getcwd() + '/' + folder_name + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)

    #Subfolder with relevant plots:
    path_bar_pl_rel = os.getcwd() + '/' + folder_name + '/' + 'Relevants' + '/' 
        #folder path for the relevant plots
    
    if not os.path.exists(path_bar_pl_rel):
        os.makedirs(path_bar_pl_rel)   
    
    
    ######### 2) plotting ###############
    #Colro of the plot:
    if Ben_type == 'BK':
        Color = Bent_color['BK']
    elif Ben_type == 'S':
        Color = Bent_color['Sard']
    elif Ben_type == 'T':
        Color = Bent_color['Tur']
    else:
        print('Wrong bentonite name, color will be black, like your sould')
        Color = 'b'
      
    '''
    This is a loop plot, so beware, will take long if you plot all the elements (280) (2-3mins!).
    I inlcude in if statement the numer of replicates, currently only 2 and 3!
    I need to do a for loop with an index, since I have several df here!

    '''
    t_start = tr.time()       #[s] start time of the plot execution
    
    #TO make it more computing efficient, the Blank if will be here. I will generalize the function, so x
    #can be a df, not only a df.series. I need an additional if, that will be here:
        
    if isinstance(x, pd.Series):                 #If x is a pd Series          
        for i in list( range(df_mean_cps.shape[0] ) ):       #Loop thorugh all rows (elements)
            if df_mean_cps.index[i] in Elem_rel:                      #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
                plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                plt.title(pre_title_plt + df_mean_cps.index[i][:-4], fontsize=22, wrap=True)     #title
                if Blank_here:                                                      ####IS BLANK HERE?
                    plt.errorbar(x[1:], df_mean_cps.loc[df_mean_cps.index[i] ][1:], df_std_cps.loc[df_mean_cps.index[i] ][1:],
                                     std_x[1:], 'o--', color = Color, markersize = 5, label = '<Samples>') 
                #[1:] not to plot sample 1, the blank, which will be a horizontal line!
                    plt.hlines(df_mean_cps.loc[df_mean_cps.index[i] ][0], min(x), max(x), color = Color, label = '<Blk>' )
                else:               #No blank
                    plt.errorbar(x, df_mean_cps.loc[df_mean_cps.index[i] ], df_std_cps.loc[df_mean_cps.index[i] ],
                                     std_x, 'o--', color = Color, markersize = 5, label = '<Samples>') 
                    #Like that you can plot the blank
                plt.ylabel(y_label, fontsize= Font)              #ylabel
                plt.xlabel(x_label, fontsize = Font)
                plt.tick_params(axis='both', labelsize= Font)              #size of axis
                if LogScale:                                                     ####LOG SCALE?
                    plt.yscale('log') 
                    plt.xscale('log') 
                plt.minorticks_on()             #enabling minor grid lines
                plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
                plt.grid(which = 'major')
                plt.legend(fontsize = Font)
                plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_'  + df_mean_cps.index[i][:-4] + '.png', format='png', bbox_inches='tight')       
            else:                                                       #if the element is not relevant
                if plot_everything == True :     #if you want to plot all the elements (may be desired?)
                #    
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + df_mean_cps.index[i][:-4], fontsize=22, wrap=True)     #title
                    if Blank_here:                                                      ####IS BLANK HERE?
                        plt.errorbar(x[1:], df_mean_cps.loc[df_mean_cps.index[i] ][1:], df_std_cps.loc[df_mean_cps.index[i] ][1:],
                                     std_x[1:], 'o--', color = Color, markersize = 5, label = '<Samples>') 
                #[1:] not to plot sample 1, the blank, which will be a horizontal line!
                        plt.hlines(df_mean_cps.loc[df_mean_cps.index[i] ][0], min(x), max(x), color = Color, label = '<Blk>' )
                    else:               #No blank
                        plt.errorbar(x, df_mean_cps.loc[df_mean_cps.index[i] ], df_std_cps.loc[df_mean_cps.index[i] ],
                                     std_x, 'o--', color = Color, markersize = 5, label = '<Samples>') 
                    #Like that you can plot the blank
                    plt.ylabel(y_label, fontsize= Font)              #ylabel
                    plt.xlabel(x_label, fontsize = Font)
                    plt.tick_params(axis='both', labelsize= Font)              #size of axis
                    if LogScale:                                                  ####LOG SCALE?
                        plt.yscale('log') 
                        plt.xscale('log')                 
                    plt.minorticks_on()             #enabling minor grid lines
                    plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
                    plt.grid(which = 'major')
                    plt.legend(fontsize = Font)            
                    plt.savefig(folder_name +'/' +  
                        pre_save_name + '_'  + df_mean_cps.index[i][:-4] +'.png', format='png', bbox_inches='tight')
                    #To save plot in folder
                plt.close()

    elif isinstance(x, pd.DataFrame):                 #If x is a pd Dataframe       
        for i in list( range(df_mean_cps.shape[0] ) ):       #Loop thorugh all rows (elements)
            if df_mean_cps.index[i] in Elem_rel:                      #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
                plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                plt.title(pre_title_plt + df_mean_cps.index[i][:-4], fontsize=22, wrap=True)     #title
                if Blank_here:                                                      ####IS BLANK HERE?
                    plt.errorbar(x.loc[x.index[i]][1:], df_mean_cps.loc[df_mean_cps.index[i] ][1:], df_std_cps.loc[df_mean_cps.index[i] ][1:],
                                     std_x.loc[std_x.index[i]][1:], 'o--', color = Color, markersize = 5, label = '<Samples>') 
                #[1:] not to plot sample 1, the blank, which will be a horizontal line!
                    plt.hlines(df_mean_cps.loc[df_mean_cps.index[i] ][0], min(x), max(x), color = Color, label = '<Blk>' )
                else:               #No blank
                    plt.errorbar(x.loc[x.index[i]], df_mean_cps.loc[df_mean_cps.index[i] ], df_std_cps.loc[df_mean_cps.index[i] ],
                                     std_x.loc[std_x.index[i]], 'o--', color = Color, markersize = 5, label = '<Samples>') 
                    #Like that you can plot the blank
                plt.ylabel(y_label, fontsize= Font)              #ylabel
                plt.xlabel(x_label, fontsize = Font)
                plt.tick_params(axis='both', labelsize= Font)              #size of axis
                if LogScale:                                                     ####LOG SCALE?
                    plt.yscale('log') 
                    plt.xscale('log') 
                plt.minorticks_on()             #enabling minor grid lines
                plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
                plt.grid(which = 'major')
                plt.legend(fontsize = Font)
                plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_'  + df_mean_cps.index[i][:-4] + '.png', format='png', bbox_inches='tight')       
            else:                                                       #if the element is not relevant
                if plot_everything == True :     #if you want to plot all the elements (may be desired?)
                #    
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + df_mean_cps.index[i][:-4], fontsize=22, wrap=True)     #title
                    if Blank_here:                                                      ####IS BLANK HERE?
                        plt.errorbar(x.loc[x.index[i]][1:], df_mean_cps.loc[df_mean_cps.index[i] ][1:], df_std_cps.loc[df_mean_cps.index[i] ][1:],
                                     std_x.loc[std_x.index[i]][1:], 'o--', color = Color, markersize = 5, label = '<Samples>') 
                #[1:] not to plot sample 1, the blank, which will be a horizontal line!
                        plt.hlines(df_mean_cps.loc[df_mean_cps.index[i] ][0], min(x), max(x), color = Color, label = '<Blk>' )
                    else:               #No blank
                        plt.errorbar(x.loc[x.index[i]], df_mean_cps.loc[df_mean_cps.index[i] ], df_std_cps.loc[df_mean_cps.index[i] ],
                                     std_x.loc[std_x.index[i]], 'o--', color = Color, markersize = 5, label = '<Samples>') 
                    #Like that you can plot the blank
                    plt.ylabel(y_label, fontsize= Font)              #ylabel
                    plt.xlabel(x_label, fontsize = Font)
                    plt.tick_params(axis='both', labelsize= Font)              #size of axis
                    if LogScale:                                                  ####LOG SCALE?
                        plt.yscale('log') 
                        plt.xscale('log') 
                
                    plt.minorticks_on()             #enabling minor grid lines
                    plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
                    plt.grid(which = 'major')
                    plt.legend(fontsize = Font)            
                    plt.savefig(folder_name +'/' +  
                        pre_save_name + '_'  + df_mean_cps.index[i][:-4] +'.png', format='png', bbox_inches='tight')
                    #To save plot in folder
                plt.close()
    else:           #Error case
        print('WTF is x bro? xD')
        
        
    '''
    Well, in the past I did the if statemetns outside the for loop, but it seems to be shorter like that! Unexpected. Anyhow its
    shorter the code, save 100lines, I needed to copy that thing 8 times xD
    
    15/11/23, the if of where x is df or series I put it outisde tohugh!
    '''
        
    ######### 3) Running time displaying ###############
    '''
    The last thing will be to see and display the time needed
    '''
    
    t_run = tr.time() - t_start     #Running time

    print('###############################################')
    print('Plotting running time: ' + str(t_run) + 's')
    print('###############################################')
    



#%%### 1.20) ICPMS plotter blank Average of replicates, 3 bentonites! #############
#####################################

def ICPMS_Plotter_mean_3_blk (x_T, std_x_T, df_mean_cps_T, df_std_cps_T,
        x_BK, std_x_BK, df_mean_cps_BK, df_std_cps_BK,
        x_S, std_x_S, df_mean_cps_S, df_std_cps_S,
        x_label, y_label, folder_name = 'Plots', Logscale = False,
        Blank_here = False, plot_everything = False, 
        pre_title_plt = "Concentration of ", 
        pre_save_name = 'Conc', Elem_rel = Isot_rel ):
    '''
    Function that will plots of the data from the ICPMS (cps) vs another variable, 
    initially time, the cps and the rstd, for the 3 bentonites, plotting the average 
    values ideally (output of average computer). It will also allow to plot the 
    blank as a horizontal line, if desired.
    
    *Inputs:
        .x_S/BK/T: x axis variable in the plot (mean values) for each bentonite. 
        This should be a df series
        .std_x_T/BK/S: df Series with the std of the x variable for each bentonite
        .df_mean_cps_T/BK/S: dataframes containing the cps, but average values 
            (1,2,3,4, etc) for the 3 benotnites
            .From run of the mean and std calc
        .df_std_cps_T/BK/S: df containing the std of the mean values of the cps 
            of the 3 benotnites
        .x_label: string that will be the x label for the plot (for math stuff, 
                                    use $$. eg: '$\Delta t[h]$')
        .y_label: string that will be the y label for the plot
        .folder_name: string defining the name of the folder to create to store 
            the plots default value: 'Plots'
        .Blank_here: True if the df contain the blank, to plot it. Default: False
        . plot_everything: string defining if you want to plot all the elements 
            or only the relevant ones. Default value: False (only plot relevants)
        .pre_title_plt : title of the graph, part that appears before the name 
            of the elements (thats why pre title). Detault value: 
            "Concentration of " (note the space after of, so the element is not
            together with that!)
        . pre_save_name: name of the graph files to save. Default: 'Conc', 
            giving Conc_Mg24.png for ex    
        .Elem_rel: array containing the name of the relevant elemtns, which are
            the elements that will be saved in a specific folder. Default 
            value: (see above in the script)   
        .Logscale: string to say if you want the x and y axis in logscale or not.
            Default: False
                                    
    *Outputs:
        .Plots (saving them) of the x and df_mean_cps data, cps vs x!
    
    
    ### TO DO: ####
	.Plot 2 blk lines, <>+- std? OPtional, since for Qe I do not have blank, but for conce I do!
    '''
    
    
    ############# 1) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder a subfolder
    with the relevant elements, to be given, will be created
    '''
    
    path_bar_pl = os.getcwd() + '/' + folder_name + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)

    #Subfolder with relevant plots:
    path_bar_pl_rel = os.getcwd() + '/' + folder_name + '/' + 'Relevants' + '/' 
        #folder path for the relevant plots
    
    if not os.path.exists(path_bar_pl_rel):
        os.makedirs(path_bar_pl_rel)   
    
    
    ######### 2) plotting ###############
    '''
    This is a loop plot, so beware, will take long if you plot all the elements 
    (280) (2-3mins!).I inlcude in if statement the numer of replicates, currently
    only 2 and 3!I need to do a for loop with an index, since I have several df here!

    28/3/24 I generalize this making x variable (and std) a df df, not only a 
    df.series
    '''
    t_start = tr.time()       #[s] start time of the plot execution
    
    ######if all x and std are pd.series
    if isinstance(x_T, pd.Series):          #if x, std_x are df.series (base case)
            #Note I check only one, but I assume the 3 are the same type
        for i in list( range(df_mean_cps_T.shape[0] ) ): 
                            #Loop thorugh all rows (elements)

            if df_mean_cps_T.index[i] in Elem_rel: #if the element is relevant
            #note the -4 is so that that element contain only name and number, 
            #like Mg26, not Mg26 (MR), in order to check with the list!
                plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
                plt.title(pre_title_plt + df_mean_cps_T.index[i][:-4], fontsize=22,
                          wrap=True)     #title
                if Blank_here:  #If Blank is here, to plot it
                    plt.hlines(df_mean_cps_T.loc[df_mean_cps_T.index[i]][0],
                               min(x_T), max(x_T), color = Bent_color['Tur'],
                               label = 'MS Tur')
                    plt.errorbar(x_T[1:], df_mean_cps_T.loc[df_mean_cps_T.index[i]][1:], 
                             df_std_cps_T.loc[df_mean_cps_T.index[i] ][1:],
                         std_x_T[1:], 'o--', markersize = 5, color = Bent_color['Tur'],
                         label = '<Tur>')    #Tur bentonite
                    plt.hlines(df_mean_cps_BK.loc[df_mean_cps_T.index[i]][0],
                               min(x_BK), max(x_BK), color = Bent_color['BK'],
                               label = 'MS BK')
                    plt.errorbar(x_BK[1:], df_mean_cps_BK.loc[df_mean_cps_BK.index[i]][1:], 
                             df_std_cps_BK.loc[df_mean_cps_BK.index[i] ][1:],
                         std_x_BK[1:], 'ro--', markersize = 5, color = Bent_color['BK'],
                         label = '<BK>')    #BK bentonite
                    plt.hlines(df_mean_cps_S.loc[df_mean_cps_S.index[i]][0],
                               min(x_T), max(x_T), color = Bent_color['Sard'],
                               label = 'MS Sar')
                    plt.errorbar(x_S[1:], df_mean_cps_S.loc[df_mean_cps_S.index[i]][1:], 
                             df_std_cps_S.loc[df_mean_cps_S.index[i] ][1:],
                         std_x_S[1:], 'ro--', markersize = 5, color = Bent_color['Sard'],
                         label = '<Sar>')    #Sar bentonite
                else:               #no blank plot
                    plt.errorbar(x_T, df_mean_cps_T.loc[df_mean_cps_T.index[i]], 
                             df_std_cps_T.loc[df_mean_cps_T.index[i] ],
                         std_x_T, 'o--', markersize = 5, color = Bent_color['Tur'],
                         label = '<Tur>')    #Tur bentonite
                    plt.errorbar(x_BK, df_mean_cps_BK.loc[df_mean_cps_BK.index[i]], 
                             df_std_cps_BK.loc[df_mean_cps_BK.index[i] ],
                         std_x_BK, 'ro--', markersize = 5, color = Bent_color['BK'],
                         label = '<BK>')    #BK bentonite
                    plt.errorbar(x_S, df_mean_cps_S.loc[df_mean_cps_S.index[i]], 
                             df_std_cps_S.loc[df_mean_cps_S.index[i] ],
                         std_x_S, 'ro--', markersize = 5, color = Bent_color['Sard'],
                         label = '<Sar>')    #Sar bentonite
                #[1:] not to plot sample 1, the blank, which will be a horizontal line!
            #
                plt.ylabel(y_label, fontsize= Font)              #ylabel
                plt.xlabel(x_label, fontsize = Font)
                plt.tick_params(axis='both', labelsize= Font)    #size of axis
                if Logscale:            #if True
                    plt.yscale('log') 
                    plt.xscale('log') 
                plt.minorticks_on()             #enabling minor grid lines
                plt.grid(which = 'minor', linestyle=':', linewidth=0.5) 
                        #which both to plot major and minor grid lines
                plt.grid(which = 'major')
                plt.legend(fontsize = Font)
                plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_'  + df_mean_cps_T.index[i][:-4] +
                        '.png', format='png', bbox_inches='tight')
            #
            else:        #if the element is not relevant
                if plot_everything == True :     #if you want to plot all the elements (may be desired?)
                #    
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + df_mean_cps_T.index[i][:-4], fontsize=22, wrap=True)     #title
                    if Blank_here:  #If Blank is here, to plot it
                        plt.hlines(df_mean_cps_T.loc[df_mean_cps_T.index[i]][0],
                                   min(x_T), max(x_T), color = Bent_color['Tur'],
                                   label = 'MS Tur')
                        plt.errorbar(x_T[1:], df_mean_cps_T.loc[df_mean_cps_T.index[i]][1:], 
                                 df_std_cps_T.loc[df_mean_cps_T.index[i] ][1:],
                             std_x_T, 'o--', markersize = 5, color = Bent_color['Tur'],
                             label = '<Tur>')    #Tur bentonite
                        plt.hlines(df_mean_cps_BK.loc[df_mean_cps_T.index[i]][0],
                                   min(x_BK), max(x_BK), color = Bent_color['BK'],
                                   label = 'MS BK')
                        plt.errorbar(x_BK[1:], df_mean_cps_BK.loc[df_mean_cps_BK.index[i]][1:], 
                                 df_std_cps_BK.loc[df_mean_cps_BK.index[i] ][1:],
                             std_x_BK, 'ro--', markersize = 5, color = Bent_color['BK'],
                             label = '<BK>')    #BK bentonite
                        plt.hlines(df_mean_cps_S.loc[df_mean_cps_S.index[i]][0],
                                   min(x_T), max(x_T), color = Bent_color['Sard'],
                                   label = 'MS Sar')
                        plt.errorbar(x_S[1:], df_mean_cps_S.loc[df_mean_cps_S.index[i]][1:], 
                                 df_std_cps_S.loc[df_mean_cps_S.index[i] ][1:],
                             std_x_S, 'ro--', markersize = 5, color = Bent_color['Sard'],
                             label = '<Sar>')    #Sar bentonite
                    else:               #no blank plot
                        plt.errorbar(x_T, df_mean_cps_T.loc[df_mean_cps_T.index[i]], 
                                 df_std_cps_T.loc[df_mean_cps_T.index[i] ],
                             std_x_T, 'o--', markersize = 5, color = Bent_color['Tur'],
                             label = '<Tur>')    #Tur bentonite
                        plt.errorbar(x_BK, df_mean_cps_BK.loc[df_mean_cps_BK.index[i]], 
                                 df_std_cps_BK.loc[df_mean_cps_BK.index[i] ],
                             std_x_BK, 'ro--', markersize = 5, color = Bent_color['BK'],
                             label = '<BK>')    #BK bentonite
                        plt.errorbar(x_S, df_mean_cps_S.loc[df_mean_cps_S.index[i]], 
                                 df_std_cps_S.loc[df_mean_cps_S.index[i] ],
                             std_x_S, 'ro--', markersize = 5, color = Bent_color['Sard'],
                             label = '<Sar>')    #Sar bentoniteplt.ylabel(y_label, fontsize= Font)                #ylabel
                    plt.xlabel(x_label, fontsize = Font)
                    plt.tick_params(axis='both', labelsize= Font)              #size of axis
                    if Logscale:        #if True, do it
                        plt.yscale('log') 
                        plt.xscale('log') 
                    plt.minorticks_on()             #enabling minor grid lines
                    plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        
                    plt.grid(which = 'major')
                    plt.legend(fontsize = Font)            
                    plt.savefig(folder_name +'/' +  
                        pre_save_name + '_'  + df_mean_cps_T.index[i][:-4] +
                        '.png', format='png', bbox_inches='tight')
                    #To save plot in folder
        
            plt.close()             #to clsoe the plot not to consume too much resources
            
    elif isinstance(x_T, pd.DataFrame):       #######if x, std_x are df DF!!
        for i in list( range(df_mean_cps_T.shape[0] ) ): #Loop thorugh all rows (elements)

            if df_mean_cps_T.index[i] in Elem_rel: #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
                plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                plt.title(pre_title_plt + df_mean_cps_T.index[i][:-4], 
                          fontsize=22, wrap=True)     #title
                if Blank_here:      #To plot the blank
                    plt.errorbar(x_T.loc[x_T.index[i]][1:], 
                             df_mean_cps_T.loc[df_mean_cps_T.index[i]][1:], 
                             df_std_cps_T.loc[df_mean_cps_T.index[i]][1:],
                         std_x_T.loc[std_x_T.index[i]][1:], 'o--', markersize = 5, 
                         color = Bent_color['Tur'], label = '<Tur>') 
                    plt.hlines(df_mean_cps_T.loc[df_mean_cps_T.index[i]][0],
                               min(x_T.loc[x_T.index[i]]), max(x_T.loc[x_T.index[i]]),
                               color = Bent_color['Tur'], label = 'MS Tur')
                    plt.errorbar(x_BK.loc[x_BK.index[i]][1:], 
                             df_mean_cps_BK.loc[df_mean_cps_BK.index[i] ][1:], 
                             df_std_cps_BK.loc[df_mean_cps_BK.index[i] ][1:],
                         std_x_BK.loc[std_x_BK.index[i]][1:], 'o--', markersize = 5, 
                         color = Bent_color['BK'], label = '<BK>')    #BK bentonite
                    plt.hlines(df_mean_cps_BK.loc[df_mean_cps_BK.index[i]][0],
                               min(x_BK.loc[x_BK.index[i]]), max(x_BK.loc[x_BK.index[i]])
                               ,color = Bent_color['BK'], label = 'MS BK')
                    plt.errorbar(x_S.loc[x_S.index[i]][1:], 
                             df_mean_cps_S.loc[df_mean_cps_S.index[i] ][1:], 
                             df_std_cps_S.loc[df_mean_cps_S.index[i] ][1:],
                         std_x_S.loc[std_x_S.index[i]][1:], 'o--', 
                         markersize = 5, color = Bent_color['Sard'], label = '<Sar>')
                    plt.hlines(df_mean_cps_S.loc[df_mean_cps_S.index[i]][0],
                               min(x_S.loc[x_S.index[i]]), max(x_S.loc[x_S.index[i]])
                               ,color = Bent_color['Sard'], label = 'MS S')
                else:       #no blank plot
                    plt.errorbar(x_T.loc[x_T.index[i]], 
                             df_mean_cps_T.loc[df_mean_cps_T.index[i] ], 
                             df_std_cps_T.loc[df_mean_cps_T.index[i] ],
                         std_x_T.loc[std_x_T.index[i]], 'o--', markersize = 5, 
                         color = Bent_color['Tur'], label = '<Tur>')    #Tur bentonite
                    plt.errorbar(x_BK.loc[x_BK.index[i]], 
                             df_mean_cps_BK.loc[df_mean_cps_BK.index[i] ], 
                             df_std_cps_BK.loc[df_mean_cps_BK.index[i] ],
                         std_x_BK.loc[std_x_BK.index[i]], 'o--', markersize = 5, 
                         color = Bent_color['BK'], label = '<BK>')    #BK bentonite
                    plt.errorbar(x_S.loc[x_S.index[i]], 
                             df_mean_cps_S.loc[df_mean_cps_S.index[i] ], 
                             df_std_cps_S.loc[df_mean_cps_S.index[i] ],
                         std_x_S.loc[std_x_S.index[i]], 'o--', 
                         markersize = 5, color = Bent_color['Sard'], label = '<Sar>')    #Sar bentonite
            #
                plt.ylabel(y_label, fontsize= Font)              #ylabel
                plt.xlabel(x_label, fontsize = Font)
                plt.tick_params(axis='both', labelsize= Font)   #size of axis
                if Logscale:            #if True
                    plt.yscale('log') 
                    plt.xscale('log') 
                plt.minorticks_on()             #enabling minor grid lines
                plt.grid(which = 'minor', linestyle=':', linewidth=0.5) #which both to plot major and minor grid lines
                plt.grid(which = 'major')
                plt.legend(fontsize = Font)
                plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_'  + df_mean_cps_T.index[i][:-4] + '.png', format='png', bbox_inches='tight')
            #
            else:        #if the element is not relevant
                if plot_everything == True :     #if you want to plot all the elements (may be desired?)
                #    
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + df_mean_cps_T.index[i][:-4], 
                              fontsize=22, wrap=True)     #title
                    plt.errorbar(x_T.loc[x_T.index[i]], 
                                 df_mean_cps_T.loc[df_mean_cps_T.index[i] ], 
                                 df_std_cps_T.loc[df_mean_cps_T.index[i] ],
                         std_x_T.loc[std_x_T.index[i]], 'o--', markersize = 5, 
                         color = Bent_color['Tur'], label = '<Tur>')    #Tur bentonite
                    plt.errorbar(x_BK.loc[x_BK.index[i]], 
                                 df_mean_cps_BK.loc[df_mean_cps_BK.index[i] ], 
                                 df_std_cps_BK.loc[df_mean_cps_BK.index[i] ],
                         std_x_BK.loc[std_x_BK.index[i]], 'o--', 
                         markersize = 5, color = Bent_color['BK'], label = '<BK>')    #BK bentonite
                    plt.errorbar(x_S.loc[x_S.index[i]], 
                                 df_mean_cps_S.loc[df_mean_cps_S.index[i] ], 
                                 df_std_cps_S.loc[df_mean_cps_S.index[i] ],
                         std_x_S.loc[std_x_S.index[i]], 'o--', markersize = 5,
                         color = Bent_color['Sard'], label = '<Sar>')    #Sar bentonite
                    #
                    plt.ylabel(y_label, fontsize= Font)                #ylabel
                    plt.xlabel(x_label, fontsize = Font)
                    plt.tick_params(axis='both', labelsize= Font)  #size of axis
                    if Logscale:        #if True, do it
                        plt.yscale('log') 
                        plt.xscale('log') 
                    plt.minorticks_on()             #enabling minor grid lines
                    plt.grid(which = 'minor', linestyle=':', linewidth=0.5) #which both to plot major and minor grid lines
                    plt.grid(which = 'major')
                    plt.legend(fontsize = Font)            
                    plt.savefig(folder_name +'/' +  
                        pre_save_name + '_'  + df_mean_cps_T.index[i][:-4] +
                        '.png', format='png', bbox_inches='tight')
                    #To save plot in folder
        
            plt.close()     #to clsoe the plot not to consume too much resources
            
        

    else:           #errro case
        print('Wrong type of x, std (T, BK, S), what did u put bro? xD')
        
        
    ######### 3) Running time displaying ###############
    '''
    The last thing will be to see and display the time needed
    '''
    
    t_run = tr.time() - t_start     #Running time

    print('###############################################')
    print('Plotting running time: ' + str(t_run) + 's')
    print('###############################################')

    
        

    
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@ End ICPMS shit xD @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@q    
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
 
####################################################
#%% ######### 2.1) PSO fit #############################
###################################################
def PSO_fit(t, Q, delta_t=0, delta_Q =0, folder_name = 'Fits', x_label = 'x', 
            y_label = 'y', Color = 'b', save_name = '', post_title = ' '):    
    '''
    Function to do and compute some variables relevant to the PSO (pseudo second 
    order) kinetic model. THis model comes from
    d(Q(t))7dt = K * (Q_e - Q(t))**2,
    
    where Q_e = Q(t ==> \infty) , the equilibirum sorbed quantity. THe solution of 
    that can be casted in linear form:
        t/Q(t) = 1/KQ_e**2 + t/Q_e
        
    so, plotting t/Q(t) vs t is linear (y = ax + b), being 1/Q_e the slope, 
    and 1/KQ_e**2 the intercept. 
    For the fit I will use my fit function.
    
    You may need to select certain time intervals, and not all of them. Note 
    the units are defined by t and t/Qt!
    
    Note I add + (LR) to the column name in the fit serie!! Watch out, maybe 
    you need to modify it in the future??????
    
    *Inputs
        .t, Q: df series containing the time and Q(t) data. Expected the averaged 
        values. They Must have same index as the df columns in order to plot them!!
        .delta_t/Q: df with the errors of t and Q. Default value = 0, since I do
        not use them!
        .x_label, y_label= x and y label, for the plot. Default value: 'x' and 'y'
        .post_title = '' : title to add after 'Linear fit '
        .save_name = filename of the fit plot, if it wants to be save. Default 
        value = '' ==> no saving.this variable is followed by .png for savinf
        .Color = 'b': color for the plot
        .Folder_name: folder name, where to store the fit plots
    
    
    *Outputs
        .df series with all the relevant info, from the fit and computed 
        quantities, errors (quadratic propagation) included The input units 
        define those units!! Remember saltpepper!
    
    
    '''    
    ############# 0) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder 
    a subfolder with the relevant elements, to be given, will be created
    '''
    
    path_bar_pl = os.getcwd() + '/' + folder_name + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)


    ############## 1) Calcs #################
    #I need to compute t/Q(t) to do the PSO fit!
    
    t__Q = t / Q          #t/Q(t) for S
    Delta_t__Q = np.abs(t__Q) * np.sqrt((delta_Q / Q )**2 + (delta_t /t )**2 )  
                #error, unused!!
    
    ############# 2)Fit ######################
    
    fit = Fits.LinearRegression(t, t__Q, delta_t, Delta_t__Q,
                                   x_label = x_label, y_label = y_label, 
                                   Color = Color, 
                save_name = folder_name +'/' + save_name, post_title = post_title)       
                            #Fit (i dont use npo variable, fit variable)
    
    ################ 3) Model parameters ################
    '''
    From that I can also get Qe and K easy:
        y = ax + b;
         a = 1/Qe ==> Qe = 1/a
         b = 1/KQe**2 == > K = 1/bQe**2 = a**2 /b
    '''
    fit['Q_e'] = 1 / fit['a']         #Qe = 1/a, y= ax + b
    fit['\Delta(Q_e)'] = fit['\Delta(a)'] /fit['a']**2     #Delta(Qe)
    fit['K'] = fit['a']**2 /fit['b']         #K = 1/b * Qe**2 = a**2/b
    fit['\Delta(K)'] = np.abs(fit['K']) * np.sqrt( 2*(fit['\Delta(a)'] / fit['a'] )**2 + (fit['\Delta(b)'] / fit['b'])**2 )  
                            #Delta(K) np.abs() so its always >0

    ########## 4) Return ###########
    
    return fit
    

#%% ######### 2.2) PFO fit #############################
###################################################
def PFO_fit(t, Q, delta_t=0, delta_Q =0, p_0 = None, folder_name = 'Fits', x_label = 'x',
            y_label = 'y', Color = 'b', save_name = '', post_title = ' ', npo=100):   
    '''
    Function to do and compute some variables relevant to the PFO (pseudo first order) kinetic model. THis model
    comes from
    d(Q(t))/dt = K * (Q_e - Q(t))
    
    where Q_e = Q(t ==> \infty) , the equilibirum sorbed quantity. THe solution of that is:
        Q(t) = Q_e* (1- exp(-t*K_1) ) 
        
    I could fit the data to that equation. Note that usually that is casted into linear form:
            ln (Qe - Q) = ln Qe - K_1t,
    And the Qe-Q you compute by using the experimetnal Qe, from the graph, and from the fit
    you get the other. Bro, WFT??????
    
    So I do the fit. You may need to select certain time intervals, and not all of them. 
    Note the units are defined by t and Q!
    
    Note I add + (LR) to the column name in the fit serie!! Watch out, maybe you need to modify it in the future??????
    
    *Inputs
        .t, Q: df series containing the time and Q(t) data. Expected the averaged values
            . Must have same index as the df columns in order to plot them!!
        .delta_t/Q: df with the errors of t and Q. Default value = 0, since I do not use them!
        .p_0 = None: initial stimation of the fit parameters. 
        .x_label, y_label= x and y label, for the plot. Default value: 'x' and 'y'
        .post_title = '' : title to add after 'Linear fit '
        .save_name = filename of the fit plot, if it wants to be save. Default value = '' ==> no saving.
                    this variable is followed by .png for savinf
        .Color = 'b': color for the plot
        .Folder_name: folder name, where to store the fit plots
        .npo=100: number of points for the fit plot
    
    
    *Outputs
        .df series with all the relevant info, from the fit and computed quantities, errors (quadratic propagation) included
            The input units define those units!! Remember saltpepper!
    
    
    '''    
    ############# 0.1) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder a subfolder
    with the relevant elements, to be given, will be created
    '''
    
    path_bar_pl = os.getcwd() + '/' + folder_name + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)

    
    ######## 0.2) Fit eq #########
    def PFO_eq(x, C_1, C_2):            #equation to fit!
        return C_1 * (1 - np.exp(-C_2 * x))


    ############# 1)Fit ######################
    param, cov = curve_fit(PFO_eq, t, Q, p0 = p_0)        #easy, in principle! No initial stim of parameters given!
    std = np.sqrt(np.diag(cov))     #std, since a cov amtrix is returned (read in the raw code)
    
    ################ 2) Model parameters ################
    Qe = param[0]
    K1 = param[1]
    Delta_Qe = std[0]
    Delta_K1 = std[1]
    
    #Storing them in a df Series
    values = {'Qe' : Qe, '\Delta(Qe)' : Delta_Qe,
              'K1' : K1, '\Delta(K1)' : Delta_K1}
    Ser_values = pd.Series(values, name = post_title)      #gathering output in a df Series
            #naming the column like the post_title variable, since this variable is an isotope: U238    
    
    
    ############# 3) Plot of the fit##########
    t_vector = np.linspace(min(t),max(t),npo)         #for the fit plotting
    
    fig = plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
    ax = fig.add_subplot(111)
    ax.errorbar(t, Q, delta_Q, delta_t, 'o', color = Color, markersize = 5, label = 'Data')
    ax.plot(t_vector, PFO_eq(t_vector, Qe, K1),'--', color = Color,
            label= 'Fit: ' + y_label + f' = {Qe:.1e} ' + '$\cdot$ [1- exp(-'+ x_label + '$\cdot$' +f'+{K1:.1e} )]')      #fit
            #.2f to show 2 decimals on the coefficients!
            #2e for scientific notation with 2 significative digits
    ax.set_title('PFO fit ' + post_title, fontsize=22)          #title
    ax.set_xlabel(x_label, fontsize= Font)                                    #xlabel
    ax.set_ylabel(y_label, fontsize= Font)                                    #ylabel
    ax.tick_params(axis='both', labelsize= Font)            #size of tick labels  
    ax.grid(True)                                              #show grid
    ax.legend(fontsize = Font)             #legend
                    #Plot of the fit equation. (0,0) is lower-left corner, and (1,1) the upper right
    plt.savefig(folder_name +'/' + save_name +'.png', format='png', bbox_inches='tight')                
                    ###This require some thoughts!!!!! to automatize the show of the equation!!!!!!!!!!!
    
    
    ########## 4) Return ###########
    
    return Ser_values


#%% ######### 2.3) Freundlich isot fit #############################
###################################################
def Fre_fit(Ce, Qe, delta_Ce=0, delta_Qe =0, p_0 = None, folder_name = 'Fits',
            x_label = 'log($C_e [ng/g]$)', y_label = 'log($Q_e [ng/g_{be}]$)',
            Color = 'b', save_name = '', post_title = ' ', npo=100):   
    '''
    Function to do and compute some variables relevant to the Freundlich fit
    of an adsorption isotherm Q_e = f (C_e), the sorbed quantity as a function
    of the equilibrium concentration. Its equation is
    
        Q_e = K_F * C_e**n,           K_F, n constants
    That can be linearized into (10 log is fine)
        loq Q_e = log K_F + n log C_e
    
    The lineal fit of that is trivial:  y= ax + b
                y= log Q_e
                x= Log C_e
                a = n
                b= log K_F
    (easier than if I use 1/n as constant!)
    
    Note I add + (LR) to the column name in the fit serie!! Watch out, maybe you need 
    to modify it in the future??????
    
    *Inputs
        .Ce, Qe: df series containing C_e and Q_e data. Expected the averaged values
            . Must have same index as the df columns in order to plot them!!
        .delta_Ce, delta_Qe: df with their uncertainties. Default value = 0, since I do not use 
        them!
        .p_0 = None: initial stimation of the fit parameters. 
        .x_label, y_label= x and y label, for the plot. Default value: 'log($C_e [ng/g]$)' and 
                    log($Q_e [ng/g_{be}]$)' respectively!
        .post_title = '' : title to add after 'Linear fit '
        .save_name = filename of the fit plot, if it wants to be save. Default value = '' ==> no saving.
                    this variable is followed by .png for savinf
        .Color = 'b': color for the plot
        .Folder_name: folder name, where to store the fit plots
        .npo=100: number of points for the fit plot
    
    
    *Outputs
        .df series with all the relevant info, from the fit and computed quantities, 
        errors (quadratic propagation) included
            The input units define those units!! Remember saltpepper!
    
    
    '''    
    ############# 0.1) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder a subfolder
    with the relevant elements, to be given, will be created
    '''
    
    path_bar_pl = os.getcwd() + '/' + folder_name + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)

    
    ############## 1) Calcs #################
    #I need to compute the logarithms!
    logCe = np.log10(Ce)
    logQe = np.log10(Qe)
    
    delta_logCe = delta_Ce / np.abs(Ce)         #error of the log!
    delta_logQe = delta_Qe / np.abs(Qe)

    ############# 2)Fit ######################
    
    fit = Fits.LinearRegression(logCe, logQe, delta_logCe, delta_logQe,
                                   x_label = x_label, y_label = y_label, 
                                   x_legend = 'log($C_e$)', y_legend = 'log($Q_e$)',
                                   Color = Color, 
                                   save_name = folder_name +'/' + save_name, 
                                   post_title = post_title)       
                            #Fit (i dont use npo variable, fit variable)
                #note that for the legnd I delete the units!!
    
    
    ################ 3) Model parameters ################
    '''
    From that I can also get the constants; applying log10 == log, no ln!!! :
        y= ax + b : loq Qe = n log Ce + log KF
                    y= log Q_e
                    x= Log C_e
                    a = n
                    b= log K_F ==> K_F = 10**b
                    delta_K_F = delta_b * ln(10) * K_F
                    delta_n = delta_a
                    (partial log(Kf) / partial Kf = 1/ (KF ln(10)) )
    
    And the units?
    n is adimensional, ofc
    K_F not, since Q_e = K_F * C_e**1/n ==> K_F = Q_e * C_e*n ==>
    K_F is ng/g_be * (ng/g_tot)**n = ng**n+1/(g_be * g_tot**n)
    '''
    fit['n'] = fit['a']         
    fit['\Delta(n)'] = fit['\Delta(a)'] 
    fit['K_F[ng^(1-n)*g_tot^n/(g_be)]'] = 10**fit['b']        
    fit['\Delta(K_F[ng^(1-n)*g_tot^n/(g_be)])'] = fit['K_F[ng^(1-n)*g_tot^n/(g_be)]'] *fit['\Delta(b)'] * np.log(10)
    
    
    
    # ############# 4) Plot of the fit##########
    # logCe_vector = np.linspace( min(logCe),max(logCe),npo )         #for the fit plotting
    
    # fig = plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
    # ax = fig.add_subplot(111)
    # ax.errorbar(logCe, logQe, delta_logCe, delta_logQe, 'o', color = Color, markersize = 5,
    #             label = 'Data')
    # ax.plot(logCe_vector, logCe_vector*fit['a'] + *fit['b'],'--', color = Color,
    #         label= 'Fit: ' + y_label + f' = {Qe:.1e} ' + '$\cdot$ [1- exp(-'+ x_label + '$\cdot$' +f'+{K1:.1e} )]')      #fit
    #         #.2f to show 2 decimals on the coefficients!
    #         #2e for scientific notation with 2 significative digits
    # ax.set_title('Freundlich fit ' + post_title, fontsize=22)          #title
    # ax.set_xlabel(x_label, fontsize= Font)                                    #xlabel
    # ax.set_ylabel(y_label, fontsize= Font)                                    #ylabel
    # ax.tick_params(axis='both', labelsize= Font)            #size of tick labels  
    # ax.grid(True)                                              #show grid
    # ax.legend(fontsize = Font)             #legend
    #                 #Plot of the fit equation. (0,0) is lower-left corner, and (1,1) the upper right
    # plt.savefig(folder_name +'/' + save_name +'.png', format='png', bbox_inches='tight')                
    #                 ###This require some thoughts!!!!! to automatize the show of the equation!!!!!!!!!!!
    
    
    ########## 4) Return ###########
    
    return fit



#%% ######### 2.4) Langmuir isot fit #############################
###################################################

def Lang_fit(Ce, Qe, delta_Ce=0, delta_Qe =0, p_0 = None, Fit_type = 1,
             folder_name = 'Fits', x_label = '$C_e [ng/g]$', 
             y_label = '$C_e/Q_e [g_{be}/g_{tot}]$',
            Color = 'b', save_name = '', post_title = ' ', npo=100):   
    '''
    Function to do and compute some variables relevant to the 
    Langmuir fit of an adsorption isotherm Q_e = f (C_e), the sorbed 
    quantity as a function of the equilibrium concentration. 
    Its equation is
        Q_e = Q_max * K_L/(1+K_L* C_e) * C_e
        
    Being K_L the sorption rate. Note it has a maximum!
    
    That function is usually linearized in order to fit it. Two common 
    linearization are:
        1) C_e/Q_e = C_e/Q_max + 1/(Q_max *K_L)
        The best linearization is (Guo2019, in terms, of minimizing 
                                   uncertainties):
        2) Q_e/C_e = -K_L* Q_e + Q_max * K_L
    
    I will do these function so that you can do one of the 2 linearizations
    
    *Inputs
        .Ce, Qe: df series containing C_e and Q_e data. Expected the averaged
                values. Must have same index as the df columns in order to 
                plot them!!
        .delta_Ce, delta_Qe: df with their uncertainties. Default value = 0, 
        since I do not use them!
        .p_0 = None: initial stimation of the fit parameters. 
        .x_label, y_label= x and y label, for the plot. Default value: 
            'log($C_e [ng/g]$)' and log($Q_e [ng/g_{be}]$)' respectively!
        .post_title = '' : title to add after 'Linear fit '
        .save_name = filename of the fit plot, if it wants to be save. 
            Default value = '' ==> no saving. this variable is followed 
            by .png for savinf
        .Color = 'b': color for the plot
        .Folder_name: folder name, where to store the fit plots
        .npo=100: number of points for the fit plot
        .Fit_type: linearization type. Default: 1, meaning linearization 1
    
    
    *Outputs
        .df series with all the relevant info, from the fit and 
                computed quantities, errors (quadratic propagation) included
            The input units define those units!! Remember saltpepper!
    
    
    '''    
    ############# 0.1) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder
    a subfolder with the relevant elements, to be given, will be created
    '''
    
    path_bar_pl = os.getcwd() + '/' + folder_name + '/'
        #Note os.getcwd() give current directory. With that structure we 
        #are able to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)

    
    ############## 1) Calcs #################
    Ce_Qe = Ce / Qe
    delta_Ce_Qe = Ce_Qe * np.sqrt((delta_Ce / Ce)**2 + 
                                  (delta_Qe / Qe)**2 )
    Qe_Ce = Qe / Ce
    delta_Qe_Ce = Qe_Ce * np.sqrt((delta_Ce / Ce)**2 + 
                                  (delta_Qe / Qe)**2 )

    ############# 2) Fit and parameters ######################
    
    if Fit_type ==1:         #Linearizaiton 1
        fit = Fits.LinearRegression(Ce, Ce_Qe, delta_Ce, delta_Ce_Qe,
                                   x_label = x_label, y_label = y_label, 
                                   x_legend = '$C_e$', 
                                   y_legend = '$C_e/Q_e$',
                                   Color = Color, 
                                   save_name = folder_name +'/' + save_name, 
                                   post_title = post_title)       
                            #Fit (i dont use npo variable, fit variable)
                #note that for the legnd I delete the units!!
            #
        print("############################/n")
        print("Linearization 1 done, Ce/Qe vs Ce! /n")
        print("############################/n")
                ####### Parameters obtention #############
        """
        C_e/Q_e = C_e/Q_max + 1/(Q_max *K_L)   
        y = ax+b
     
        then:
            a = 1/Qmax ==> Qmax = 1/a
            b = 1/(Qmax*K_L) ==> K_L = 1/(b*Qmax)
     
        and their uncertainty:
            delta Q_max / Qmax = delta a/a ==> delta Qmax = delta a /a**2
            delta KL /KL = sqrt( (delta b/b)**2 + (delta Qmax/Qmax)**2 ) ==>
            delta KL = Kl * sqrt( (delta b/b)**2 + (delta Qmax/Qmax)**2 )
         
        """          
        fit["Q_max[ng/g_be]"] = 1/ fit["a"]
        fit["\Delta(Q_max[ng/g_be])"] = fit['\Delta(a)']/fit['a']**2
        fit["K_L[g/ng]"] = 1/ (fit["b"] * fit["Q_max[ng/g_be]"] )
        fit["/Delta(K_L[g/ng])"] = fit["K_L[g/ng]"] * np.sqrt(
             (fit["\Delta(Q_max[ng/g_be])"]/fit["Q_max[ng/g_be]"])**2 + 
             (fit['\Delta(b)']/fit["b"])**2 )
         
    else:       #Linearization 2!
        fit = Fits.LinearRegression(Qe, Qe_Ce, delta_Qe, delta_Qe_Ce,
                                   x_label = x_label, y_label = y_label, 
                                   x_legend = '$Q_e$', 
                                   y_legend = '$Q_e/C_e$',
                                   Color = Color, 
                                   save_name = folder_name +'/' + save_name, 
                                   post_title = post_title)       
                            #Fit (i dont use npo variable, fit variable)
                #note that for the legnd I delete the units!!
            #
        print("############################/n")
        print("Linearization 2 done, Qe/Ce vs Qe! /n")
        print("############################/n")
    
                ####### Parameters obtention #############
        """
        Q_e/C_e = -K_L* Q_e + Q_max * K_L 
        y = ax+b
     
        then:
            a = - K_L
            b = Qmax*K_L ==> Qmax = b/K_L
     
        and their uncertainty:
            delta K_L = delta a
            delta Qm /Qm = sqrt( (delta b/b)**2 + (delta KL/KL)**2 ) ==>
            delta Qm = Qm * sqrt( (delta b/b)**2 + (delta KL/KL)**2 )
         
     """          
        fit["K_L[g/ng]"] = -fit["a"]
        fit["Q_max[ng/g_be]"] = fit["b"] / fit["K_L[g/ng]"]
        
        fit["/Delta(K_L[g/ng])"] = fit['\Delta(a)']        
        fit["\Delta(Q_max[ng/g_be])"] = fit["Q_max[ng/g_be]"] * np.sqrt(
             (fit["/Delta(K_L[g/ng])"]/fit["K_L[g/ng]"])**2 + 
             (fit['\Delta(b)']/fit["b"])**2 )

    
    # ############# 4) Plot of the fit##########
    # logCe_vector = np.linspace( min(logCe),max(logCe),npo )         #for the fit plotting
    
    # fig = plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
    # ax = fig.add_subplot(111)
    # ax.errorbar(logCe, logQe, delta_logCe, delta_logQe, 'o', color = Color, markersize = 5,
    #             label = 'Data')
    # ax.plot(logCe_vector, logCe_vector*fit['a'] + *fit['b'],'--', color = Color,
    #         label= 'Fit: ' + y_label + f' = {Qe:.1e} ' + '$\cdot$ [1- exp(-'+ x_label + '$\cdot$' +f'+{K1:.1e} )]')      #fit
    #         #.2f to show 2 decimals on the coefficients!
    #         #2e for scientific notation with 2 significative digits
    # ax.set_title('Freundlich fit ' + post_title, fontsize=22)          #title
    # ax.set_xlabel(x_label, fontsize= Font)                                    #xlabel
    # ax.set_ylabel(y_label, fontsize= Font)                                    #ylabel
    # ax.tick_params(axis='both', labelsize= Font)            #size of tick labels  
    # ax.grid(True)                                              #show grid
    # ax.legend(fontsize = Font)             #legend
    #                 #Plot of the fit equation. (0,0) is lower-left corner, and (1,1) the upper right
    # plt.savefig(folder_name +'/' + save_name +'.png', format='png', bbox_inches='tight')                
    #                 ###This require some thoughts!!!!! to automatize the show of the equation!!!!!!!!!!!
    
    
    ########## 5) Return ###########
    
    return fit



    
#%% ######### 4) TGA reader ##################### 
##################################################


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
    
    
 
#%% ######### 5.1) XRD reader F141 (Cold lab) #################### 
######################################################
 
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


#%% ######### 5.2) XRD reader, F130 (hot XRD) #################### 
######################################################

def Read_XRD_F130 (name, t_paso = 10, Skip_rows = 266):
    '''
    Function that reads the .ras file from the XRD in F130 (Olaf Walter), returning a df
    with the relevant info. It will also plot it and save it, in the same folder!
    
    The .ras file contain lot of text at the beginning (1st 266 lines), then 3 variables:
            2Theta, Counts, Unknown
        The unknown seem to be always 1, so I delete it.
        
    Some of the text I skip contian relevant info. Eg:
            *MEAS_SCAN_START_TIME "04/10/24 15:21:09" (near the last lines of text)
    
    *Inputs:
        .name: name of the file. Ej: 'file.ras'
        .t_paso: time per step, in seconds. Ej: 10 [s]. Thi assumes the operation 
        mode in step, the usual
        .Skip_rows: number of rows of the .ras file to skip. Default: 266. Spotted 
        from opening the file
        
    *Output
        .df with the 2Theta, Counts and cps
    
    '''
    
    ##### 1. Reading #########
    #Reading was not trivial, but I came up with the way of doing it:
    aux = pd.read_csv(name, skiprows = Skip_rows, sep = ' ', 
                      names = ['2Theta[°]', 'Counts', 'npi'], encoding='latin')
    
    '''That worked, but 2 last columns are text, that I can delete easily. Also I have 3 column, 
    being the 3rd always 1, so I will also delte it:'''
    
    df = aux.iloc[:-2,:-1]        #Removing 3rd column, and last 2 rows
    
    #Since the 2theta data is not numeric because of these 2 rows, I need to make them numeric:
    df = df.apply(pd.to_numeric)  
    
    #The cps are easy to get, since I define the time measuring each step:
        
    df['CPS'] = df['Counts']/t_paso         #Computing Counts Per Secons (CPS)


    '''
    I will also return the cps nromalized, min value 0 and max 1. This can be accomplished
    by appling the operation:
            (x-min(x)) /(max(x)-min(x))         para todo x
    '''
    
    df['CPS_norm'] = ( df['CPS'] -df['CPS'].min() ) / (df['CPS'].max() - df['CPS'].min())
                        #normalized cps, min value 0, and max 1
    
    
    ##Prints useful for command line sytling (also for interlaminar space)
    print('\n##############')
    print('Currently working with:\n') 
    print( name[:-4])
    
    ########### 2. Plotting ##
    plt.figure(figsize=(11, 8))  # width, heigh 6.4*4.8 inches by default (11,8)
    # I need to enlarge since because of the title is big, I think
    plt.title('XRD: ' + name[:-4], fontsize=22,
              wrap=True, loc='center')  # title
    plt.plot(df['2Theta[°]'], df['CPS_norm'], label='data')
    plt.xlabel(" $2 \Theta [°]$", fontsize=18)  # ylabel
    plt.ylabel('cps', fontsize=18)
    plt.tick_params(axis='both', labelsize=18)  # size of axis
    plt.minorticks_on()             #enabling minor grid lines
    plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
    plt.grid(which = 'major')
    #plt.legend(fontsize=18)
    plt.savefig( name[:-4]+ '.png', format='png',
                bbox_inches='tight')  # To save plot, same name as file, change extensio
    plt.show()
    
    
    #Finally, return the df:
    return df


#%% ######### 5.3) XRD, Get interlaminar space bentonites #################### 
######################################################

def XRD_Get_interl_sp (XRD_df, DosTheta_inter, Kalpha = 1.5401):
    '''
    Function that will get the interlaminar space of the bentonites, based on the first basal
    refelction 001, which will satisfy:
            lambda = 2d sin (theta)    n = 1
            
    For this, we need:
        1. Perform the gaussian fit of the peak (Fits module used!)
        2. Compute the interlaminar space d
    
    The .ras file contain lot of text at the beginning (1st 266 lines), then 3 variables:
            2Theta, Counts, Unknown
        The unknown seem to be always 1, so I delete it.
        
    Some of the text I skip contian relevant info. Eg:
            *MEAS_SCAN_START_TIME "04/10/24 15:21:09" (near the last lines of text)
    
    *Inputs:
        .XRD_df: df containing the info from the XRD. Important the names, must have
            '2Theta[°]' and 'CPS_norm' (I refer to those
        .DosTheta_inter: array with min and maximum values of the 2Theta variable from the XRD df,
        to do the fit. Ej: [5,30]
        .kalpha: Wavelength of the Kalpha1 radiation of the Cu, the element of the XRD = 1.5401Angs
        .name: name of the file. Ej: 'file.ras'

    *Output
        .df with the interlaminar space!
    
    '''
    
    
    ####### 1. Gaussian fit #########
    #We need to get the desired interval to perform the fit, eye spotted from the XRD diagram
    
    inter_x =XRD_df['2Theta[°]'].loc[ (XRD_df['2Theta[°]'] < DosTheta_inter[1]) & 
                                        (XRD_df['2Theta[°]'] > DosTheta_inter[0]) ]
        #getting desired interval
    #Now I need the cps values (y) of those x values. Since the index is preserve, I could use it
    inter_y = XRD_df['CPS_norm'][inter_x.index]
    
    Fit = Fits.Gaussian_fit(inter_x.values, inter_y.values)            #Gaussian fit

    ########## 2. Interlaminas space calc
    #Note that np.sin works in radian, so I convert the angle from degree to radians
    
    d = Kalpha / (2*np.sin(np.deg2rad(Fit['mean'][0] / 2) ) )      #[Ams]
    Delta_d = d * np.cos(np.deg2rad(Fit['mean'][0] / 2) ) * Fit['\Delta(mean)'][0]/2/ np.sin(
                                                np.deg2rad(Fit['mean'][0] / 2) ) #[Ams]

    aux = {'d[A]' : d, '\Delta(d)[A]' : Delta_d }    #variable containing everyting
        
        #Let´s print that so it appears in the command line:
    print('\n#######################\n')
    print('Interlaminar space of  (see above)')
    print(aux)
    print('\n ###############')
        
    
    #Finally, return d and its error in a df, recycling the one for the fit
    Fit['d[A]'] = d
    Fit['\Delta(d)[A]'] = Delta_d
    
    return Fit


