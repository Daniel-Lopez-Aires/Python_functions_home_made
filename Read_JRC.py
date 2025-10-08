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

For printing a break of line, \n goes at the end of the print, not at the 
beginning!
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
import re #regex, to group by
from scipy.stats import pearsonr, spearmanr #statistical correlation tests

# ----------- Useful variables -----------------------------

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
            #General list of nuclides of interest!
            #Reserve: 'Ti46', 'Ti47', 'Ti48', 'Ti49', 'Ti50',
            #Eu151 less abundant as Eu153, but Eu153 sufffer interferences from
            #Ba oxides,so for low Eu concentrations, Eu151 better!! [Stef]   
Elem_rel = ['Si(MR)', 'Al(MR)', 'Mg(MR)', 'Mn(MR)', 'Fe(MR)', 
            'Ca(MR)', 'Na(MR)', 'P(MR)', #bentonite elements
            'Sr(LR)', 'Cs(LR)', 'Eu(LR)','La(LR)', 'U(LR)', 'Mo(LR)'] 
        #General list of elements of interest (merging isotopes)

"""
Isot rele Cs are from the Cs sep
"""
Font = 18               #Fontsize, for the plots (labels, ticks, legends, etc)           
Markersize = 7                  #markersize for all the plots

#For my personal laptop
#At_we_rad = pd.read_excel('/home/dla/Python/at_wt_natural_elements_SVW.xlsx',
                       #   'To_read_atom_weight', index_col=0)   #read of the
                    #excel with the atomic weights
#For guest laptop at JRC:  
# At_we_rad = pd.read_excel('C:/Users/localadmin/Desktop/Python/at_wt_natural_elements_SVW.xlsx', 
#                       index_col=0)
At_we_rad = pd.read_excel('C:/Users/Administrator/Desktop/Python/at_wt_natural_elements_SVW.xlsx', 
                      index_col=0, sheet_name = 'To_read_atom_weight_rad')  
            #containing radioactive info, manually added by me 
At_we_nat = pd.read_excel('C:/Users/Administrator/Desktop/Python/at_wt_natural_elements_SVW.xlsx', 
                      index_col=0, sheet_name = 'To_read_atom_weight')      
        #Path from guest laptop from JRC)       
              
Rad_dat = pd.read_excel('C:/Users/Administrator/Desktop/Python/Rad_dat_DLA.xlsx', 
                      index_col=0)      #half lifes and masses of radionuclides

        
#############################################################            
#%%## ## 1.1) ICPMS excel reader #############
#####################################

def Read_ICPMS_excel (excel_name, Sheet_name = 'To_read', Numbering_row = 1,
                      Is_wash_inside = 0,
                      return_debug = False):
    '''
    Function that will read a excel sheet from the excel file from ICPMS and 
    will return a df with the relevant information, for easier handling /plotting. 
    The sheet could have any format, raw cps from Stef could be read.
    Assumptions for the function:
        .Numbering row right before the hearde row (row with names)
        .4 rows after the header, data starts!
        
    
    You could use this function to get the raw data (output from ICPMS) or to 
    correct them for the ICPMS dilution factor.     

    Note sometimes some random NaN data from excel can be added. Easy solution, 
    go to the excel sheet, and delete those rows/columns. No clue why this 
    happens, but that solves it (:   
                                                                                 
    The wash solution will be removed, if desired. BEware, there should not be 
    any extra info at the bottom, so for example the "Flags" thing that Rafa 
    sometimes includes at the end, move it to the top, or simply delete!!

    
    *Inputs:
        .excel_name: string with the name of the excel, with the .xlsx. note if 
        you select the file on the FIle viewer and do copy paste, the name will 
        be there. But dont forget the '', it must be an string! Eg: 'Excel.xlsx'
        .Sheet_name: string with the name of the sheet with the data to read 
            future maybe also concentration values?). Default value: 'To_read' 
        (from acid vs no acid test)
        .Numbering_now: number to indicate the starting row (from excel), which 
        contains the sample number: 1,2,3,4... Default: 1 (in the sheet to_read  
        is the 1st row, in the raw outcome from icpms might not be)
        .Is_wash_inside: boolean to indicate whether the wash is inside or not.
            Default: 0 (is not inside)
        .return_debug: if you want to get some extra df for debug 
        (raw data, without cleaning, so like the excel). Default value = False

        
    *Outputs:
        .A pandas df with the cps/%rsd or whatever it is reading. Depending whether
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
    Dat = pd.read_excel(excel_name, Sheet_name, header = Numbering_row, 
                        index_col=0)
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
    df_cps = Dat.iloc[4:,:] #get data from row 4 on
    
    #Another cleaning will be putting the df in numeric format. It is in object 
    #format, which gives problems
    #(YeroDivisionerror) with divisions, while for numbers there is no problem.

    df_cps = df_cps.apply(pd.to_numeric)    
                
    ###Removing the wash data, if sinde
    if Is_wash_inside:
        df_cps = df_cps.iloc[:,:-1]     #Removing last columns
    
    
    ############## 3) Ouptuts ######################
        
    if return_debug == True:    #return the debug
        #
        return df_cps, Dat              #return of values
    
    else: #no return debug, but yes corrections
            return df_cps

    '''
    Sucessfully debbugged, enjoy bro! Note I return the last the Df corrected 
    because  I may not be interested in that, for ex if I want to plot the 
    raw data.
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
    .m_s, V_s: array with the mass (g) and volumes (L) of the samples, or single 
        values. When it is an array, the column names must be the sample names,
        as in the df_ppb. Beware for the MS case, I was not doing that!
        Default values:
        m_S = 1000g, V_s = 1L, to have desntiy 1kg/L
    .Delta_m_s, Delta_V_s: values of the uncertainties. Defaults:
            Delta_m_s = 0.0001g, Delta_V_s = 0.001L
            
    *Outputs
    .df with the concentration in M
    .df with the uncertainties in M, ASSUMING no uncertainty in the atomic weights
    '''

        
            
    ############### 1) Data cleaning #####
    '''
    Some data, specially the masses and volumes may not be numeric, so we must
    convert them to numeric, to perform the operations. However, this can only
    be made if the variables are not a float (number). I need to check for that
    using if statements.
    
    Also, if we are speaking about series, some data of m_s, V_s may be missing.
    This is particularly True for the MS case, since the CL of BIC I did not
    record any mass, since all the solution I put was either BIC, or the CL.
    Then, I will replace 0 and NaN, by the average value.
    '''
    if isinstance(m_s, float) == False and isinstance(m_s, int) == False: 
                    #Case m_s is not a float/integer, a series
        m_s = m_s.apply(pd.to_numeric) 
        #
        #To replace 0 and NaN by the average, 1st I need to replcae 0 to NaN,
        #since NaN do not affect in the average, and then do it:
        m_s.replace(0, np.nan, inplace = True)
        m_s.replace(np.nan, m_s.mean(), inplace = True)
    
    if isinstance(V_s, float) == False and isinstance(V_s, int) == False:
                        #Case V_s is not a float
        V_s = V_s.apply(pd.to_numeric) 
        #
        V_s.replace(0, np.nan, inplace = True)
        V_s.replace(np.nan, V_s.mean(), inplace = True)
    
    
    ################ 2) Conversion from ppb to M #########3
    '''
    If I am reading a ppb data, I could convert it to Mol, I would only need
        1) Atomic weights
        2) Density of the samples, which I could have, since I measured both
    
    The operation is easy:
        ng/gtot *gtot/Vtot *mol/g = M ==> ppb * rho * At = M
    '''
                    #
    
    rho = m_s/V_s               #density of the samples
    Delta_rho = rho*np.sqrt((Delta_m_s/m_s)**2 + (Delta_V_s/V_s)**2)
                #uncertainty of the density
                
    #Creation of empty df to store the data
    df_M = pd.DataFrame( index =df_ppb.index, columns = df_ppb.columns )
                    #empty df, but with defined columns and rows
    
    print('------------- Converting from ppb to M...-------------------------')
    print('Beware with the atomic weight df, if unproperly read will give problems! ')
    print('--------------------------------\n')
    
    #Now to apply it I would need a loop
    for isotope in df_ppb.index:
        elem = isotope[:-4]  # remove isotope number suffix (MR) or (LR)
        if elem not in At_we_rad.index:
            print(f" ! Skipping {isotope}: not found in At_we_rad")
            continue  # skip this isotope
        
        At_mass = At_we_rad.loc[elem, "atomic mass (u)"]
        df_M.loc[isotope, :] = (
            df_ppb.loc[isotope, :] * 1e-9 / At_mass * rho)
    
    #It is not numeric, so lets convert it:
    df_M = df_M.apply(pd.to_numeric)            #conversion to numeric
    
    #The uncertainty calc can be outside the loop
    df_M_std = df_M * np.sqrt((Delta_rho/rho)**2 + (df_ppb_std/df_ppb)**2)
            #ASSUMING no error in the atomic weights!!!!!!!      


    ############# 2) Output ##########

    Output= {'dat' : df_M, 'dat_std': df_M_std, 'dat_%rsd': df_M_std/df_M*100}
    
    return Output             
                    
                    
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
            1) Dlution actor for the ICPMS sample preparation (simeq 50). In this
                        case you add 
                                    .8.8mL HNO3 1M
                                    .1mL IS (2IS; 0.5mL each)
                                    .0.2mL sample
            2) Dilution factor for the sample you use for the ICPMS sample prep 
                        (simeq 1). 
                In this case you add some HNO3 conc (65% w/w)to the sample, 
                to stabilize it.
                                                                            
    Note the excel should be a bit preprocessed:
            
        1) Including sample preparation (neccesary to get the Dilution factor). 
            In this excel, the sample names must be the same as the names in the top
               of the ICPMS data (do it manually, lazy spaniard, less siesta and
            more work!). This sheet must be called 'Sample_prep'. Note the 
            structure is like ICPMS results, column for sample. This is needed
            for further operations (corrections and so)
               
    *Inputs:
        .excel_name: string with the name of the excel, with the .xlsx. note if 
        you select the file on the FIle viewer and do copy paste, the name will
        be there. But dont forget the '', it must be an string! Eg: 'Excel.xlsx'
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

    D_f = Dat_sa_prep.iloc[D_f_data[1]-1, 1: len(Dat_sa_prep)]   
                #Dilution factor (pandas Series)
    #The -1 is because python start in 0 while excel in 1 for counting rows and columns
    
    
    #I need to put the correct index names for the operations, that can be done like:
    D_f.index = Dat_sa_prep.iloc[D_f_data[0]-1, 1: len(Dat_sa_prep)]     
                            #proper index name (to operate)
    
    
    ######## 3) cleaning ###############
    '''
    If there is no dil 2, in the dilution 2 table, you might encounter
    Div#0, so we will remove it, replcaing it by 1:
    '''
    D_f = D_f.replace('#DIV/0!',1)
    
    #Also, if there are NaN rows and columns (empty in the excel), I will remove:
    D_f= D_f.loc[D_f.index.dropna()]
    
    '''
    Since its object type, I will make it numeric, since everything will be easier
    with it (and strictly its true)
    '''
    D_f = D_f.apply(pd.to_numeric)              #conversion to numeric data type
    
    ########### 4) Return #############
    '''
    Finally we return Df, which is a pandas serie
    '''
    
    return D_f            #return


#%%######## 1.4) ICPMS Sens finder #############
#####################################

def ICPMS_Sens_finder(Excel_name, Excel_sheet = 'Calib', Sens_column = 12):
    '''
    Function to get the senstivity of each mass from the excel. It is an easy
    task to do, you just read in thee xcel where you did the calibration.
    
    *Inputs:
        .Excel_name
        .Excel_sheet
        .Sens_colum: number indicating the column where the Sens is.Spotted from
        excel, so column A is number 1
    
    *Outputs:
        .df with the sens, its std and rstd [cps/ppb]
    '''
    
    
    ########## Loading #################
    Sheet = pd.read_excel(Excel_name, Excel_sheet, header = 5, index_col =0)

    i = 12 #excel count, of Sens column!
    #from that the relevant data is:
    df_sens = Sheet[Sheet.columns[i-2:i+1]]

    #Now I need to just rename then and there they are!
    df_sens.columns = ['Sens [cps/ppb]', 'std [cps/ppb]', 'rsd [cps/ppb]']
    
    
    return df_sens

#%% ##### 1.4) ICPMS Dilution factor corrector #############
#####################################

def ICPMS_Df_corrector (df_dat, Df):
    '''
    Function that will apply the correction for the dilution factor to the ICPMS results.
    Note There is 2 dilution factors involved:
            1) Dlution factor for the ICPMS sample preparation (simeq 50). In this
                    case you add 
                                    .8.8mL HNO3 1M
                                    .1mL IS (2IS; 0.5mL each)
                                    .0.2mL sample
            2) Dilution factor for the sample you use for the ICPMS sample prep (simeq 1). 
                    In this case you add some HNO3 conc (65% w/w)to the sample, 
                    to stabilize it.
                                                                            
                                                                            
    The correction is essentially scalatin for that factor, so the results takes 
    into account that only a portion was measuring. So:
            df_dat * Df (Df >=1)
    
    Its fundamental the labelling is appropiate! The data_exp sheet AND 
    Sample_prep sheet must have same labels as the ICPMS output!!!!! Once the 
    ICPMS resutls are there, change name to both, otherwise will return NaN!!!!

    *Inputs:
        .df_dat: dataframe containing the cleaned data, to which the correction 
        should be applied.  the isotopes are the index, so all columns are data!
        .D_f: pandas series containing the dilution factor to apply. Note the 
        labelling is crutial, both inputs should have same labels (remember 
        that you change the name in the exp sheet to match the names that Stefaan used)

    *Outputs:
        .df with the correction factor (Df) applied

	CAn the [:,:] be removed ?? I would say yes, but check!!!
        '''
    
    
    ########### 1) Calcs ###########
    '''Now, I can apply the Dilution factor to both the std and the cps, should 
    be straightforward, same fashion than above. Well, not as simple, since for 
    the multiplication the indexes should be the same so, I redefined (above) 
    the Df indexes so they matched the Df ones, and then that calc is straightforward
    '''
    
    df_corrected = pd.DataFrame(df_dat * Df,
                                columns = df_dat.columns, 
                                index = df_dat.index)  #computing the correction
    
    #Just in case that df is not numeric (depends on Df mostly), I will convert
    #it to numeric
    
    df_corrected = df_corrected.apply(pd.to_numeric)
    
    ########### 2) Return #############
    return df_corrected             #return



#%% ######## 1.5) ICPMS Sample blank substraction #############
#####################################

def ICPMS_Sample_Blk_corrector (df_dat, Nrepl = 2):
    '''
    Function that will apply the sample blank correction to the other samples in
    the ICPMS results df. This version is the 2 replicates version. Remember we
    already applied in the excel the IS correction and the ICPMS blanks corrections.
    Now its time for the blanks you did in your experiment.

    The correction is essentially substracting the blank (number 1) to the rest
    of the samples, but involved some operations since we have dataframes. So 
    those will be here. Necessary that the data contain no Div0, ensure in the
    excel by deleting those! that only a portion was measuring. 
    
    Only for 2 and 3 replicates, but could be generalized for N replicates easily.


    *Inputs:
        .df_dat: dataframe containing the data, the full data, with the 2 replicates.
        This is the  output fromt he reader function. You could apply this before
        or after Df corrections. Format: isotopes as index, columns the samples, 
        1st 1st replicate, then 2nd replicate. 2 replicates assume this function!!!!
        .Nrepl: number of replicates. Default value: 2 (2 replicates). Can also be 3

    *Outputs:
        .df with the correction factor (Df) applied
        '''
    
    
    ########### 1) Calcs ###########
    '''
    I must treat the 2 experiments are different, I should substract the blank 
    to the 1st emasurements and the 2 to the others. Since I ordered it in the
    right way (1st replicacte 1, then replicate 2, I could) split it easily :D
            df.shape gives the shape of the df, n_rows, n_columns
    
    Note the df have number of samples * 2 replicates columns.
    
    Then, I will create a new dataframe substracting that data. To do so, I need 
    to get rid of the isotopes column, since is text, and then add it again. 
    Watch, the substraction is easy with a pandas mehotd.

    I shuold then remove those columns
    from there, and replace negatives values for 0, for a good plot
    '''

    if Nrepl ==2:               #2 replicates, standard case
        df_1 = df_dat.iloc[ :, 0: round( ( df_dat.shape[1] ) / 2 ) ]   #1st replicate
        df_2 = df_dat.iloc[ :, round( ( df_dat.shape[1] ) / 2 ) :  ]   #replicate 2
    
        #The next step is blank substraction
        df_1_blk = df_1.subtract(df_1.iloc[:,0], axis = 0 ) 
                #0 since 1st columns is the blank  
        df_1_blk.drop( [df_1.iloc[:,0].name], axis = 1, inplace = True) 
                #removing the value I use to substract
    
    #Finally, for plotting purposes, we will replace negative values with 0:
        df_1_blk[df_1_blk < 0] = 0          #substituyin negative values with 0! 
                                    #This needs that no Div0 in excel!    
    #And now we do the same for replicate 2   
    
        df_2_blk = df_2.subtract(df_2.iloc[:,0],  axis = 0 )        #substraction
        df_2_blk.drop( [df_2.iloc[:,0].name], axis = 1, inplace = True)     
        df_2_blk[df_2_blk < 0] = 0      #replcaing negative values with 0
            #Those 3 lines would be needed for a loop, so sohuld be easy, if needed
    
    #Finally, lets store it

        df_blk = pd.concat( [df_1_blk, df_2_blk], axis = 1)   
                #mergind the 2 little df ina  huge one
    
    elif Nrepl == 3:         #3 replicates case
        #Copy paste the before code, now for 3 replicates
        df_1 = df_dat.iloc[:, : round(df_dat.shape[1] / 3)]
        df_2 = df_dat.iloc[:, round(df_dat.shape[1] / 3): 2*round(df_dat.shape[1] / 3)]
        df_3 = df_dat.iloc[:, 2*round(df_dat.shape[1] / 3) :]

        df_1_blk = df_1.subtract(df_1.iloc[:,0], axis = 0 ) 
        df_2_blk = df_2.subtract(df_2.iloc[:,0], axis = 0 ) 
        df_3_blk = df_3.subtract(df_3.iloc[:,0], axis = 0 ) 

        df_1_blk.drop( [df_1.iloc[:,0].name], axis = 1, inplace = True)
        df_2_blk.drop( [df_2.iloc[:,0].name], axis = 1, inplace = True)
        df_3_blk.drop( [df_3.iloc[:,0].name], axis = 1, inplace = True)
        
        df_1_blk[df_1_blk < 0] = 0          #replcaing negative values with 0
        df_2_blk[df_2_blk < 0] = 0     
        df_3_blk[df_3_blk < 0] = 0      
        
        df_blk = pd.concat( [df_1_blk, df_2_blk, df_3_blk], axis = 1)         
                        #mergind the 2 little df ina  huge one
    
    else:                   #Otherwise, error
        print('Wrong number of replicates Nrepl!, nothing has been done!')
        df_blk = 0                      #Error case, not doing anything
    
    
    ########### 2) Return #############
    return df_blk             #return



#%%######## 1.6) ICPMS: Get Mass number #############
#####################################

def Get_A_Resol(isotope_name):
    '''
    Function to get the mass number A and the resolution type from the isotope 
    name of the ICPMS excel. The name is in the format:
            .BB111(MR)
    where BB can be 1 or 2 leters indicating the chemical symbol, and 111 can 
    be 2 numbers also, the mass number. (MR) or (LR) is low or high resolution. 
    There are 2 excepctions, Ar40Ar40(LR/MR) and U238O16(LR/MR).
    
    * Input
        .isotope_name: string with the ICPMS excel format, eg 'U238(LR)', 
        'Co59(MR)', etc
        
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
        1) Name with 1 or 2 letters? with isalpha() you get True if letter. 
            isdigit() for numbers)
    I need to introduce the 2 expcetions as different cases also. All in a loop, 
    since if i get the exceptions, I do not want to continue
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
        IS_meas = ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)',
                   'Co59(MR)', 'In115(MR)'],
        name_IS_sens_LR_plot = 'IS_sensLR_plot', 
        name_IS_sens_MR_plot = 'IS_sensMR_plot'):
    '''
    Function part of the ICPMS Data processing! Function that will compute the 
    IS sens (cps/ppb) for a df containing the cps and ppb, and another df with
    its std. The format is like Stefaans excel, his first sheet. It IS needed
    the ppb table below the cps data. How below? Do not care, I find the data
    using locate functions ;) Note std is computed using the squared error 
    propagation methods!  

    
    !!!!!!ASSUMPTION!!!! : Delta(ppb IS)/ppb IS = 1%, completely arbitrary, 
    to simplify everything, as a 1st approx! it may be changed after! Those ppb
    come from perkin Elmer multielementa solutions, and derived calcs 
    (logbooks, calculations of concentrations, etc)
        !!!!!!!!!!!!!!!!!!
    
    Note the plots appear in the plots pannel in Spyder, but are also saved ;)

    * Input
        .df_cps_ppb: df containing the cps and ppb data. THe fashion is like 
        Stefaans raw sheet. Isotopes are index. Take care of the names of the
        columns (like in IS conc table or so, the values you find), if they 
        dont exist will give error! For the ppb data, the names are just 
        Co-59, Ho-165, In-115, Th-232.This is extremely important. Otherwise it
        will not work!
        . df_std: df containng the std of the cps. Default value: None, so I
        dont need to give it, since in the past I was not giving it, not to 
        obtian error with the old scripts
        .IS_meas: array containing in a list the measured Internal Standards, 
        containing its resolution, like how they appear in the isotopes column. 
        Default value:
        ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)']
            That mean those isotopes were measured. If Ho165(MR) also measured, 
            just included it, and fine ;)
        .name_IS_sens_LR_plot: name for the plot of the IS sens for LR case. 
        Similar for MR. Default values: 'IS_sensLR_plot' and 'IS_sensMR_plot'. 
        To that is added the .png to create the file.
              
    #Output
        .df with the IS sens data
        .df with the std of the IS sens.
        .Plot with the IS sens plot, 2 plots, 1 for LR and other for MR, in .png file
    
    Note that if this only one argument for output written, then the output 
    will be a tuple with 2 df        
        
    ###### TO Do
        .Function to compute the ppb data, from Sum isobars? May be complex, 
            possibly too time consuming ?
        . Think a way so that std can be optionally given, so my old scripts 
            are still working? The alternative is modifying them, inclduding the
            std, will not be fatal though xD
    '''
    
    # if df_std == 'No':      #if no std passed
    #     df_std = df_cps_ppb * 0 
    #Creating a ceros matrix if no error provided, so we avoid doing an if loop
    
    
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
            
        cps_IS = df_cps_ppb.loc[df_cps_ppb.index == value_to_find]  
                            #cps of the IS (give full row)
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
    the column to index:we need to insert the isotopes column, giving as values 
    the isotopes names or so:
    '''
    #df_IS_sens['Isotopes'] = ['Co LR', 'In LR', 'Ho LR', 'Th LR', 'Co MR', 'In MR']
    df_IS_sens['Isotopes'] = IS_meas        #setting isotopes name from the loop variable, better
    df_IS_sens.set_index('Isotopes', inplace = True)
        
    df_IS_sens_std['Isotopes'] = IS_meas
    df_IS_sens_std.set_index('Isotopes', inplace = True)

    '''
    As a final stylish, I can put the same column names as the df used to create 
    the IS df:
    '''
    df_IS_sens.columns = df_cps_ppb.columns
    df_IS_sens_std.columns = df_cps_ppb.columns
    
    print('############################################ \n')
    print('Attention: Here, to compute the std of the sens, assuming error in '+
          'theoretical ppb data of 1% Take into account, or at least remember !!!')
    print('############################################ \n')
    
    'I can print the %rstd values = std/mean * 100 of the IS sens, useful to spot fluctuations!'
    rstd = df_IS_sens.std(axis =1) / df_IS_sens.mean(axis =1) * 100  
            #%rstd of the IS sens. ofc agrees with excel!
    print('################################################################# \n')
    print('%rstd of the IS sens:')
    print(rstd)
    print('Values >= 10/15% start to be suspicious, something happened! (Th oxidation for ex?) \n')
    print('############################################################### \n')
    
    
    #################################################
    ############# IS sens plotter ####################
    ########
    '''
    I need to do 2 plots, 1 for the LR and other for the MR. I should sort them 
    somehow. I found how xD, explicit plotting, recommended by python, better 
    than implicit(what you normally do xD)
    '''
    
    pltL = plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default  LR plot!
    axL = pltL.subplots()
    axL.set_title("IS sens LR along measuring sequence", fontsize=22, wrap=True) 
            #title
   
    pltM = plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default   MR plot!
    axM = pltM.subplots()
    axM.set_title("IS sens MR along measuring sequence", fontsize=22, wrap=True)           
    
    for i in range(0, df_IS_sens.shape[0]):   #looping through rows of the df
        if df_IS_sens.index[i][-4:] == '(LR)':      #low resolution
            axL.plot(list(range(0, df_IS_sens.shape[1])), df_IS_sens.iloc[i,:],
                     '-o' ,label = df_IS_sens.index[i]) 
            
        else:   #MR
            axM.plot(list(range(0, df_IS_sens.shape[1])), df_IS_sens.iloc[i,:],
                     '-o' ,label = df_IS_sens.index[i])  
    
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
    pltL.savefig(name_IS_sens_LR_plot + '.png', format='png', bbox_inches='tight')    
            #note I call plt, not ax!
    pltL.show()                         #showing the plot
    pltM.savefig(name_IS_sens_MR_plot + '.png', format='png', bbox_inches='tight')
    pltM.show()
         
    ############## Return of values ##########################
    return df_IS_sens, df_IS_sens_std            #mass is an string!


#%% ########## 1.8) ICPMS IS sens correction #############
#####################################

def IS_sens_correction(df_raw, df_raw_std, df_IS_sens, df_IS_sens_std,
        IS_meas = ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)',
                   'Co59(MR)', 'In115(MR)']):
    '''
    PART OF THE ICPMS Data processing!
    
    Function that will apply the IS sens correction to the cps data, from the raw 
    ICPMS data. This is aprt of the ICPMS data analysis. You need to correct 
    for IS sens, then apply ICPMS blanks, and then calibrate with 2+4 and 3+5 IS.
    
    'IS conc [ppb]' should be the box that indicates where the ppb chart start. 
    That is, in A column in excel, that should be written, since that is used to
    find the ppb data!!!
    
    Note the plots appear in the plots pannel in Spyder, and that the IS 
    measured decide the correction, I defined them with if statements, have not
    yet come up with a better idea...

    
    * Input
        .df_raw: df containing the raw cps and ppb data. THe fashion is like Stefaans
                raw sheet
        .df_raw_std: df containing the std from the cps data. 
        .df_IS_sens: df containing the IS sens data (cps/ppb). Ensure the size 
        and indexing are appropiates!
        .df_IS_sens_std: df containing the IS sens data std (cps/ppb)
        .IS_meas: array containing in a list the measured Internal Standards, 
        containing its resolution, like  how they appear in the isotopes column. 
        Default value:
        ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)']
        
    #Output
        .df with the IS corrected data, containing the cps and ppb data
        .df with the std of the IS corrected data
        
        
    ###### TO Do #############
        .Switch cases if other corrections needed (avoid one IS because the 
                                                 measuring was wrong?)
    #########################
    '''
    
    ############## 1) Calcs #########################
    '''
    This is simple. If the 4 IS are fine, I do the correction in the following fashion:
            .From 59 to 80 apply Co59
            .From 81 to 138 apply In115
            .From 139 to 209 Ho165
            .From 210 to 248 Th232

    If there are less IS, then I will make less divisions. Ex, if for MR there is
    only Co and In,
            .Co from 59 to 80
            .Rest In115
            
    The correction is:
        new cps = cps / sens * <sens>
        
    Being sens = cps/ppb teor of the IS used to correct, and <sens> te average
    of all the sensitivies.
        
        
    I shold be able to detech the name of th eisotope (row), get number, and apply
    those limits. The correction is data / IS sens * <IS sens>, being sens = cps/ppb

    The name format is:
        .A11(LR)
        .AB11(LR)
        .A111(LR)
        .AB111(MR)
        .And some random elements, like Ar40Ar40(L/MR), U238O16(L/MR)

    So, avoiding the weird stuff, the chemical name can have 1 or 2 letters, 
    and the number can be 2 or 3 ciphers. I created a function to get the mass 
    number. Note the operations should be line by line. The structure would be 
    loop like, at the beginning I check the mass number, and apply. Lets do a
    single one
    
    Now I can try to do a loop. I could give where the IS sens start from the 
    df inspection (240), and I know the order:
        .Co LR
        .In LR
        .Ho R
        . Th LR
        .CO MR
        . In MR
        
    Next step is to create a function to get the mass and apply one or other 
    correction. We already have the function, so now we need to create a loop 
    (after a function) that get the mass, and apply one or other correction
    '''   
    
    df_IS_co = df_raw.copy()       
        #Thats the proper way to create it, copy so it is not asigned to it!
    df_IS_co_std = df_raw.copy() 
    
    '''
    The loop should go until the alst isotope, which I can find by finding IS 
    conc ppb, the first thing for the ppb chart! so, THIS data is mandatory that
    exist like that!! 
    To generalize, I will do if statement for the different IS measured cases. 
    By far only 2, the sequence done in the past, the 4 in LR and 2 in MR; and 
    now (8/23) 4 in MR also, the most detailed case. 11/23, only LR elements!
    '''
    
    if IS_meas == ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 
                   'Co59(MR)', 'In115(MR)']:   
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
    
    Function to be run after the IS sensitivity correction, to apply the next step
    in the ICPMS data analysis process, the ICPMS blanks correction (std 0ppt and 
                blanks). 
    
    Note that in the excel, int he part of the ppb table, you should delete the 
    names, that stefaan write twice, for this to work!!!! 
    Also needed that the ppb data table is introduced by:  "IS conc [ppb]". beware,
    sometimes Stef writes (ppb), not [ppb]!!
    
    
    *Inputs:
        .columns_blks: np.array([]) indicating the number of the columns contaning 
        blanks (std 0ppt and blank std).
            Ex: columns_blanks = np.array([1, 2]). Those columns number are from excel!
        .df_IS_co: df containing the cps (also ppb, not needed but there it is), 
        output from the IS sens correction funciton. The formatting is like the excel. 
        .df_IS_co_std: df containing the std of the cps corrected from the IS
    *Outputs:
        .df containinng the cps data corrected for the ICPMS blanks. Note it also
        contains the ppb data, but now modified so they are random numbers. 
        .df containing the std of the cps corrected for the ICPMS blanks. 
            Quadratic error propagation used!
            
            
    ##### To DO ###########
        *Improve style and remove ppb data (would need modifications for the IS 
                sens function)
    
    '''
    
    ###################################################################
    '''
    So, for the Blank correction, I need:
        1) To put ina  df all of the blanks. I could indicate column numbers, KISS
        2) Compute average from the cps of the blank
        3) Substract either the average or the ICPMS blanks, or the minimum cps
            value of the samples (from IS corrected, no ICPMS blanks), to ensure
            no negative values. 

    I would really need the IS corrected data without the sens also bro, but I 
    could get it easily I think
    '''
    
    columns_blks = columns_blks - 1 
        #to adapt to the system in the df, isotopes are index, not columns!
    
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
    df_blanks['<cps>'] = df_blanks.mean(axis = 1) 
                #I append the mean values. Gives a warning, 
                                                            #but nothing is wrong!
    df_blanks_std ['<cps>'] = df_blanks.std(axis = 1) 

    #3) The substraction
    '''
    THe first step is getting the minimun values of the cps data
    Those values are for the columns of the samples except the blanks!
    '''
    df_Is_co_no_blk = df_IS_co.drop(df_IS_co.columns[columns_blks ], axis = 1)
                        #df initial without the ICPMS blanks!! 

    '''
    I could do the min to that, but since also contain the cps data there are 
    some NaN, which make things not work. If I do fill nan with the mean values,
    that could fix it. Lets try bro! 

    '''

    df_Is_co_no_blk.fillna(999999.999, inplace = True) 
            #filling NaN values with 99999 (easily recognizable)

    #Now the min values are, storing them in the same df:
    df_Is_co_no_blk['min cps'] = df_Is_co_no_blk.min(axis = 1) 
                #axis = 1 for columns!
    '''
    For this I would need a loop to check for each row which value to substract.
    I would say we always substract the mean, and Stefaan also think that. adn 
    I probe that matematically, so it is like that xD

    So, now the problem is perform that operation here on python, since the df 
    are different.

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
                       name_plot_LR_bef = 'IS_sensLR_plotBEF', 
                       name_plot_MR_bef = 'IS_sensMR_plotBEF',
                       name_plot_LR_aft = 'IS_sensLR_plot', 
                       name_plot_MR_aft = 'IS_sensMR_plot',
  IS_meas = ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)'],
                       excel_name = 'df_IS_Blks_corr.xlsx'):
    '''
    SUITE of ICPMS Data Processing!
    
    Function that will apply all the steps for the automatization of the ICPMS data 
    processing:
        1) Read cps amd ppb data and %rsd 
        2) Compute sens = cps/ppb and plot it
        3) Apply the sensitivity correction
        4) Plot the new sensitivtiy plots (should be straight lines)
        5) Substract ICPMS blanks
        6) Save that data (to calibrate after, by hand, for the moment..)
    This function computes as well the std using quadratic uncertainty propagation. 
    Note the excel sheets containing the std also contain info on the ppb, but 
    that info is just nonsense, I didnt erase it not to have problem with
    dimensions!
    
    Important notes:
        . Ensure no blank is here!!
        . To do 2), its needed that the ppb data table begins with IS conc [ppb]! 
        Befpre there must be only the cps data, nothing else, no text nor anything!!
        . To do 3), you define the IS cases (which IS are to be used), so beware, 
                maybe your case its not (yet) defined!!  
        .For 6), you need to have installed xlsxwriter (conda install xslxwriter
                                                        for anaconda isntall)
        
    *Inputs:
        .df_cps: df containing the cps data and also the ppb data, in the classic 
        format. Must not contain the wash,  will give errors (divide by zero). 
        Take care of the names (like for the cps table), they are crutial 
        for the sens calc and  correction! Also about the format, not rows with 0s
        etc.  Take a lot of care!!!
        
        .IS_meas: array containing in a list the IS and its resolution, like how
            they appear in the isotopes column. 
            Default value:
                ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)',
                 'In115(MR)']
            That mean those isotopes were measured. If Ho165(MR) also measured, 
            just included it, and fine ;)
        .ICPblk_columns: np array containing the columns numbers where the ICPMS
        blanks are (blank and std 0ppt). Numbers from the excel (A = 0, B = 1, etc)
        .name_plot_L(M)R_bef(aft): name of the plots of the IS sensitivity before 
        and after doing the IS correction. 
        .excel_name: string that will be the name of the file containing the 
        output. Default: 'df_IS_Blks_corr.xlsx'
        
    *Output:
        .Dictionary containing the data applying the ICPMS blanks and sensitivity 
        corrections :)
    
    
    ###### To Do: #######
            .Automatize more stuff? such as the sens corrections, not by case 
                        something better, more general??
            .Output more things? Or less??
            :Delete the ppb data when not needed, so in excels it doesnt exist neither?
    '''
    
    ##### 0 ) Std calc ###########
    '''From the %rsd and the cps the std is trivial, %rsd = std/cps * 100 ==> 
    'std = cps*%rstd/100
    '''
    
    df_std = ICPMS_std_calculator(df_cps, df_rsd)           #Std of the cps
    
    ###### 1) IS sens calc ######
    df_IS_sens, df_IS_sens_std = IS_sens_calculator_plotter(df_cps, df_std, IS_meas, 
        name_IS_sens_LR_plot = name_plot_LR_bef, 
        name_IS_sens_MR_plot = name_plot_MR_bef)         
                #calculation the IS correction
                #I define names of the plot, so the other ones have the default name
    
    print('\n ########################')
    print('Step 1. done, computed IS sens :) ')
    print('#############################')
        #printing that tha step went good, would be good for debuggind!
        
        
    ###### 2) IS sens correction and new sens calc ##########
    df_IS_co, df_IS_co_std = IS_sens_correction(df_cps, df_std, df_IS_sens, 
                                        df_IS_sens_std, IS_meas)
                               #applying the IS correction to the cps data
    
    print('\n ########################')
    print('Step 2.1 done, applyed the IS correction :)) ')
    print('#############################')    
    
    df_IS_sens_co, df_IS_sens_co_std = IS_sens_calculator_plotter(df_IS_co, 
                            df_IS_co_std, IS_meas, 
                name_IS_sens_LR_plot= name_plot_LR_aft, 
                name_IS_sens_MR_plot= name_plot_MR_aft)    
                                #getting and plotting new IS sens
    
    print('\n ########################')
    print('Step 2.2. donce, copmuted the new IS sens :)) ')
    print('#############################')
    
    
    ##### 3)ICPMS Blk correction #########
    df_Blks_co, df_Blks_co_std = ICPMS_ICPMSBlanks_corrector(df_IS_co, 
                                    df_IS_co_std, ICPblk_columns)    
                                        #correcting for ICPMS blanks
    print('\n ########################')
    print('Step 3 (final). done, applyed the Blk correction :))) ')
    print('#############################')

    ##### 4) Saving and Output #########
    '''
    Here I want to save the df after IS correction, and after the Blk correction. 
    Both steps would be nice, for debugging!
    
    This requires installings xlsxwriter
    '''
    writer = pd.ExcelWriter(excel_name, engine = 'xlsxwriter')      #excel writer

    #########Blk correction data
    df_Blks_co.to_excel(writer, sheet_name = 'Blk_correction', startrow = 6, 
                        header = False, freeze_panes = (6, 1))            
    #saving to excel. I freeze row 6 and column 1, so I can see all the data 
    #in a good way :)                   
    #Chatgpt helped me to get the format I want, the one that Stefaan uses :)

    excel_sheet = writer.sheets['Blk_correction']
    bold_format = writer.book.add_format({'bold': True})      
    excel_sheet.write_row('B2', df_Blks_co.columns, bold_format)     
    #Write from cell B2 with the numer of columns. Note B2 is column B, row 2
    excel_sheet.write_row('B1', range(1, len(df_Blks_co.columns) + 1), bold_format)  
                                #2nd row with columns names
    excel_sheet.write_row('B3', [None] * len(df_Blks_co.columns))            
                                #row 3 empty
    excel_sheet.write_row('B4', ['Net <Int>'] * len(df_Blks_co.columns))         
                                #row 4 with a value repeated
    excel_sheet.write_row('B5', ['[cps]'] * len(df_Blks_co.columns))     
                                #row 5 with a value repeated

    #########std Blk correction data
    df_Blks_co_std.to_excel(writer, sheet_name = 'Blk_correction_std', startrow = 6, 
                            header = False, freeze_panes = (6, 1))            

    excel_sheet = writer.sheets['Blk_correction_std']
    bold_format = writer.book.add_format({'bold': True})      
    excel_sheet.write_row('B2', df_Blks_co_std.columns, bold_format)
    excel_sheet.write_row('B1', range(1, len(df_Blks_co_std.columns) + 1), bold_format)  
                        #2nd row with columns names
    excel_sheet.write_row('B3', [None] * len(df_Blks_co_std.columns))            
                        #row 3 empty
    excel_sheet.write_row('B4', ['Net Int std'] * len(df_Blks_co_std.columns))        
                            #row 4 with a value repeated
    excel_sheet.write_row('B5', ['[cps]'] * len(df_Blks_co_std.columns))     
                    #row 5 with a value repeated
    
    
    ##########Is correction data
    '''
    In this sheets I will also include the sensitivity calc, both, corrected and
    before the correction!
    '''
    df_IS_co.to_excel(writer, sheet_name = 'IS_correction', startrow = 6, 
                      freeze_panes = (6, 1), header = False)        
            #saving to excel in another sheet
            #THe start trow make that Co59 is on row 7, as it should be!
    df_IS_sens_co.to_excel(writer, sheet_name = 'IS_correction', 
                           startrow = 6 + df_IS_co.shape[0] + 2,header = False)
                        #putting the new IS sensitivity below the IS corrected data!!
    
    df_IS_sens.to_excel(writer, sheet_name = 'IS_correction', 
            startrow = 6 + df_IS_co.shape[0] + df_IS_sens_co.shape[0] + 4,header = False)
                        #putting the old sens below the new iS sensitivity!
    
    excel_sheet2 = writer.sheets['IS_correction']
    excel_sheet2.write_row('B2', df_IS_co.columns, bold_format)      
                #1st row with numbers of columns
    excel_sheet2.write_row('B1', range(1, len(df_IS_co.columns) + 1), bold_format) 
                    #2nd row with columns names
    excel_sheet2.write_row('B3', [None] * len(df_IS_co.columns))    #row 3 empty
    excel_sheet2.write_row('B4', ['<Int>'] * len(df_IS_co.columns))  
                        #row 4 with a value repeated
    excel_sheet2.write_row('B5', ['[cps]'] * len(df_IS_co.columns))     
                                    #row 5 with a value repeated
    ##Trial, to write stufff before the iS sens
    excel_sheet2.write_row('B' + str(6 + df_IS_co.shape[0] +2), 'IS sens, corrected',
                            bold_format)
    excel_sheet2.write_row('B' + str(6 + df_IS_co.shape[0] + df_IS_sens_co.shape[0] +4),
                           'IS sens', bold_format)
    #Does not work how I wanted, but enough, I do not need perfect, its for me xDDDD
    #Perfect is the enemy of good enough
    
    ############ IS correction std data
    df_IS_co_std.to_excel(writer, sheet_name = 'IS_correction_std', startrow = 6,
                          freeze_panes = (6, 1), header = False)        
                        #saving to excel in another sheet
    df_IS_sens_co_std.to_excel(writer, sheet_name = 'IS_correction_std', 
                               startrow = 6 + df_IS_co_std.shape[0] + 2, header = False)
                        #putting the new IS sensitivity below the IS corrected data!!
    
    excel_sheet2 = writer.sheets['IS_correction_std']
    excel_sheet2.write_row('B2', df_IS_co_std.columns, bold_format)      
                    #1st row with numbers of columns
    excel_sheet2.write_row('B1', range(1, len(df_IS_co_std.columns) + 1), bold_format)  
                    #2nd row with columns names
    excel_sheet2.write_row('B3', [None] * len(df_IS_co_std.columns))            
                    #row 3 empty
    excel_sheet2.write_row('B4', ['std'] * len(df_IS_co_std.columns))         
                    #row 4 with a value repeated
    excel_sheet2.write_row('B5', ['[cps]'] * len(df_IS_co_std.columns))     
                    #row 5 with a value repeated    
    
    
    #writer.save()                 #critical step, save the excel xD 
    writer.close()      #save was deprecated bro xD

    ################ Return ####################
    '''
    I return all the df:
        .ppb and ppb std IS corrected
        .ppb and std Blk corrected
    Althought I mainly need the excel to continue with the calib
    '''
    
    To_save = {'ppb_Blk_co' : df_Blks_co, 'ppb_std_Blk_co' : df_Blks_co_std,
               'ppb_Is_co' : df_IS_co, 'ppb_std_co': df_IS_co_std} #to return
    
    
    print('------------------------------------------------------')
    print('Finished the Data Process succesfully :) ')
    print('------------------------------------------------------')
    
    
    return To_save



#%%######## 1.11) ICPMS Isotope selector #############
#####################################
def ICPMS_Isotope_selector(df_cps, Isotopes):
    '''
    Function to choose from all the elements measured in the DF some selected 
    elements. This funciton will be used mainly to plot all those elements in 
    a single plot, since for that having one df witht he info is the best.
    
    
    *Inputs:
        .df_cps: df with the cps/ppb from the ICPMS measurement. Same format as
        "always" (jsut impleemnted xD), index are the isotopes then samples, 
        1streplicate, then 2nd
        .Isotopes: np.array containing the isotopes that I want, including the
        reoslution: np.array[('Sr84(LR)', 'Si28(MR)')]
        
        columns_blks: np.array([]) indicating the number of the columns contaning 
        blanks (std 0ppt and blank std). Ex: columns_blanks = np.array([1, 2])
        .df_IS_co: df containing the cps (also ppb, not needed but there it is),
        output from the IS sens correction funciton. The formatting is like the excel. 
    *Outputs:
        .df like the input one, but containing only the results for the desired isotopes
            

    '''
    #
    '''
    TO get the columns, given the relevant elements, I could do what I do
    in the cps/ppb function:
        I have a np aray with the elemnts to find (IS_meas), and I find them like that:
            
    value_to_find = IS_meas[i]
    matching_rows = df_cps_ppb_dat.loc[df_cps_ppb_dat.iloc[:,0] == value_to_find]  
    #give the full row!
    '''
    
    df_Nucl_rel = pd.DataFrame()         #Empty df to store the values

    for i in range(0, len(Isotopes)):     #Loop through all the relevant elements

       matching_rows = df_cps.loc[df_cps.index == Isotopes[i]]  #give the full row!
                   #rows whose elements are the relevant elements
           
           #TO store temporarily those values I create an auxiliary df
       df_aux = pd.DataFrame(data = matching_rows)

           #And I add that to the storing df
       df_Nucl_rel = df_Nucl_rel.append(df_aux, ignore_index= False)
           #ignore index false to store the index, key!

    ########### Return ###########
    return df_Nucl_rel        #return of the data





#%%############ 1.12) Kd calculaor #############
#####################################

def ICPMS_KdQe_calc (df_dat, df_dat_std, df_VoM_disol, df_m_be, 
                     df_VoM_disol_std = 0.0001, df_m_be_std = 0.0001, Nrepl = 2, 
                     ret_Co__Ceq = False):
    '''
    Function that will compute the distribution constant Kd and the adsorption 
    quantity Q_e from the ppb data obtained with ICPMS. 
    
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
    equilibirum since its kinetic, so Qe, Ke ==> Q(t), K(t). 
    
    ICPMS data sorted in replicates order: 1_1, 1_2,... 2_1, 2_2,...
    being sample 1 the blank (C_0)

    Necessary that the data contain no Div0, ensure in the excel by using the 
    iferror(operation, 0) function!
    
    Note this requires a df series with the volume, that you were not measuring
    in the first exp
    (up to 8/23). Note ten that units are involved!. If measuring mass ing and 
    volumes in L, Q


    *Inputs:
        .df_dat: dataframe containing the ICPMS data, the full data, with the 2 or 3 
            replicates. Should be Dfs corrected. Format: isotopes as index, 
            columns the samples, 1st 1st replicate, then 2nd replicate. Note 
            the 1st sample in each replicate must be the sample blank! This
            samples will be removed from the Qe and Kd df!
        .d_dat_std: similar df but with the std
        .df_VoM_disol: pd series containing the volume [mL] added to the bottle
            of the solution, BIC, or whatever. normally 50ml or the total 
            mass of the solution [g]. If df_dat in ppb, this must be the total
            mass so that Q_e is in g/g !
        .df_VoM_disol_std: same but with the std. Default value: 0.0001g
            (balance uncertainty)
        .df_m_bent: pd series contaning the mass of bentonite [g] in the bottle
        (normally 250mg)
        .df_m_bent_std: same but witht he std. Default: 0.0001g
        .Nrepl: number of replicates. Default value = 2. 3 also accepted
        ret_Co__Ceq: if True, returns a df with C_0 - C_eq = False
    
    *Outputs: dictionary with:
        .df with the Kd data
        .df with q_e data
        .Df qith rsq
        .Df with their std
        
        '''
    
    
    ########## 0) Precalcs ##########
    
    #1st lets convert the df to numeric, in case they are not:
    df_m_be = df_m_be.apply(pd.to_numeric)
    df_VoM_disol = df_VoM_disol.apply(pd.to_numeric)
    
    #
    '''
    To avoid the div0 error, which occurs when I have 0 as a values in the df, 
    which I have for all the elements that were not found by ICPMS, I can just 
    put NaN instead, since that will not give the Div0 error when computing Kd
    
    13/2/25, To do this, I was doing:
        df_dat.replace(0, np.nan, inplace=True)  #replace 0 values with NaN, 
                                                  #to avoid Div0 error!
        
    But that introduced a problem, NaN in the ppb values, since I do inplace =
    True. Best option is just to do that replcamente for the operation!
    '''
    
    
    #Splitting replicates
    N_sa = df_dat.shape[1]      #number of samples
    samples_per_repl = N_sa // Nrepl       #samples per repl = number of MS
    
    #lets get the replicates data in a list:
    repl_dat = []           #ppb/M data per repl
    repl_dat_std = []           #ppb/M data per repl
    repl_VoM_disol = []     #V disol (or mas)
    repl_m_be = []      #bentonite mass per replc
    
    for r in range(Nrepl):
        start = r * samples_per_repl        #1st sample of the replicate
        end = (r + 1) * samples_per_repl    #Last sample of the replicate
        repl_dat.append(df_dat.iloc[:, start:end])
        repl_dat_std.append(df_dat_std.iloc[:, start:end])
        repl_VoM_disol.append(df_VoM_disol.iloc[start+1:end]) #is a series, so less index
                #start+1 not to count the blank, which has 0 as value (same for mbe)
        repl_m_be.append(df_m_be.iloc[start+1:end]) #same 
    
    
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
    
    #storing variables
    repl_C0__Ceq = []    
    repl_C0__Ceq_std = []           #std
    repl_Qe = []
    repl_Qe_std = []                        #std
    repl_Kd = []
    repl_Kd_std = []
    repl_rsq = []
    repl_rsq_std = []
    
    for r in range(Nrepl):             #calcs per replicate
        '''
        Note that the 1st sample per replicate is the blank. There masses = 0,
        so the calc can not be made. Hence, I omit the blank in Qe, Kd calc
        '''
        C_0 = repl_dat[r].iloc[:,0]          #C0, blank concentration
        C_0_std = repl_dat_std[r].iloc[:,0]         #std of C_0
        Ceq__C0 = repl_dat[r].subtract(C_0, axis = 0 )
        Ceq__C0_std = np.sqrt( (repl_dat_std[r]**2).add(C_0_std**2, axis = 0) )
                    #std, using the pd functions, better and easier to understand!
        C0__Ceq = -Ceq__C0      #Obtaining C0-C_eq
        C0__Ceq_std = Ceq__C0_std 
        #
        #With that I can get Qe, Kd
        Qe = C0__Ceq.iloc[:,1:] * repl_VoM_disol[r] / repl_m_be[r]  #Qe
        Qe_std = np.abs(Qe)* np.sqrt( (df_VoM_disol_std / repl_VoM_disol[r])**2 + 
                              ( df_m_be_std / repl_m_be[r])**2 )
        Kd =  Qe.div(repl_dat[r].iloc[:,1:]).replace(0, np.nan)
                    #changing 0s in Conc for nan, not to have div0 errrors!
        Kd_std = np.abs(Kd) * np.sqrt( (Qe_std/Qe)**2 + 
        (repl_dat_std[r].iloc[:,1:] /(repl_dat[r].iloc[:,1:]).replace(0, np.nan) )**2 )
        
        rsq = C0__Ceq.div(C_0.values, axis = 0)*100           #rsq [%]
        rsq_std = np.abs(rsq) * np.sqrt( ( (C0__Ceq_std/ C0__Ceq)**2).add( 
                                        (C_0_std/C_0)**2, axis = 0) ) 
        
        #Before storing them, I will remova back the NaN, by doing Qe to 9999999999,
        #since it
        #could give errors (like that the number is easy noticeable!)
        Qe.fillna(99999999999, inplace = True)
        Qe_std.fillna(99999999999, inplace = True)
        Kd.fillna(99999999999, inplace = True)
        Kd_std.fillna(99999999999, inplace = True)
        
        #Finally, lets store it
        repl_C0__Ceq.append(C0__Ceq)   
        repl_C0__Ceq_std.append(C0__Ceq_std) 
        repl_Qe.append(Qe)
        repl_Qe_std.append(Qe_std)
        repl_Kd.append(Kd)
        repl_Kd_std.append(Kd_std)
        repl_rsq.append(rsq)
        repl_rsq_std.append(rsq_std)
        
    #Now we create a df out of them, concatenating them
        #we convert to numeric, in case it is needed
    df_Qe = pd.concat(repl_Qe, axis=1)
    df_Qe_std = pd.concat(repl_Qe_std, axis=1)
    df_Kd = pd.concat(repl_Kd, axis=1)  
    df_Kd_std = pd.concat(repl_Kd_std, axis=1)
    df_C0__Ceq = pd.concat(repl_C0__Ceq, axis = 1)
    df_C0__Ceq_std = pd.concat(repl_C0__Ceq_std, axis = 1)
    df_rsq = pd.concat(repl_rsq, axis = 1)
    df_rsq_std = pd.concat(repl_rsq_std, axis = 1)
    
    df_Qe = df_Qe.apply(pd.to_numeric)
    df_Qe_std = df_Qe_std.apply(pd.to_numeric)
    df_Kd = df_Kd.apply(pd.to_numeric)
    df_Kd_std = df_Kd_std.apply(pd.to_numeric)
    df_C0__Ceq =df_C0__Ceq.apply(pd.to_numeric) 
    df_C0__Ceq_std =df_C0__Ceq_std.apply(pd.to_numeric) 
    df_rsq =df_rsq.apply(pd.to_numeric) 
    df_rsq_std =df_rsq_std.apply(pd.to_numeric)    


    ########### 2) Return #############
    #Here the if for returning or not C_0 - C(t) applies
    
    results = {'Qe': df_Qe, 'Qe_std': df_Qe_std, 'Qe_%rsd' : df_Qe_std/df_Qe*100, 
               'Kd': df_Kd, 'Kd_std': df_Kd_std,'Kd_%rsd' : df_Kd_std/df_Kd*100, 
               'rsq' : df_rsq, 'rsq_std': df_rsq_std, 'rsq_%rsd' : df_rsq_std/df_rsq*100, }
    if ret_Co__Ceq:            #if you want to retreieve it
        results['C0-Ceq'] = df_C0__Ceq
        results['C0-Ceq_std'] = df_C0__Ceq_std
        results['C0-Ceq_%rsd'] = df_C0__Ceq_std/df_C0__Ceq*100
    return results
    
   


#%% ########## 1.13) Kd calculaor, Adsorption version #############
#####################################
def ICPMS_KdQe_calc_Ad (df_MS, df_MS_std, df_dat, df_dat_std, df_VoM_disol, 
        df_m_be, df_VoM_disol_std = 0.001, df_m_be_std = 0.0001, 
        ret_Co__Ceq = False, Nrepl = 3):
    '''
    Function that will compute the distribution constant Kd and the adsorption 
    quantity Q_e from the ppb data obtained with ICPMS. Their std deviation will
    also be computed, following the quadratic order propagation. Initial vs
    final concentration will also be returned, if desired. 
    
    The sorbed quantity Q_e is:
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
    Note that in the case of no equilibrium (eg Kinetics exp), Q_e, Kd ==> 
        Q_e(t), K(t). 
        
    Improvement made thanks to Espe (3/2024):
    Note that after the solid liquid sepparation,
    you retrieve a bit less of liquid, since the bentonite adsorbs some. This
    means that the concentration you measure with the ICP-MS is the concentration
    of this liquid, with mass m_liq_e <= m_liq_0. Then the Q_e should be 
    calculated like this, using that final mass!.
    
    
    Here we have different moher solutions, whihc is the main different
    from the other function, where all the samples had a common mother solution 
    (C_0). Here we have several C_0!

    Necessary that the data contain no Div0, ensure in the excel by using the 
    iferror(operation, 0) function!
    
    Note that Qe, Kd will be calculated for all the samples, even for 1 (blanks)!
    
    Note this requires a df series with the volume, that you were not measuring
    in the first exp (up to 8/23). Note ten that units are involved!. If 
    measuring mass ing and volumes in L, Q.


    *Inputs:
        .df_dat: dataframe containing the data, the full data, with the 3 
        replicates. Should be Dfs corrected
            Format: isotopes as index, columns the samples, 1st 1st replicate, 
            then 2nd replicate.
        .df_dat_std: similar df but with their std
        .df_MS: df containng the mother solution data, C_0
        .df_MS_std= similar df but with their std
        .df_VoM_disol: pd series containing the volume [mL] added to the bottle 
        of the solution, BIC, or whatever. normally 50ml OR the total mass of 
        the solution [g]. If df_dat in ppb, this must be the total mass
        so that Q_e is in g/g !
        .df_VOM_disol_std: value of the error of V or mL. Default: 0.001 [L] 
                Liquid assumed!!!
        .df_m_bent: pd series contaning the mass of bentonite [g] in the bottle 
        (normally 250mg)
        .df_m_bent_std: value of the error of m_be. Default: 0.0001g
        ret_Co__Ceq: if True, returns a df with C_0 - C_eq = False
        .Nrepl: number of replicates. Default: 3
    
    *Outputs:
        .Dictionary with Qe, Qe_std, Kd, Kd_std, rsq, rsq_std
        .If desired, dictionary including also Co-Ce and its std
        '''
    
    ########## 0) Precalcs ##########
    #1st, in case of the data is not numeric, lets convert it to numeric
    df_VoM_disol =  df_VoM_disol.apply(pd.to_numeric) 
    df_m_be = df_m_be.apply(pd.to_numeric) 

    '''
    To avoid the div0 error, which occurs when I have 0 as a values in the df, 
    which I have for all the elements that were not found by ICPMS, I can just
    put NaN instead, since that will not give the Div0 error when computing Kd.
    This is done but only in the Kd calc, to avoid any error (in C0-Ceq for ex)
    '''
    
    
    '''
    We will make it more efficient by doing the calcs on a replicate base, so
    in a for loop per replicate, the calcs will be done (gpt based)
    '''
    print('############# Qe calc comment: #########################')
    print('Are all the variable in order? It is mandatory since we will split them!')
    print('######################################\n')
    
    #Splitting replicates
    N_sa = df_dat.shape[1]      #number of samples
    samples_per_repl = N_sa // Nrepl       #samples per repl = number of MS
    
    #lets get the replicates data in a list:
    repl_dat = []           #ppb/M data per repl
    repl_dat_std = []           #[ppb/M] std of the data
    repl_VoM_disol = []     #V disol (or mas)
    repl_m_be = []      #bentonite mass per replc
    
    for r in range(Nrepl):
        start = r * samples_per_repl            #1st sample of the replicate
        end = (r + 1) * samples_per_repl        #Last sample of the replicate
        repl_dat.append(df_dat.iloc[:, start:end])          #dat
        repl_dat_std.append(df_dat_std.iloc[:, start:end])      #dat std
        repl_VoM_disol.append(df_VoM_disol.iloc[start:end]) 
            #it is a series, so less index
        repl_m_be.append(df_m_be.iloc[start:end])                   #same 
    
    #Now we can proceed with the calcs
    ########### 1) Calcs ###########
    '''
    The operations to perform are:
        1) C_0 - C_eq (>0)
        2) 1) * V/m = q_e
        3) 2)	1/C_eq = Kd
    
    I must treat the 3 experiments are different, I should substract the blank 1
    to the 1st emasurements and the 2 to the others. Since I ordered it in the 
    right way (1st replicacte 1, then replicate 2, I could) split it easily :D
            df.shape gives the shape of the df, n_rows, n_columns
    
    Note the df have number of samples * 3 replicates columns.
    
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
    
    #storing variables
    repl_C0__Ceq = []    
    repl_C0__Ceq_std = []           #std
    repl_Qe = []
    repl_Qe_std = []                        #std
    repl_Kd = []
    repl_Kd_std = []
    repl_rsq = []
    repl_rsq_std = []
    
    for r in range(Nrepl):             #calcs per replicate
        '''
        Note that I have N different mother solutions, which also means N different
        samples. I could do that with a for loop, but I found a better version. 
        I can ubstrcat df ignoring their indexes by doing df.values!
        '''
        Ceq__C0 = pd.DataFrame(repl_dat[r].values - 
            df_MS.values, index = df_MS.index, columns = repl_dat[r].columns)
        Ceq__C0_std = pd.DataFrame( np.sqrt(repl_dat_std[r].values**2 + 
            df_MS_std.values**2), index = df_MS.index, columns = repl_dat[r].columns)
        C0__Ceq = -Ceq__C0      #Obtaining C0-C_eq
        C0__Ceq_std = np.abs(Ceq__C0_std)      #Obtaining C0- std
        #With that I can get Qe, Kd
        Qe = C0__Ceq * repl_VoM_disol[r] / repl_m_be[r]
        Qe_std = np.abs(Qe) * np.sqrt( (df_VoM_disol_std / df_VoM_disol[r])**2 + 
                              ( df_m_be_std / df_m_be[r])**2 ) 
        Kd =  Qe.div(repl_dat[r] ).replace(0, np.nan)
                    #changing 0s in Conc for nan, not to have div0 errrors!
        Kd_std = np.abs(Kd) * np.sqrt( (Qe_std/Qe)**2 + 
        (repl_dat_std[r] /(repl_dat[r] ).replace(0, np.nan) )**2 )
        
        
        rsq = C0__Ceq / df_MS.values * 100
                #They have different row names, thats why the .values
        rsq_std = np.abs(rsq) * np.sqrt( (C0__Ceq_std/ C0__Ceq)**2 + 
                                        (df_MS_std/df_MS).values**2)
        
        #Before storing them, I will remova back the NaN, by doing Qe to 9999999999,
        #since it
        #could give errors (like that the number is easy noticeable!)
        Qe.fillna(99999999999, inplace = True)
        Qe_std.fillna(99999999999, inplace = True)
        Kd.fillna(99999999999, inplace = True)
        Kd_std.fillna(99999999999, inplace = True)
        
        #Finally, lets store it
        repl_C0__Ceq.append(C0__Ceq)
        repl_C0__Ceq_std.append(C0__Ceq_std)
        repl_Qe.append(Qe)
        repl_Qe_std.append(Qe_std)
        repl_Kd.append(Kd)
        repl_Kd_std.append(Kd_std)
        repl_rsq.append(rsq)
        repl_rsq_std.append(rsq_std)
        
    #Now we create a df out of them, concatenating them
        #we convert to numeric, in case it is needed
    df_Qe = pd.concat(repl_Qe, axis=1)
    df_Qe_std = pd.concat(repl_Qe_std, axis=1)
    df_Kd = pd.concat(repl_Kd, axis=1)  
    df_Kd_std = pd.concat(repl_Kd_std, axis=1)
    df_C0__Ceq = pd.concat(repl_C0__Ceq, axis = 1)
    df_C0__Ceq_std = pd.concat(repl_C0__Ceq_std, axis = 1)
    df_rsq = pd.concat(repl_rsq, axis = 1)
    df_rsq_std = pd.concat(repl_rsq_std, axis = 1)
    
    df_Qe = df_Qe.apply(pd.to_numeric)
    df_Qe_std = df_Qe_std.apply(pd.to_numeric)
    df_Kd = df_Kd.apply(pd.to_numeric)
    df_Kd_std = df_Kd_std.apply(pd.to_numeric)
    df_C0__Ceq =df_C0__Ceq.apply(pd.to_numeric) 
    df_C0__Ceq_std =df_C0__Ceq_std.apply(pd.to_numeric) 
    df_rsq =df_rsq.apply(pd.to_numeric) 
    df_rsq_std =df_rsq_std.apply(pd.to_numeric)     

    ########### 2) Return #############
    #Here the if for returning or not C_0 - C(t) applies
    
    results = {'Qe': df_Qe, 'Qe_std': df_Qe_std, 'Qe_%rsd' : df_Qe_std/df_Qe*100, 
               'Kd': df_Kd, 'Kd_std': df_Kd_std,'Kd_%rsd' : df_Kd_std/df_Kd*100, 
               'rsq' : df_rsq, 'rsq_std': df_rsq_std, 'rsq_%rsd' : df_rsq_std/df_rsq*100, }
    if ret_Co__Ceq:            #if you want to retreieve it
        results['C0-Ceq'] = df_C0__Ceq
        results['C0-Ceq_std'] = df_C0__Ceq_std
        results['C0-Ceq_%rsd'] = df_C0__Ceq_std/df_C0__Ceq*100
    return results
        
    
    return results


#%%####### 1.14) Mean/std of replicates calculator #############
#####################################


def ICPMS_MeanStd_calculator (df_dat, df_dat_std = None, Nrepl = 2, Debug = 0):
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
    
    A modification will be done, since, if replicates are consitent, the std
    of the average will be low. But maybe each replciate have a non negligible
    std, making the average std wrong. In this case, the std of the average will
    be the average of the std of the replicates. That will be done for each
    element.

    *Inputs:
        .df_dat: dataframe containing the data, the full data, with the 2 
        or 3 replicates. Should be Dfs corrected
            Format: isotopes as index, columns the samples, 1st 1st replicate, 
            then 2nd replicate.
        .df_dat_std: same but with the std. Default: None (so older code works)
        .Nrepl: number of replicates. Default value = 2. 3 also accepted
    
    *Output:
        .Dictionary containing df with
            .< >
            .std
            .%rsd
        The df have same column names: Data + number. Data and not sample, 
        not to mix it, since the Data 1 may be from replicates 2 for ex, in 
        case of Qe, kd. But could be from replicates 1 for variables like 
        Conc for ex
        
        '''
    
    #df_dat.replace(0, np.nan, inplace=True)     
                #replace 0 values with NaN, to avoid Div0 error!
    
    ######### 0) Data processing ##############
    #We are converting the data into a df, so if it is a series, will be converted
    #to a DF:
    is_series = isinstance(df_dat, pd.Series)
    if is_series:
        df_dat = df_dat.to_frame(name="value").T  # make it a 1-row DataFrame
        if df_dat_std is not None:  #do the same for df_std if it exist
            df_dat_std = df_dat_std.to_frame(name="value").T

    ncols = df_dat.shape[1]             #number of samples
    if ncols % Nrepl != 0:
        raise ValueError(f"Number of columns ({ncols}) not divisible by Nrepl={Nrepl}")

    # Reshape into (n_rows, N_sa, Nrepl)
    N_sa = ncols // Nrepl           #number of samples per repl
    colnames = [f"Data {i+1}" for i in range(N_sa)]     #names of the columns
        #for the new df
        
    ########### 1) compute and collect mean and std ###########
    #GPT based
    
    means, stds = [], []            #to store the values, mean and std of repl
    std_repl = []       #to accumulat the average of the std between replicates
    
    for i in range(N_sa):
        cols = [i + j*N_sa for j in range(Nrepl)] #column numbers for the 3 replicates
                    #of a given sample
        df_temp = df_dat.iloc[:, cols]  #df with the 3 replicates for each sample
        means.append(df_temp.mean(axis=1) )
        stds.append(df_temp.std(axis=1) ) #store the std of the mean
        
        if df_dat_std is not None:      #if dat_std is provided
            df_temp_std = df_dat_std.iloc[:, cols]      #same but for the std
            std_repl.append(df_temp_std.mean(axis = 1) )    #doing mean of stds
        else:      #no dat_std provided, so using the other value
             std_repl.append(df_temp.std(axis=1) )      #using the value of the < > calc
             
             
    #------ Final storing in df ------------
    df_mean = pd.concat(means, axis=1)
    df_mean.columns = colnames
    df_std = pd.concat(stds, axis=1)
    df_std.columns = colnames
    df_std_repl = pd.concat(std_repl, axis=1)
    df_std_repl.columns = colnames
    
    #Now I will chose the higher one, or std (from mean calc), or <stds (repl)>
    
    df_final_std = pd.DataFrame( np.maximum(df_std.values, df_std_repl.values),
        index=df_std.index, columns=colnames )      #the final std to be used
    
    # If original was Series, return as Series
    if is_series:
        df_mean = df_mean.iloc[0]
        df_final_std = df_final_std.iloc[0]


    ## Returning a dictionary
    # Compute %RSD safely (NaN where mean == 0)
    df_rsd = df_final_std / df_mean.replace(0, np.nan) * 100 
    
    if Debug:           #return both stds to compare them!
        Output = {'< >': df_mean, 'std': df_final_std, '%rsd' : df_rsd,
                  'std_repls' : df_std_repl, 'std_from_mean' : df_std}

    else:
        Output = {'< >': df_mean, 'std': df_final_std, '%rsd' : df_rsd}


    return Output
        
    
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
    

    
    
    
#     if type(df_dat) == pd.core.frame.DataFrame :         #If True data is a df

#         df_mean = pd.DataFrame()         #Empty df to store the values
#         df_std = pd.DataFrame()        #Empty df to store the std
    
#         if Nrepl == 2:          #Standard case, 2 replicates
#               df_1 = df_dat.iloc[ :, 0: round( ( df_dat.shape[1] ) / 2 ) ]      
#                           #1st replicate
#               df_2 = df_dat.iloc[ :, round( ( df_dat.shape[1] ) / 2 ) :  ]       
#                           #replicate 2
#         #
#               for i in range(df_1.shape[1]):   
#                   #loop thorugh all elements, but with index to work with 2 df
#                   df_temp = df_dat.iloc[:, [i, i+ df_1.shape[1] ] ]        
#                           #df containing the 2 replicates of the number i
            
#                 #TO store temporarily those values I create an auxiliary df
#                   df_aux1 = pd.DataFrame(data = df_temp.mean(1), columns = [i] )
#                   df_aux2 = pd.DataFrame(data = df_temp.std(1), columns = [i] )

#                 #And I add that to the storing df (add as columns)
#                   df_mean['Data ' + str(i +1) ] = df_aux1
#                   df_std['Data ' + str(i +1) ] = df_aux2
        

#         elif Nrepl == 3:            #3 replicates
#         #Gathering the replicates sepparately    
#             df_1 = df_dat.iloc[:, : round(df_dat.shape[1] / 3)]
#             df_2 = df_dat.iloc[:, round(df_dat.shape[1] / 3): 2*round(df_dat.shape[1] / 3)]
#             df_3 = df_dat.iloc[:, 2*round(df_dat.shape[1] / 3) :]
        
#             for i in range(df_1.shape[1]):         #loop thorugh all elements, but with index to work with 2 df
#                 df_temp = df_dat.iloc[:, [i, i+ df_1.shape[1], i+ 2 * df_1.shape[1] ] ]        
#                                                 #df containing the 3 replicates of the number i
            
#                 #TO store temporarily those values I create an auxiliary df
#                 df_aux1 = pd.DataFrame(data = df_temp.mean(1), columns = [i] )      #< >. 1 indicates compute by columns
#                 df_aux2 = pd.DataFrame(data = df_temp.std(1), columns = [i] )       # Std 

#                 #And I add that to the storing df (add as columns)
#                 df_mean['Data ' + str(i +1) ] = df_aux1
#                 df_std['Data ' + str(i +1) ] = df_aux2

#     elif type(df_dat) == pd.core.frame.Series :      #if data is a serie
#         df_mean = pd.Series()         #Empty df to store the values
#         df_std = pd.Series()        #Empty df to store the std
        
#         if Nrepl == 2:          #Standard case, 2 replicates
#               df_1 = df_dat.iloc[ 0: round( ( df_dat.shape[0] ) / 2 ) ]      #1st replicate
#               df_2 = df_dat.iloc[ round( ( df_dat.shape[0] ) / 2 ) :  ]       #replicate 2
#             #
#               for i in range(df_1.shape[0]):         #loop thorugh all elements, but with index to work with 2 df
#                 df_temp = df_dat.iloc[[i, i+ df_1.shape[0] ] ]        #df containing the 2 replicates of the number i
                
#                 #TO store temporarily those values I create an auxiliary df
#                 #df_temp.mean and .std gives mean and std
                
#                 df_mean['Data ' + str(i +1 ) ] = df_temp.mean()      #<>
#                 df_std['Data ' + str(i +1 ) ] = df_temp.std()         #std

        
#         elif Nrepl == 3:            #3 replicates
#         #Gathering the replicates sepparately    
#             df_1 = df_dat.iloc[: round(df_dat.shape[0] / 3)]
#             df_2 = df_dat.iloc[round(df_dat.shape[0] / 3): 2*round(df_dat.shape[0] / 3)]
#             df_3 = df_dat.iloc[ 2*round(df_dat.shape[0] / 3) :]
            
#             for i in range(df_1.shape[0]):         #loop thorugh all elements, but with index to work with 2 df
#                 df_temp = df_dat.iloc[ [i, i+ df_1.shape[0], i+ 2 * df_1.shape[0] ] ]        
#                                                     #df containing the 3 replicates of the number i

#                 #And I add that to the storing df (add as columns)
#                 df_mean['Data ' + str(i +1) ] = df_temp.mean()
#                 df_std['Data ' + str(i +1) ] = df_temp.std()    
    
    
#     else:       #weird data
#         print('############################')
#         print('\n Bro WTF? Wrong data type! pd.Series or pd.DataFrame only!')
#         print('############################')
        
#         df_mean = False
#         df_std = False
    
# ########### 2) Return #############
#     #Returning a dictionary
    
#     Output = {'< >': df_mean, 'std': df_std, '%rsd' : df_std/df_mean*100}

#     return Output    
           #return



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
                             return_leached = False, Nucl_rel = Isot_rel,
                             Nrepl = 3):
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
        .Nucl_rel: array containing the name of the relevant nuclides.
        Eg: np.array(['U238(LR)', 'Sr88(LR)'])
        .Nrepl: number of replicates. Default: 3
        
    *OUTPUTS
        .df_ppb_br: df with the ICPMS data with the bentonite contribution 
            removed (br). Note the procedural blank will still be there. Though
            they may not be useful, at least for the Qe function they will
        .df_ppb_br_std: df of the std of the ppb_br
    '''

    if Nrepl == 3:         #Case of 3 replicates 
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
    elif Nrepl == 2:           #2 replc
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
    print(df_leached.loc[Nucl_rel])
    print('######################################\n')

    print('########## Uncertainty of those: ########')
    print(df_leached_std.loc[Nucl_rel])
    print('######################################\n')
    
    ########### Output ###################
    
    if return_leached:      #True, so return it
        return df_ppb_br, df_ppb_std_br, df_leached, df_leached_std
    else:       #False, dont return it
        return df_ppb_br, df_ppb_std_br
    
    
#%% ######## 1.17 Ratio based correction Bentonite leached ##############
#################################################################    

def ICPMS_Removal_Bent_leach_ratio(data_dict, return_leached=False,
    Nucl_rel=None, Nrepl=3):
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
                    2_1-0_2
                    3_1-0_3
    And it will perform the following correction to the other samples:
        C_f corr = C_f * C_0/(C_0+C_leach)
    which would be like:
        1_1 corr = 1_1 * 0_1/(0_1+C_leach 1)
        1_2 corr = 1_2 * 0_2/(0_1+C_leach 1), 1_3 corr = 1_3 * 0_3/(0_3+C_leach 1),
        ...
        2_2 corr = 2_2 * 0_2/(0_2+C_leach 2), 2_3 cor = 2_3 * 0_3/(0_3+C_leac)
        ...
        
    *INPUTS
        .data_dict: dictionary containig the input df in the order:
            ppb/M data
            ppb/M std data
            MS ppb/M data
            MS ppb/M std data
            (keys not read)
        All DataFrames must have isotopes as index and samples as columns. They
        should be in order: 1,2,3,...
        .return_leached: boolean, to indicate wether df with the ppb (and 
            their std) leached by bentonite should be returned or not. Default:
            False
        .Nucl_rel: array containing the name of the relevant nuclides.
        Eg: np.array(['U238(LR)', 'Sr88(LR)']). For printing the df of Cleach
        .Nrepl: number of replicates. Default: 3
        
    *OUTPUTS
        .Dictionary containing the ppb, std, and %rsd values
    '''

# === 1. Extract input DataFrames ===
    df_dat, df_dat_std, df_MS, df_std_MS = list(data_dict.values())
    
    # === 2. Split replicates ===
    n_cols = df_dat.shape[1]
    samples_per_repl = n_cols // Nrepl
    repl_dat = []
    repl_std_dat = []

    for r in range(Nrepl):
        start = r * samples_per_repl
        end = (r + 1) * samples_per_repl
        repl_dat.append(df_dat.iloc[:, start:end])
        repl_std_dat.append(df_dat_std.iloc[:, start:end])

    # === 3. Compute leached concentration per replicate ===
    C_leach = []
    C_leach_std = []

    for r in range(Nrepl):
        leach = repl_dat[r].iloc[:, 0] - df_MS.iloc[:, 0]
        std_leach = np.sqrt(
            repl_std_dat[r].iloc[:, 0] ** 2 + df_std_MS.iloc[:, 0] ** 2
        )
        C_leach.append(leach)
        C_leach_std.append(std_leach)

    df_C_leach = pd.concat(C_leach, axis=1)
    df_C_leach.columns = [f'Repl {r+1}' for r in range(Nrepl)]

    df_C_leach_std = pd.concat(C_leach_std, axis=1)
    df_C_leach_std.columns = [f'Repl {r+1}' for r in range(Nrepl)]

    # === 4. Compute correction factors and apply them ===
    df_dat_br_list = []
    df_dat_std_br_list = []

    for r in range(Nrepl):
        leach = C_leach[r].values[:, None]
        leach_std = C_leach_std[r].values[:, None]

        corr = df_MS / (df_MS + leach)

        # Error propagation for correction factor
        delta_corr = corr * np.sqrt(
            (df_std_MS / df_MS) ** 2 +
            (df_std_MS**2 + leach_std**2) /
            (df_MS**2 + leach**2)
        )

        # Apply correction
        df_corr = repl_dat[r] * corr.values
        df_corr_std = df_corr * np.sqrt(
            (repl_std_dat[r] / repl_dat[r]) ** 2 +
            (delta_corr / corr).values ** 2
        )

        df_dat_br_list.append(df_corr)
        df_dat_std_br_list.append(df_corr_std)

    # === 5. Concatenate all corrected data ===
    df_dat_br = pd.concat(df_dat_br_list, axis=1)
    df_dat_std_br = pd.concat(df_dat_std_br_list, axis=1)

    # === 6. Reporting ===
    if Nucl_rel is not None:
        print('####### Leached concentrations (dat) of relevant nuclides:')
        print(df_C_leach.loc[Nucl_rel])
        print('########## Their uncertainties:')
        print(df_C_leach_std.loc[Nucl_rel])
        print('######################################\n')

    # === 7. Return ===
    result = {
        'dat_br': df_dat_br,
        'std_br': df_dat_std_br,
        '%rsd': df_dat_std_br/df_dat_br*100 #returning %rsd also!
    }

    if return_leached:      #add the df with the Cleached and its std
        result['C_leach'] = df_C_leach
        result['std_C_leach'] = df_C_leach_std
        result['%rsd_C_leach'] = df_C_leach_std/df_C_leach*100
    return result


    
#%% ######## 1.18 Cs correction ##############
#################################################################    

def ICPMS_Cs_correction(df_ppb, df_ppb_std, df_sens, 
                        columns_blks = np.array([1,2,3]),
                            Cs133_fis_ab = 41.34, Ba134_fis_ab = 19.07,
                            Ba136_fis_ab = 1.72, Ba137_fis_ab = 15.17,
                            Ba138_fis_ab = 64.04):
    '''
    This function will apply the corrections to compute, for a SNF leachate ICPMS:
        .Natural and fission Ba 136, 138 determination (Ba136,138 stable!)
        .Calc of Ba134,5,6,7 based on previous calcs
        .Fission-produced Cs133,4,5,7 (no natural Cs, it was removed with the
            ICPMS Blk corr)
    The data should be Xe corrected, since Xe interfere with Ba and Cs. ORIGEN
    calcs from the fuel are also need. Some rsd will be also printed throughout
    the function to check what makes the uncertainty rises a lot
    
    Note that for the std calcs, it was ASSUMED:
        i) NO std to natural abundances
        ii) NO std to ORIGEN calcs (fission abundances, etc)
        
    *INPUTS
        .df_ppb: df with the ppb data. 
        .df_ppb_std: df with the ppb std data
        .df_sens: df with the sensitivity (cps/ppb) of each nuclide(mass), as well
        as with its std
        .columns_blks: np.array([]) indicating the number of columns containing
            blks. From excel, so column A = 1
        . .._ab: fission abundances of those nuclei [%]
        
    *OUTPUTS
        .dicitonary containing 3 df:
            .df_ppb: including Ba, Cs 
            .df_ppb_std: including the rstd
            .df with the Cs abundances [%]
        
        
        
            #### To Do:
                .Link with abundance excel??
                .Problem with uncertainties, from fis ratio on, std> variable,
                since we do variable = 1- var 2, delta var = delta var 2, and 
                there the problem arises!
    '''

    ######## 0) 
    columns_blks = columns_blks - 1 
        #to adapt the counting system, python starts at0, not 1!
    print('\n------------------------------------------------')
    print('------------ Start fo the Cs calibration function ------------')
    print('Did you calibrate properly Ba masses? is mandatory!!!')
    print('-------------------------------------------------\n')


    ######### 1) Ratio Ba 138/Ba 136 calc
    rat_Ba138136_nat = 71.698/7.854        #nat abundance ratio of Ba138/ba 136
    
    rat_Ba138136 = df_ppb.loc['Ba138(LR)'] / df_ppb.loc['Ba136(LR)']
                                #measured ratio of Ba138/Ba136
    rat_Ba138136_std = rat_Ba138136 * np.sqrt(
        (df_ppb_std.loc['Ba138(LR)']/df_ppb.loc['Ba138(LR)'])**2 + 
        (df_ppb_std.loc['Ba136(LR)']/df_ppb.loc['Ba136(LR)'])**2 ) #std of the 
                #measured ratio
                
    print('------- %rsd of the Ba138/Ba136 ratio-----------')
    print((rat_Ba138136_std / rat_Ba138136 *100).round())
    print('-------------------------------------------')
    
    
    rat_Ba138136_fis = Ba138_fis_ab / Ba136_fis_ab  #fission ratio (ORIGEN)
        #no std, assumed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        
    #Excel #Checked!
    
    print('##############################################')
    print('Natural ratio Ba138/Ba136:' + f'{rat_Ba138136_nat: .3f}. Measured ratio:')
    #f for decimals, e for scientific notation
    print(rat_Ba138136.round(3))
    print('#################################')
    print('If they are close (9.something), indicates that little fission-produced Ba present')
    print('More change, such as 11.something indicates fission-produced Ba present')
    print('#####################################################\n')
    
    
    ########## 2) Calc of the fraction of natural/fission Ba:
    '''
    With the ratios I can compute the ratio of natural and fission Ba. Beware, this
    calc should only take place for radiaoct samples. For natural no, you set
    all to nat Ba!
    
    Note here, Stef set for natural samples, a ratio of 1, so all Ba is natural!
    '''
    #With the ratios I copmute the fraction of natural and fission. Beware, this
    
    frac_Ba136_nat = (rat_Ba138136 - rat_Ba138136_fis)/(
        rat_Ba138136_nat -rat_Ba138136_fis)     #fraction of nat Ba136
    frac_Ba136_nat_std = frac_Ba136_nat * rat_Ba138136_std/np.abs(rat_Ba138136 
                                                            - rat_Ba138136_fis)
                        #std
    '''
    Okay, that is not 100% true, since for the blanks, all the Ba136 is natural Ba
    (ASSUMPTION), there is no fission Ba. Hence, I will overwrite them to set
    1 as fraction, and 0 to std!
    '''
    frac_Ba136_nat.iloc[columns_blks] = 1        #no fission, is a blank!
    frac_Ba136_nat_std.iloc[columns_blks] = 0            #std
    
    
    frac_Ba136_fis = 1- frac_Ba136_nat      #Fission fraction = 1 - natural fraction
    frac_Ba136_fis_std = frac_Ba136_nat_std      #std
    
    print('--------- %rsd of fraction of Ba136 natural created:')
    print( (frac_Ba136_nat_std / frac_Ba136_nat*100).round() )
    print('%rsd of fraction of Ba136 fission created:')
    print( (frac_Ba136_fis_std / frac_Ba136_fis*100).round() )
    print('-----------------------------------------')
    '''
    Note rsd blow up here, in the fission rsd, since variable close to 0, being
    smaller than the uncertainty!
    if A = 1-B, delta A = delta B, but if A<<<, delta B = Delta A > A, making this
    results strange. At this point, I do not know what to do with this :/
    !!!!!!!!!!!!!
    !!!!!!!!!!!
    '''
    
    Ba136_nat = frac_Ba136_nat * df_ppb.loc['Ba136(LR)']    #nat Ba136!
    Ba136_nat_std = Ba136_nat * np.sqrt( (frac_Ba136_nat_std/frac_Ba136_nat)**2 +
            (df_ppb_std.loc['Ba136(LR)']/ df_ppb.loc['Ba136(LR)'])**2 ) #std
    
    Ba138_nat = Ba136_nat * rat_Ba138136_nat
    Ba138_nat_std = Ba138_nat * Ba136_nat_std/Ba136_nat  #std (rat nat no std given)
    
    Ba136_fis = frac_Ba136_fis * df_ppb.loc['Ba136(LR)']    #fiss Ba136
    Ba136_fis_std = Ba136_fis * np.sqrt((frac_Ba136_fis_std/frac_Ba136_fis)**2 + 
                     (df_ppb_std.loc['Ba136(LR)']/df_ppb.loc['Ba136(LR)'])**2 )
    
    Ba138_fis = Ba136_fis *rat_Ba138136_fis          #fiss Ba138
    Ba138_fis_std = Ba138_fis * Ba136_fis_std/Ba136_fis #std (fis ratio no std!)
    #Checked!
    
    #Printing rsd of Ba136, Ba138 natural and fission generated (both stable isotopes)
    
    print('----%rsd of natural-produced Ba136------------')
    print( (Ba136_nat_std/Ba136_nat * 100).round() )
    print('---------- %rsd of natural-produced Ba138-----------')
    print( (Ba138_nat_std/Ba138_nat * 100).round() )
    print('--------------------------------------------------')
    print('------------ %rsd of fission-produced Ba136')
    print( (Ba136_fis_std/Ba136_fis * 100).round() )
    print('------------ %rsd of fission-produced Ba138')
    print( (Ba138_fis_std/Ba138_fis * 100).round() )
            #Note they have same rsd, despite being different!
    
    ########### 3) Ba134,5,6,7, nat calc from Ba138
    #Using natural abundances data!
    
    Ba134_nat = 2.417 * Ba138_nat/71.698
    Ba134_nat_std = Ba134_nat * Ba138_nat_std/Ba138_nat      
            #std
    
    Ba135_nat = 6.592 * Ba138_nat/71.698
    Ba135_nat_std = Ba135_nat * Ba138_nat_std/Ba138_nat 
    
    # Ba136_nat = 7.854 * Ba138_nat/71.698
    # Ba136_nat_std = Ba136_nat * Ba138_nat_std/Ba138_nat 
    
    Ba137_nat = 11.232 * Ba138_nat/71.698
    Ba137_nat_std = Ba137_nat * Ba138_nat_std/Ba138_nat 
    #Checked! (std no checked, computing here for 1st time)


    #rsd printing
    print('--------- %rsd of nat Ba134 ---------')
    print( (Ba134_nat_std/Ba134_nat * 100).round() )
    print('--------- %rsd of nat Ba135 ---------')
    print( (Ba135_nat_std/Ba135_nat * 100).round() )
    # print('--------- %rsd of nat Ba136 ---------')
    # print( (Ba136_nat_std/Ba136_nat * 100).round() )
    print('--------- %rsd of nat Ba137 ---------')
    print( (Ba137_nat_std/Ba137_nat * 100).round() )
    print('----------------------------------------')
    
    
    ########## 4) Fission Ba134, 137 #####################
    #ORIGEN dependant!
    
    Ba134_fis = Ba134_fis_ab * ( Ba136_fis + Ba138_fis) / (Ba136_fis_ab + Ba138_fis_ab )
                #Ba134 fis
    Ba134_fis_std = Ba134_fis * np.sqrt(Ba136_fis_std**2 + Ba138_fis_std**2)/ (
        Ba136_fis + Ba138_fis)          #std
                
    Ba137_fis = Ba137_fis_ab * ( Ba136_fis + Ba138_fis) / (Ba136_fis_ab + Ba138_fis_ab )
            #Ba 137 fis
    Ba137_fis_std = Ba137_fis * np.sqrt(Ba136_fis_std**2 + Ba138_fis_std**2)/ (
        Ba136_fis + Ba138_fis)          #std        

    #rsd printing
    print('-------- rsd Ba134 fission generated ------------')
    print( (Ba134_fis_std/Ba134_fis *100).round() )
    print('-------- rsd Ba137 fission generated ------------')
    print( (Ba137_fis_std/Ba137_fis *100).round() )
    print('--------------------------------------------------')
    
    ############ 5) Cs (fission only, no natural) calc ########
    '''
    #Note no natural Cs since it was removed with the ICPMS Blk substraction,
    #so what remains is only radioactive!
    
    This is evident:
    Cs 134 fiss = total 134 (measured) - Ba134 nat - Ba134 fis
    and similar for Cs137, Cs135
    '''

    #
    Cs134_fis = df_ppb.loc['Ba134(LR)'] - Ba134_fis - Ba134_nat
    Cs134_fis_std = np.sqrt(df_ppb_std.loc['Ba134(LR)']**2 + Ba134_fis_std**2 
                            + Ba134_nat_std**2 )
    Cs137_fis = df_ppb.loc['Ba137(LR)'] - Ba137_fis - Ba137_nat
    Cs137_fis_std = np.sqrt(df_ppb_std.loc['Ba137(LR)']**2 + Ba137_fis_std**2 
                            + Ba137_nat_std**2 )
    #
    Cs135_fis = df_ppb.loc['Ba135(LR)'] - Ba135_nat
    Cs135_fis_std = np.sqrt(df_ppb_std.loc['Ba135(LR)']**2 + Ba135_nat_std**2)
                #no Ba 135 fis created!!
    #checked!
    
    #Cs rsd printing
    print('------------- Cs134 fission-generated rsd')
    print( (Cs134_fis_std/Cs134_fis *100).round() )
    print('------------- Cs135 fission-generated rsd')
    print( (Cs135_fis_std/Cs135_fis *100).round() )
    print('------------- Cs137 fission-generated rsd')
    print( (Cs137_fis_std/Cs137_fis *100).round() )
    
    
    '''
    Their abundances will also be computed, since they are really useful to 
    understand the data. We would expect similar Cs133 as Cs137 production, or
    slightly higher Cs137. Huge discrepancies would indicate that the data is not
    trustworthy
    '''

    Cs133_ab = df_ppb.loc['Cs133(LR)'] * 100 / (df_ppb.loc['Cs133(LR)'] + 
        Cs134_fis + Cs135_fis + Cs137_fis)
    Cs134_ab = Cs134_fis * 100 / (df_ppb.loc['Cs133(LR)'] + 
        Cs134_fis + Cs135_fis + Cs137_fis)
    Cs135_ab = Cs135_fis * 100 / (df_ppb.loc['Cs133(LR)'] + 
        Cs134_fis + Cs135_fis + Cs137_fis)
    Cs137_ab = Cs137_fis * 100 / (df_ppb.loc['Cs133(LR)'] + 
        Cs134_fis + Cs135_fis + Cs137_fis)
    
    df_Cs_ab = pd.DataFrame([Cs133_ab, Cs134_ab, Cs135_ab, Cs137_ab], 
                index = ['Cs133', 'Cs134', 'Cs135', 'Cs137']) #abundances [%]
    
    print('-------------------------------------------------------')
    print('Fission Cs abundances [%]:')
    print(df_Cs_ab)
    print('Ab of Cs133 should be <= Cs137. Huge discrepancies indicate untrustworthy data')
    print('-----------------------------------------------------------\n')
    
    
    
    ############# 6) Sens correction for Cs
    '''
    Note that the ppb of Cs were computed using the Ba sensitivity, so to correct
    for that, I can simply revert that, and multiply for the Cs133 sens:
        sens new = sens Ba 13X / sens Cs 133
        
        X = 4,5,7
        
    Fpor the std there I needed to use np.abs (what is right btw), sicne some values
    were negative!
        
    '''
    Cs134_fis_co = Cs134_fis * df_sens.loc[
        'Ba134(LR)']['Sens [cps/ppb]']/df_sens.loc['Cs133(LR)']['Sens [cps/ppb]']
    Cs134_fis_co_std = np.abs(Cs134_fis_co) * np.sqrt( (Cs134_fis_std/Cs134_fis)**2 + 
     (df_sens.loc['Ba134(LR)']['std [cps/ppb]']/
      df_sens.loc['Ba134(LR)']['Sens [cps/ppb]'])**2 + (
          df_sens.loc['Cs133(LR)']['std [cps/ppb]']
          /df_sens.loc['Cs133(LR)']['Sens [cps/ppb]'])**2 )
    
    Cs135_fis_co = Cs135_fis * df_sens.loc[
        'Ba135(LR)']['Sens [cps/ppb]']/df_sens.loc['Cs133(LR)']['Sens [cps/ppb]']
    Cs135_fis_co_std = np.abs(Cs135_fis_co) * np.sqrt( (Cs135_fis_std/Cs135_fis)**2 + 
     (df_sens.loc['Ba135(LR)']['std [cps/ppb]']/
      df_sens.loc['Ba135(LR)']['Sens [cps/ppb]'])**2 + (
          df_sens.loc['Cs133(LR)']['std [cps/ppb]']
          /df_sens.loc['Cs133(LR)']['Sens [cps/ppb]'])**2 )
    
    Cs137_fis_co = Cs137_fis * df_sens.loc[
        'Ba137(LR)']['Sens [cps/ppb]']/df_sens.loc['Cs133(LR)']['Sens [cps/ppb]']
    Cs137_fis_co_std = np.abs(Cs137_fis_co) * np.sqrt( (Cs137_fis_std/Cs137_fis)**2 + 
     (df_sens.loc['Ba137(LR)']['std [cps/ppb]']/
      df_sens.loc['Ba137(LR)']['Sens [cps/ppb]'])**2 + (
          df_sens.loc['Cs133(LR)']['std [cps/ppb]']
          /df_sens.loc['Cs133(LR)']['Sens [cps/ppb]'])**2 )
          #excel checked!

    '''
    Another useful check is to compare the total Cs obtained with the one predicted
    from the abundancies obtained from ORIGEN. They should not disagree too much
    .Huge discrepancies would indicate untrustowrthy data (measurement data)
    '''
    Cs_tot = df_ppb.loc['Cs133(LR)'] + Cs134_fis_co + Cs135_fis_co + Cs137_fis_co
    Cs_tot_std = np.sqrt(df_ppb_std.loc['Cs133(LR)']**2 + Cs134_fis_co_std**2 +
                         Cs135_fis_co_std**2 + Cs137_fis_co_std**2) #std
    
    #Origen using Cs133
    Cs_tot_OR133 = df_ppb.loc['Cs133(LR)'] / Cs133_fis_ab * 100 #Cs expected from ORIGEN
    Cs_tot_OR133_std = Cs_tot_OR133 * df_ppb_std.loc['Cs133(LR)'] / df_ppb.loc['Cs133(LR)']
    rat_Cstot_OR133 = Cs_tot/Cs_tot_OR133
    
    #Origen using Cs137
    Cs_tot_OR137 = Cs137_fis_co / Cs137_ab * 100 #Cs expected from ORIGEN
    Cs_tot_OR137_std = Cs_tot_OR137 * Cs137_fis_co_std / Cs137_fis_co
                    #Assumes no error on Cs137abundance!!!
    rat_Cstot_OR137 = Cs_tot/Cs_tot_OR137
    
    print('--------------------------------------------------------------#')
    print('Total Cs measured [ppb]: ')
    print(Cs_tot.round(3) )
    print('--------------------------------------------------\n')
    print('Total Cs from ORIGEN using Cs133 [ppb]: ')
    print(Cs_tot_OR133.round(3))
    print('--------------------------------------------------\n')
    print('Ratio total Cs measured/ORIGEN using Cs133: ')
    print(rat_Cstot_OR133.round(3))
    print('--------------------------------------------------\n')
    print('Total Cs from ORIGEN using Cs137 [ppb]: ')
    print(Cs_tot_OR137.round(3))
    print('--------------------------------------------------\n')
    print('Ratio total Cs measured/ORIGEN using Cs137: ')
    print(rat_Cstot_OR137.round(3))
    print('Huge discrepancies (>= 2/3) indicate untrustworthy measurements, use ORIGEN data!')
    print('--------------------------------------------------------------#')
    print('------------------------ End of the function ----------------')
    print('--------------------------------------------------\n')

    
    
    ############## 7) Output ##################
    '''
    Okay, I will return the df_ppb, but I will add it the info. I will add all, 
    and with time I will know if I need more or less info xD
    
    Cs tot will not be saved!!!
    
    Note I all (LR) to all the names, which will be useful for plotting, since
    I plot and save name removing 4 last digits [:-4], which removed (LR)
    
    The %rsd will also be computed and stored!
    '''
    
    ##### ppb
    df_ppb.loc['Ba134(LR)_nat'] = Ba134_nat
    df_ppb.loc['Ba134(LR)_fis'] = Ba134_fis
    df_ppb.loc['Cs134(LR)_fis'] = Cs134_fis_co
    df_ppb.loc['Ba135(LR)_nat'] = Ba135_nat
    df_ppb.loc['Cs135(LR)_fis'] = Cs135_fis_co
    df_ppb.loc['Ba136(LR)_nat'] = Ba136_nat
    df_ppb.loc['Ba136(LR)_fis'] = Ba136_fis
    df_ppb.loc['Ba137(LR)_nat'] = Ba137_nat
    df_ppb.loc['Ba137(LR)_fis'] = Ba137_fis
    df_ppb.loc['Cs137(LR)_fis'] = Cs137_fis_co
    df_ppb.loc['Ba138(LR)_nat'] = Ba138_nat
    df_ppb.loc['Ba138(LR)_fis'] = Ba138_fis
    # The total Cs will also be given as output!
    #df_ppb.loc['Cs_tot(LR)'] = Cs_tot
    #df_ppb.loc['Cs_tot_ORIGEN133(LR)'] = Cs_tot_OR133
    #df_ppb.loc['Cs_tot_ORIGEN137(LR)'] = Cs_tot_OR137
    
    #### ppb_std
    df_ppb_std.loc['Ba134(LR)_nat'] = Ba134_nat_std
    df_ppb_std.loc['Ba134(LR)_fis'] = Ba134_fis_std
    df_ppb_std.loc['Cs134(LR)_fis'] = Cs134_fis_co_std
    df_ppb_std.loc['Ba135(LR)_nat'] = Ba135_nat_std
    df_ppb_std.loc['Cs135(LR)_fis'] = Cs135_fis_co_std
    df_ppb_std.loc['Ba136(LR)_nat'] = Ba136_nat_std
    df_ppb_std.loc['Ba136(LR)_fis'] = Ba136_fis_std
    df_ppb_std.loc['Ba137(LR)_nat'] = Ba137_nat_std
    df_ppb_std.loc['Ba137(LR)_fis'] = Ba137_fis_std
    df_ppb_std.loc['Cs137(LR)_fis'] = Cs137_fis_co_std
    df_ppb_std.loc['Ba138(LR)_nat'] = Ba138_nat_std
    df_ppb_std.loc['Ba138(LR)_fis'] = Ba138_fis_std    
    #
    #df_ppb_std.loc['Cs_tot(LR)'] = Cs_tot_std
    #df_ppb_std.loc['Cs_tot_ORIGEN133(LR)'] = Cs_tot_OR133_std
    #df_ppb_std.loc['Cs_tot_ORIGEN137(LR)'] = Cs_tot_OR137_std
    
    #rsd
    df_rsd = df_ppb_std/df_ppb * 100
    
    ###### Returning #########
    #A dictionary with the 3 df will be returned, to keep it more gathered
    output = {'dat': df_ppb, 'std': df_ppb_std, '%rsd': df_rsd,
              'Cs_Ab': df_Cs_ab}
    return output



#######################################################################
# %% ######### 1.19) Isotopes to elements ###########################
######################################################################
def ICPMS_Isot_to_Elem(df, df_std, Debug = False, Radiaoct_here = False):
    '''
    Function that will take an ICP-MS datasheet (in df format) and will merge 
    all the isotopes to have elemental data. That is:
        .Merge all (LR): U(LR) = sum (U233,U234,...), Eu(LR), etc
        .Merge all (MR): Si(MR) = sum (Si28, Si29..), U(MR), etc
    
    We can simply sum; [Si] = [Si28] + [Si29] + [Si30]? Yes:
    g Si/g tot = (g Si28+gSi29+g Si30)/gtot = [Si28] + [Si29] + [Si30]
    
    ACtually not for all the elements, since some sammes have interferences. For 
    those, you should not take those masses into acccount, rather use the other
    isotopes to stimate elemental conc: for Sr, Sr87 has interefernces (Rb87), so
    if we do [Sr] = sum (Sr), the value will be high. YOu should do (Stef)
        [Sr] = sum (isot without inter) / sum (wt% of those isotopes) * 100
        
    That I will do. Hence I need to read wt%, from the excel. This funciton also
    contains the list with all the interferences, from ICPMS excels (Calib sheets)
    
    Note that calculations would also needed for the other interefering isotopes.
    Eg, Mo92. Int he excel it is writen Zr92, so for Zr is computed with the wt%. 
    But for Mo would also need to be done, since not all the isotopes are written,
    Mo92 is missing (interferni).
    
    Uncertainty computed using quadratic propagation.
        
    Created by scholargpt, adapted by me. AI works for me bro!
    Function checked with excel (F vs NF)
    
    *Inputs
        .df: df containing the icpms data. Typical format;
            -indexes are isotopes
            -Each column is a sample
        .df_std: same but with their stds
        .Radiaoct_here: boolean to indicate wheather radioactive isotopes could be 
            found or not. TO use wt% data with those or not
        .Debug: boolean to indicate if you want the output not considering the 
            interferences
    *Outputs
        .Dictionary with:
            .df
            .df_std
            .df_%rsd
        
    '''
    
    #---- 0) Pre-requisites ------------
    '''
    The most complete list of intereference masses are listed below. From that,
    we might not measure the same. WE need to select the ones we measured
    '''
    Isob_interf = ['Sr87(LR)','Zr92(LR)', 'Zr94(LR)', 'Zr96(LR)', 'Ru102(LR)',
                   'Ru104(LR)', 'Sn120(LR)', 'Sn122(LR)', 'Sn122(LR)', 
                   'Sb123(LR)','Te124(LR)', 'Te124(LR)','Ce142(LR)',
                   'Nd144(LR)','Nd148(LR)','Sm150(LR)','Sm152(LR)','Sm154(LR)',
                   'Gd156(LR)','Gd158(LR)','Gd160(LR)','Dy162(LR)','Dy164(LR)',
                   'Er168(LR)','Yb170(LR)','Yb174(LR)','Yb176(LR)',
                   'Ti48(MR)','Ti50(MR)','Cr54(MR)','Ni58(MR)','Zn64(MR)',
                   'Ge70(MR)','Ge74(MR)','Se76(MR)',
                   #Repeatitions beggin for MR
                   'Sr87(MR)','Zr92(MR)', 'Zr94(MR)', 'Zr96(MR)']
            #Element having isotopes with interferences!! In115 excluded (IS)
            #This are isotope list in icpms excel
            #From excel spot, comparing Isotpes used for calib
            
    Elem_Interf = ['Rb(LR)', 'Mo(LR)', 'Pd(LR)', 'Te(LR)', 'Rb(LR)', 'Nd(LR)',
                   'Sm(LR)', 'Gd(LR)', 'Dy(LR)','Er(LR)','Yb(LR)','Hf(LR)',
                   'Lu(LR)','Ca(MR)','V(MR)','Cr(MR)','Fe(MR)','Ni(MR)',
                   'Zn(MR)','Ge(MR)']   #Interfering elements. These are the 
        #elements with conflicintg isotopes, but that not appear in the icpms index list
        #eg: Sr87 appear and is Sr87 Rb87, hence Rb is written here.
        
    
    #Are we dfealing with radioactive material? If so use at we radioact
    
    if Radiaoct_here == True:
            excel_dat = At_we_rad       #use data with radioactive material
    else:
            excel_dat = At_we_nat   #use data with natural material
     
    #The interefernces present in the give df we can get easily:
    Interf_here = []     #to sotre the intereferences present in the df given
    for index in df.index:          #for each isotoe
        if index in Isob_interf:        #if present in the list, store it
            Interf_here.append(index)
    
    
    # Extract element and resolution using regex
    def parse_label(label):
        '''
        Intermediate function that, given ICP-MS nuclide, like U238(LR),
        will return "U(LR)", getting element and resolution (string)
        '''
        
        match = re.match(r"([A-Za-z]+)[0-9]*(?:\((LR|MR)\))(?:_.*)?", label)
            #Match pattern
            #r"[A-Za-z] to find any letter
            #r"\d for a digit
            #(?:_.*)? for _something at the end
        if match:       #If a match is found:
            element = match.group(1)
            resolution = match.group(2)
            return f"{element}({resolution})"  #This is U(LR) for ex
        else:
            return None  # Handle malformed data gracefully

    #I need to put the index as a column to do this, a new df is create not
    #to alter the original one:
    df2 = df.copy()
    df2['Element(Res)'] = df2.index.to_series().apply(parse_label)
            #create the parsed label: U(LR), Si(MR), etc as index
    df2_std = df_std.copy()    
    df2_std['Element(Res)'] = df2_std.index.to_series().apply(parse_label)    
    # Drop rows with unparsed labels
    df_clean = df2.dropna(subset=['Element(Res)'])
        #This essentially removes Ar40Ar40(MR) and the U238O16(MR) and (LR)
    df_clean_std = df2_std.dropna(subset=['Element(Res)'])

    #Grouping by it, to get elemental indexes
    df_grouped = df_clean.groupby('Element(Res)')
    df_grouped = df_clean.groupby('Element(Res)', sort = False).sum()
    #   numeric_only=True)
   #group by the index. Adding the values. With num_only you avoid errors
   #of addding string values and so!
   #sort = False not so sort them alphabetically, keeping previous order
   #thats the sum of all the isotopes. Forgetting about interf
    
    '''
    Now the magic should start. For non interferences elemtns, sum them, what I
    was doing before. For interferences, different, using the %at data.
    
    XYou need to identify interference-free isotopes for the interference cases
    
    '''
    df_Elem = pd.DataFrame(index = df_grouped.index,
                    columns = df.columns)    #Df to fill with elemental data
    df_Elem_std = pd.DataFrame(index = df_grouped.index,
                    columns = df.columns)    #Df to fill with elemental data std
    #
    for index in df_Elem.index:        #loop throught all elements in the df
        Isotopes = df_clean.loc[df_clean['Element(Res)'] == index] 
                        #df with the isotopes for the given element
        Isotopes_std = df_clean_std.loc[df_clean_std['Element(Res)'] == index] 
        #
        if any(Isotope in Interf_here for Isotope in Isotopes.index
               ) or index in Elem_Interf:   #check if any 
                #of the isotopes of that element is a interfering isotope
            #(index) in the df, OR if the element is in the interfering eleemnt list
            #for the interfering isotopes not in the icpms index.
            #the 
            #The interferences are:
            Iso_inter = [iso for iso in Isotopes.index if iso in Interf_here]
                    #that have the isotopes with interferences. We remove them:
            Isotopes_inter_free = Isotopes.drop(Iso_inter, axis = 0)
                            #isotopes removing the interferences
            Isotopes_inter_free_std = Isotopes_std.drop(Iso_inter, axis = 0)
            #
            #Then I need to do sum isotopes * 100/sum abundances
            At_mass = excel_dat.loc[[x[:-4] for x in Isotopes_inter_free.index],
                        "wt%"]          #getting wt% of those isotopes (inter free)
                #Now we can compute the elemental conc:
            df_Elem.loc[index] = Isotopes_inter_free.drop('Element(Res)', 
                            axis = 1).sum(axis = 0)/ At_mass.sum(axis = 0) * 100
                    #elemental conc = sum isotopes/wt% * 100. I remove the elem(REs) column
            df_Elem_std.loc[index] = np.sqrt(Isotopes_inter_free_std.drop('Element(Res)', 
                            axis = 1)**2).sum(axis = 0) / At_mass.sum(axis = 0) * 100
                            #quadratic std prop
        
        else:       # no intereferences, just sum the isotopes
            df_Elem.loc[index] = Isotopes.drop('Element(Res)', 
                                               axis = 1).sum(axis = 0).values
                    #summing the isotopes. I need to remove the eleemnt column!
            df_Elem_std.loc[index] = np.sqrt(Isotopes_std.drop('Element(Res)', 
                                               axis = 1)**2).sum(axis = 0).values
    
    
    ############## Return
    #df_Elm is not numeric, so lets convert to it
    df_Elem = df_Elem.apply(pd.to_numeric)  
    df_Elem_std = df_Elem_std.apply(pd.to_numeric)  
    
    if Debug == True:
        return {'dat' : df_Elem, 'dat_interferences_not_removed': df_grouped} 
        #Note df_grouped > df_Elem, which amkes makes considering the overstimation
        #for the case with intereferences. Function must be checked, but looking good!
    else:
        return {'dat': df_Elem, 'std': df_Elem_std, 
                '%rsd': df_Elem_std/df_Elem * 100}
    
    
#%% ------------- 1.19.1) Isotopes to elements, dictionary version ##############
#The same but for a dictionary as input:
#Removed, as no longer aplciable
    
#################################################################
#%% ########## 1.20 )ICPMS Homogeneizer ###########################
###############################################

    
def ICPMS_Homogenize(df_ref, df, Return_extra_mass = 0):
    '''
    Function that will homogenize ICPMS df. Imagine that for 2 experiments you
    measured different masses. The masses of 1 df will be homogenized to match
    the masses of other df, called ref.
    
    In case that in the df_ref there are masses that are not present in the df,
    a warning will be returned. NaN will be placed then.
    
    *Inputs
        -df_ref: reference df wioth ICMPMS measurements (row/mass, column per sample)
        -df: df to homogenize
        -Return_extra_mass = 0. Boolean to indicate whether to return a df with the
         extra masses (1) or not
        
    *Outputs
        -df_hom: df homogenized
        -df_extra: the extra elements, which do not appear in df_ref (if any)
    '''
    
    #Target masses to retain
    Masses_ref = df_ref.index.tolist()
    Masses_df = df.index.tolist()
    
    #Warning if missing any target nuclides
    missing = [mass for mass in Masses_ref if mass not in Masses_df]
    if missing:
        print('################# Warning ######################################')
        print('Some masses of the reference df not measured in the given df:')
        print(missing)
        print('\n Giving NaN in the new df for those masses')
        print('################################################################\n')
        #warnings.warn(f'The following masses from df_ref are missing in df:{missing}',
                     #stacklevel = 2)
        
    #Get the homogeniezed df
    df_homo = df.loc[df.index.intersection(Masses_ref)].reindex(Masses_ref)
    
    #Retrieve also the extra data (if any)
    extra_masses = [mass for mass in Masses_df if mass not in Masses_ref]
    df_extra = df.loc[extra_masses]
            #this could be empty. 
    
    ############# Returning #######
    #If df_extra is empty it will not be returned
    if Return_extra_mass:       #if True, so we want to return it
        if df_extra.empty:      #to return it, first they should exist xD  
            return df_homo, df_extra
        
    else:                   #not desired to retrieve df_extra
        return df_homo
    
    

#################################################################
#%% ########## 1.21 )ICPMS Get activity ###########################
###############################################
    
def ICPMS_Get_Activity (df_ppb, df_ppb_std ):
    '''
    Function that, starting from the df with the ICPMS concentrations, will
    compute the totala activity of the desired radionuclides. Assuming concentration
    in ppb. Note that the Activity (A) is related to the specific activity (a) as:
        A = m * a
    
    Since we do not have m, rather we have concentrations, g /gtot, we would get
    the activity in Bq/g tot! This funciton will perform the calculation for all
    the isotopes defined in the excel "Rad_dat_DLA". For the others, nothing will 
    be done, removing them from the output df
    
    The uncertainties will also be given, ASSUMING no uncertainty in specific 
    activities (reasonable since ppb std will be way bigger)
    
    This function will also print the total activity, and its std
        
    *Inputs:
        .df_ppb: df with the ICPMS data in the typical format:
            .index = masses (U238(LR), etc)
            .each column is a sample
        .df_ppb_std: df with the ppb std, in the same format
        
    *Output
        .dictionary with 2 df containing the activity and its ucnertainty
    '''
    
    
    ################ 0) Pre calcs
    N_A = 6.022e23                              #[part/mol] Avogadro Number
    
    #Gettind theisotopes of the df, removving (LR)/(MR)
    
    base_isotopes = df_ppb.index.str.extract(r"^([A-Z][a-z]?\d+)")[0]
    #base_isotopes = [a[:-4] for a in df_ppb.index]
# Map to specific activity (Rad_dat must have isotopes as index)


    a_values = Rad_dat['a (Bq/g)'].reindex(base_isotopes).to_numpy() 
                #specific activity

    # Compute activity: conc(ppb) × SA(Bq/g) × 1e-9 (ppb → g/g)
    A = df_ppb.multiply(a_values * 1e-9, axis=0)  
    
    A_std = df_ppb_std.multiply(a_values * 1e-9, axis=0)  #std of A
    
    #The columns with no radioactive data will be removed (have NaN as values)
    A.dropna(axis = 0, inplace = True, how = 'all')
            #how = all indicates only removing when all values are NaN ! 
            #(error in Astd otherwise, since for BIC NaN in ppb
    A_std.dropna(axis = 0, inplace = True, how = 'all')
    
    #Finally we will also get the total activity:
    A.loc['A_tot[Bq/gtot]'] = A.sum(axis = 0, skipna = True)
    A_std.loc['A_tot[Bq/gtot]'] = A_std.sum(axis = 0, skipna = True)
    
    #Lets print the total activity, with no decimals
    print('#-------------------------------------------------------#')
    print('Total activity/total mass [Bq/g tot]:')
    print(A.loc["A_tot[Bq/gtot]"].round() )
    print('----- And its uncertainty: -------------\n')
    print(A_std.loc["A_tot[Bq/gtot]"].round() )
    print('#------------------End of the function --------------------------#')
    ############## Returning ######################
    #THe rsd will also be computed an returned
    
    A_rsd = A_std/A * 100           #rsd
    
    #Returning a df
    result = {'Act [Bq/gtot]': A, 'Delta[Act[Bq/gtot]]': A_std, '%rsd': A_rsd}
    
    return result
    
    
    
    
#---------------------------------------------------------------------------
#%% ------------ 1.22) Statisticals correlation tests ----------------------------------
#--------------------------------------------------------------------------

def ICPMS_Correlation_test(data, element_list, type = 'Spe'):
    '''
    Funciton that will check if there is correlation betweeen ICPMS datasets.
    I will give elements to compare, and the funciton will create a matrix
    comparing all of them. Eg:
        
                U(LR) Sr(LR) Si(MR)
        U(LR)   ----   a      b
        Sr(LR)  a     ---    c
        Si(MR)  b      c    ---
        
    Note the diagonal elements are not valid, same comparison, reporting max 
    correlation. Also that matrix is symmetric, a_ij = a_ji for each i,j.
    
    The test to eprform could be either (or both)
        :Pearson. Linear test
        .Spearman: nonlinear
        
    Both are based on examininng covariances:
        cov (x,y)/ (std(x)*std(y) )
    
    Pearson states if linear relationships obtained. Spearman if they are 
    monotonic (follow similar pattern). hence, the models should kinda agree,
    at least do not have different correlation sign. IF so, you could state
    that no clear correlation obtained.
    
    THe p-value will also be returned, which indicate how likely the value. 
    The p-value is the probability of obtaining a correlation coefficient 
    (or other test statistic) at least as extreme as the one you observed, 
    if the null hypothesis were true.

    Null hypothesis (H₀): “There is no correlation between the two variables” 
            (true r = 0).
    Alternative hypothesis (H₁): “There is a correlation” (true r ≠ 0).

    So, the p-value answers:
        “If U and Si were truly unrelated, how likely is it that I would see
        an r this strong just by chance?”
    Typically p < 0.05 → statistically significant.
    
    In geochemistry / environmental sciences, people often use the following 
    guidelines (for Pearson or Spearman):
        ∣r∣<0.3 → weak or negligible
        0.3 ≤ ∣r∣<0.5 → moderate
        0.5 ≤ ∣r∣<0. strong
        ∣r∣≥0.7 → very strong
    
    
    *Inputs
        .data: df with the ICPMS data to compare. INdex are elements, 
                        columns samples
        
        .element list: list containing the elements to compare. Eg: 
            ['U(LR)', 'Sr(LR)', 'Si(MR)']
        .type: string to indicate which test:
                'Spe' ---> Spearman (non linear)
                'Pea' ---- > pearson (linear)
                'Both' ---> do both!
            Default: 'Spe'
        from the ICPMS results.
    
    *Output
        .Dictionary with 2 df
            .df with the statistical correlation r
            .df with the p-value
        If both is stated, apaprt from the 4 df, 2 per method, a df with r*r
        is also returned
    
    '''

    #--------------- Empty df creation to store the results
    df_r = pd.DataFrame([], index = element_list, columns = element_list)
                        #empty df creation, to store the results, r
    df_p = pd.DataFrame([], index = element_list, columns = element_list)
                         #empty df creation, to store the results, p
                                            
    
    if type == 'Spe':                                           #Spearson test
        for i in range(df_r.shape[0]):               #loop thoruhg rows
            for j in range(df_r.shape[0]):          #loop thourhg columns
                df_r.iloc[i,j] = spearmanr(data.loc[element_list[i] ].values, 
                                          data.loc[element_list[j] ].values )[0]
                df_p.iloc[i,j] = spearmanr(data.loc[element_list[i] ].values, 
                                          data.loc[element_list[j] ].values )[1]
    elif type == 'Pea':                                     #Pearson test
        for i in range(df_r.shape[0]):       #loop thoruhg rows
            for j in range(df_r.shape[0]):          #loop thourhg columns
                df_r.iloc[i,j] = pearsonr(data.loc[element_list[i] ].values, 
                                          data.loc[element_list[j] ].values )[0]
                df_p.iloc[i,j] = pearsonr(data.loc[element_list[i] ].values, 
                                          data.loc[element_list[j] ].values )[1]
                #
    elif type == 'Both':             #2 both analysis
        #Need to create the storing df for the 2nd test
        df_r2 = pd.DataFrame([], index = element_list, columns = element_list)
                            #empty df creation, to store the results, r
        df_p2 = pd.DataFrame([], index = element_list, columns = element_list)
        #
        for i in range(df_r.shape[0]):                  #loop thoruhg rows
            for j in range(df_r.shape[0]):          #loop thourhg columns
                df_r.iloc[i,j] = pearsonr(data.loc[element_list[i] ].values, 
                                          data.loc[element_list[j] ].values )[0]
                df_p.iloc[i,j] = pearsonr(data.loc[element_list[i] ].values, 
                                          data.loc[element_list[j] ].values )[1]
                df_r2.iloc[i,j] = spearmanr(data.loc[element_list[i] ].values, 
                                          data.loc[element_list[j] ].values )[0]
                df_p2.iloc[i,j] = spearmanr(data.loc[element_list[i] ].values, 
                                          data.loc[element_list[j] ].values )[1]                
                
    else:
        print('-------------------------------------------')
        print('Wrong type writen! Only Spe, Pea of Both available!')
        print('-------------------------------------------\n')


    #-------- Return -------------------------
    
    if type == 'Both':                                  #return both tests!
        return {'Pea, r' : df_r.apply(pd.to_numeric) , 
        'Pea, p-value': df_p.apply(pd.to_numeric),
        'Spe, r' : df_r2.apply(pd.to_numeric) , 
        'Spe, p-value': df_p2.apply(pd.to_numeric),
        'r_Spe * r_Pea': df_r.apply(pd.to_numeric) * df_r2.apply(pd.to_numeric)}  
                #r*r showed to see if they predict different relations!
        print('------------------------------------------------------')
        print('The product of r will be returned. When r*r <0, the correlations' +
              ' methods are disagreeing !')
        print('-------------------------------------------\n')
    #
    else:                                             #return only one
        return {'r' : df_r.apply(pd.to_numeric) , 
            'p-value': df_p.apply(pd.to_numeric) }        #returning a dictionary
    
    
    
    
    
    
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
                      pre_save_name = 'Conc_rsd', Nucl_rel = Isot_rel, 
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
        .Nucl_rel: array containing the name of the relevant elemtns, which 
        are the elements that will be saved in a specific folder. D
        efault value: (see above in the script)   
        .Plot_everything: boolean stating if we plot everything or only the 
        relevant elements (in Nucl_rel). 
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
        if df_1.index[i] in Nucl_rel or plot_everything == True: 
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
            if df_1.index[i] in Nucl_rel:  #if the element is relevant
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
                      pre_save_name = 'Conc_rsd', Nucl_rel = Isot_rel, 
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
        .Nucl_rel: array containing the name of the relevant elemtns, which 
        are the elements that will be saved in a specific folder. D
        efault value: (see above in the script)   
        .Plot_everything: boolean stating if we plot everything or only the 
        relevant elements (in Nucl_rel). 
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
        if df_1.index[i] in Nucl_rel or plot_everything == True: 
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
            if df_1.index[i] in Nucl_rel:  #if the element is relevant
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
                      pre_save_name = 'Conc_rsd', Nucl_rel = Isot_rel, 
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
        .Nucl_rel: array containing the name of the relevant elemtns, which 
        are the elements that will be saved in a specific folder. D
        efault value: (see above in the script)   
        .Plot_everything: boolean stating if we plot everything or only the 
        relevant elements (in Nucl_rel). 
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
        if df_1.index[i] in Nucl_rel or plot_everything == True: 
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
            if df_1.index[i] in Nucl_rel:  #if the element is relevant
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
                   Nucl_rel = Isot_rel, 
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
        .Nucl_rel: array containing the name of the relevant elemtns, which 
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
            if index in Nucl_rel:  #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
                plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                plt.title(pre_title_plt + index[:-4], fontsize=22, wrap=True)           #title
                plt.plot(x[:int(len(x)/2)], row[:int(len(x)/2)], 'bo--', 
                     markersize = Markersize, label = 'Repl_1') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x[int(len(x)/2):], row[int(len(x)/2) :], 'ro--', 
                     markersize = Markersize, label = 'Repl_2') 
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
                         markersize = Markersize, label = 'Repl_1') 
                    plt.plot(x[int(len(x)/2):], row[int(len(x)/2):], 'ro--', 
                         markersize = Markersize, label = 'Repl_2') 
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
            if index in Nucl_rel:  #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
                plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                plt.title(pre_title_plt + index[:-4], fontsize=22, wrap=True)           #title
                plt.plot(x[:int(len(x)/3)], row[:int(len(x)/3)], 'bo--', 
                     markersize = Markersize, label = 'Repl_1') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x[int(len(x)/3): 2* int(len(x)/3)], row[int(len(x)/3) :2* int(len(x)/3)], 'ro--', 
                     markersize = Markersize, label = 'Repl_2') 
                plt.plot(x[2* int(len(x)/3):], row[2* int(len(x)/3) :], 'go--', 
                     markersize = Markersize, label = 'Repl_3') 
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
                     markersize = Markersize, label = 'Repl_1') 
                    #+1 needed since the df contain a row with the column names!
                    plt.plot(x[int(len(x)/3): 2* int(len(x)/3)], row[int(len(x)/3) :2* int(len(x)/3)], 'ro--', 
                     markersize = Markersize, label = 'Repl_2') 
                    plt.plot(x[2* int(len(x)/3):], row[2* int(len(x)/3) :], 'go--', 
                     markersize = Markersize, label = 'Repl_3') 
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
                    pre_save_name = 'Conc', Nucl_rel = Isot_rel ):
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
        .Nucl_rel: array containing the name of the relevant elemtns, which
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
        if df_cps['Sard'].index[i] in Nucl_rel:        #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
            plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
            plt.title(pre_title_plt + df_cps['Sard'].index[i][:-4], fontsize=22, wrap=True)           #title
            #PLot bentonite 1, Sard
            plt.plot(x['Sard'][:int(len(x['Sard'])/2)], df_cps['Sard'].loc[df_cps['Sard'].index[i] ][:int(len(x['Sard'])/2)], 'o--', color = Bent_color['Sard'],
                     markersize = Markersize, label = 'Repl_1 S') 
                    #+1 needed since the df contain a row with the column names!
            plt.plot(x['Sard'][int(len(x['Sard'])/2):], df_cps['Sard'].loc[df_cps['Sard'].index[i] ][int(len(x['Sard'])/2):], 'o--', color = Bent_color['Sard'],
                     markersize = Markersize, label = 'Repl_2 S') 
            #PLot bentonite 2, T
            plt.plot(x['Tur'][:int(len(x['Tur'])/2)], df_cps['Tur'].loc[df_cps['Sard'].index[i] ][:int(len(x['Tur'])/2 )], 'o--', color = Bent_color['Tur'],
                     markersize = Markersize, label = 'Repl_1 T') 
                    #+1 needed since the df contain a row with the column names!
            plt.plot(x['Tur'][int(len(x['Tur'])/2):], df_cps['Tur'].loc[df_cps['Sard'].index[i] ][int(len(x['Tur'])/2 ):], 'o--', color = Bent_color['Tur'],
                     markersize = Markersize, label = 'Repl_2 T') 
            #PLot bentonite 3, BK
            plt.plot(x['BK'][:int(len(x['BK'])/2)], df_cps['BK'].loc[df_cps['Sard'].index[i] ][:int(len(x['BK'])/2 )], 'o--', color = Bent_color['BK'],
                     markersize = Markersize, label = 'Repl_1 BK') 
                    #+1 needed since the df contain a row with the column names!
            plt.plot(x['BK'][int(len(x['BK'])/2):], df_cps['BK'].loc[df_cps['Sard'].index[i] ][int(len(x['BK'])/2 ):], 'o--', color = Bent_color['BK'],
                     markersize = Markersize, label = 'Repl_2 BK') 
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
                     markersize = Markersize, label = 'Repl_1 S') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x['Sard'][int(len(x['Sard'])/2):], df_cps['Sard'].loc[df_cps['Sard'].index[i] ][int(len(x['Sard'])/2):], 'o--', color = Bent_color['Sard'],
                     markersize = Markersize, label = 'Repl_2 S') 
                    #PLot bentonite 2, T
                plt.plot(x['Tur'][:int(len(x['Tur'])/2)], df_cps['Tur'].loc[df_cps['Sard'].index[i] ][:int(len(x['Tur'])/2 )], 'o--', color = Bent_color['Tur'],
                     markersize = Markersize, label = 'Repl_1 T') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x['Tur'][int(len(x['Tur'])/2):], df_cps['Tur'].loc[df_cps['Sard'].index[i] ][int(len(x['Tur'])/2 ):], 'o--', color = Bent_color['Tur'],
                     markersize = Markersize, label = 'Repl_2 T') 
                    #PLot bentonite 3, BK
                plt.plot(x['BK'][:int(len(x['BK'])/2)], df_cps['BK'].loc[df_cps['Sard'].index[i] ][:int(len(x['BK'])/2 )], 'o--', color = Bent_color['BK'],
                     markersize = Markersize, label = 'Repl_1 BK') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x['BK'][int(len(x['BK'])/2):], df_cps['BK'].loc[df_cps['Sard'].index[i] ][int(len(x['BK'])/2 ):], 'o--', color = Bent_color['BK'],
                     markersize = Markersize, label = 'Repl_2 BK') 
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
            Nucl_rel = Isot_rel, Logs = 0 ):
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
        .Nucl_rel: array containing the name of the relevant elemtns, which 
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
                if y_1.index[i] in Nucl_rel or plot_everything == True:      #if the element is relevant
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + y_1.index[i][:-4], fontsize=22, wrap=True)     #title
                    if Blank_here:      #if Blank here ==> 1st colum is blk
                        plt.plot(x_1[1:], y_1.loc[y_1.index[i] ][1:], 'bo--', 
                                 markersize = Markersize, label = 'Repl_1')          #repl 1
                        plt.plot(x_2[1:], y_2.loc[y_2.index[i] ][1:], 'ro--', 
                                 markersize = Markersize, label = 'Repl_2')          #repl 2                                             
                        plt.hlines(y_1.loc[y_1.index[i] ][0], min(x_1), max(x_1) , label = 'Blk_1', color = 'b' )
                        plt.hlines(y_2.loc[y_2.index[i] ][0], min(x_2), max(x_2) , label = 'Blk_2' , color = 'r' )                                                                           
                    else:                   #No blank!
                        plt.plot(x_1, y_1.loc[y_1.index[i] ], 'bo--', 
                                 markersize = Markersize, label = 'Repl_1')          #repl 1
                        plt.plot(x_2, y_2.loc[y_2.index[i] ], 'ro--', 
                                 markersize = Markersize, label = 'Repl_2')          #repl 2
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
                if y_1.index[i] in Nucl_rel or plot_everything == True:      #if the element is relevant
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + y_1.index[i][:-4], fontsize=22, wrap=True)     #title
                    if Blank_here:      #if Blank here ==> 1st colum is blk
                        plt.plot(x_1[1:], y_1.loc[y_1.index[i] ][1:], 'bo--', 
                                 markersize = Markersize, label = 'Repl_1')          #repl 1
                        plt.plot(x_2[1:], y_2.loc[y_2.index[i] ][1:], 'ro--', 
                                 markersize = Markersize, label = 'Repl_2')          #repl 2
                        plt.plot(x_3[1:], y_3.loc[y_2.index[i] ][1:], 'go--', 
                                 markersize = Markersize, label = 'Repl_3')          #repl 3                                               
                        plt.hlines(y_1.loc[y_1.index[i] ][0], min(x_1), max(x_1) , label = 'Blk_1', color = 'b' )
                        plt.hlines(y_2.loc[y_2.index[i] ][0], min(x_2), max(x_2) , label = 'Blk_2', color = 'r' )
                        plt.hlines(y_3.loc[y_3.index[i] ][0], min(x_3), max(x_3) , label = 'Blk_3', color = 'g' )                                                                             
                    else:                   #No blank!
                        plt.plot(x_1, y_1.loc[y_1.index[i] ], 'bo--', 
                                 markersize = Markersize, label = 'Repl_1')          #repl 1
                        plt.plot(x_2, y_2.loc[y_2.index[i] ], 'ro--', 
                                 markersize = Markersize, label = 'Repl_2')          #repl 2
                        plt.plot(x_3, y_3.loc[y_2.index[i] ], 'go--', 
                                 markersize = Markersize, label = 'Repl_3')          #repl 3 
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
                if y_1.index[i] in Nucl_rel or plot_everything == True:      #if the element is relevant
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + y_1.index[i][:-4], fontsize=22, wrap=True)     #title
                    if Blank_here:      #if Blank here ==> 1st colum is blk
                        plt.plot(x_1.loc[x_1.index[i]][1:], y_1.loc[y_1.index[i] ][1:], 'bo--', 
                                 markersize = Markersize, label = 'Repl_1')          #repl 1
                        plt.plot(x_2.loc[x_2.index[i]][1:], y_2.loc[y_2.index[i] ][1:], 'ro--', 
                                 markersize = Markersize, label = 'Repl_2')          #repl 2                                           
                        plt.hlines(y_1.loc[y_1.index[i] ][0], min(x_1.loc[x_1.index[i]]), 
                                   max(x_1.loc[x_1.index[i]]) , label = 'Blk_1', color = 'b' ) 
                        plt.hlines(y_2.loc[y_2.index[i] ][0], min(x_2.loc[x_2.index[i]]), 
                                   max(x_2.loc[x_2.index[i]]) , label = 'Blk_2', color = 'r' )                                                                            
                    else:                   #No blank!
                        plt.plot(x_1.loc[x_1.index[i]], y_1.loc[y_1.index[i] ], 'bo--', 
                                 markersize = Markersize, label = 'Repl_1')          #repl 1
                        plt.plot(x_2.loc[x_2.index[i]], y_2.loc[y_2.index[i] ], 'ro--', 
                                 markersize = Markersize, label = 'Repl_2')          #repl 2
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
                if y_1.index[i] in Nucl_rel or plot_everything == True:      #if the element is relevant
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + y_1.index[i][:-4], fontsize=22, wrap=True)     #title
                    if Blank_here:      #if Blank here ==> 1st colum is blk
                        plt.plot(x_1.loc[x_1.index[i]][1:], y_1.loc[y_1.index[i] ][1:], 'bo--', 
                                 markersize = Markersize, label = 'Repl_1')          #repl 1
                        plt.plot(x_2.loc[x_2.index[i]][1:], y_2.loc[y_2.index[i] ][1:], 'ro--', 
                                 markersize = Markersize, label = 'Repl_2')          #repl 2                                           
                        plt.plot(x_3.loc[x_3.index[i]][1:], y_3.loc[y_3.index[i] ][1:], 'go--', 
                                 markersize = Markersize, label = 'Repl_3')          #repl 3       
                        plt.hlines(y_1.loc[y_1.index[i] ][0], min(x_1.loc[x_1.index[i]]), 
                                   max(x_1.loc[x_1.index[i]]) , label = 'Blk_1', color = 'b' )
                        plt.hlines(y_2.loc[y_2.index[i] ][0], min(x_2.loc[x_2.index[i]]), 
                                   max(x_2.loc[x_2.index[i]]) , label = 'Blk_2', color = 'r' )                                                                           
                        plt.hlines(y_3.loc[y_3.index[i] ][0], min(x_3.loc[x_3.index[i]]), 
                                   max(x_3.loc[x_3.index[i]]) , label = 'Blk_3', color = 'g' )                                                                          

                    else:                   #No blank!
                        plt.plot(x_1.loc[x_1.index[i]], y_1.loc[y_1.index[i] ], 'bo--', 
                                 markersize = Markersize, label = 'Repl_1')          #repl 1
                        plt.plot(x_2.loc[x_2.index[i]], y_2.loc[y_2.index[i] ], 'ro--', 
                                 markersize = Markersize, label = 'Repl_2')          #repl 2                                           
                        plt.plot(x_3.loc[x_3.index[i]], y_3.loc[y_3.index[i] ], 'go--', 
                                 markersize = Markersize, label = 'Repl_3')          #repl 3       
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
    
    

##############################################################################
#%%### 1.19) ICPMS plotter blank Average of replicates, N datasets ############
##############################################################################

def ICPMS_Plotter_mean_blk_N (
    x_list, std_x_list,
    y_list, std_y_list,
    element_index,
    x_label, y_label,
    labels=None, colors=None, Fmts = None,
    folder_name='Plots', pre_title_plt='Concentration of ',
    pre_save_name='Conc', Nucl_rel=Elem_rel,
    Logscale=False, Blank_here=False, plot_everything=False, font_size= Font ):
    '''
    Function that will plots of the data from the ICPMS (cps) vs another variable, 
    initially time, the cps and the rstd, for the 2 bentonites, plotting the average 
    values ideally (output of average computer). It will also allow to plot the 
    blank as a horizontal line, if desired.
    
    *Inputs:
        .x_list: list with the x variable (df or series)
        .std_x_list: list with the x std variables (df or series). If None is provided,
        no error will be assigned
        .y_list: list with the y variable
        .std_y_list: list with the y std variables
        .x_label: string that will be the x label for the plot (for math stuff, 
                                    use $$. eg: '$\Delta t[h]$')
        .y_label: string that will be the y label for the plot
        .label: list with the labels. The hline will be: Label name + MS (if 
            there is blank)
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
        .Nucl_rel: array containing the name of the relevant elemtns, which are
            the elements that will be saved in a specific folder. Default 
            value: Elem_rel
        .Logscale: string to say if you want the x and y axis in logscale or not.
            Default: False
        .font_size: font_size for the plot. Default: 18
        .colors: array with the colors to use. Eg: ['red', 'green', 'blue',
                'orange', 'olive', 'pink', 'purple', 'yellow']. Defalut: None (random)
        .Fmt: array defining markers and connecting line. Default: None (random
                defined in the function). Eg of linesytles:
            ['-', '--', '-.', ':', '-.'].Eg of markers:
                ['o','s', '+', 'v', '^', '*', 'D', 'p']
        
                                    
    *Outputs:
        .Plots (saving them) of the x and df_mean_cps data, cps vs x!
    
    
    ### TO DO: ####
	.Plot 2 blk lines, <>+- std? OPtional, since for Qe I do not have blank, 
    but for conce I do!
    '''
    
    
   # === 1. Folder Setup ===
    path_main = os.path.join(os.getcwd(), folder_name)
    path_relevant = os.path.join(path_main, 'Relevants')
    os.makedirs(path_relevant, exist_ok=True)

    # === 2. Dataset Setup ===
    N = len(x_list)
    labels = labels or [f'Data {i+1}' for i in range(N)]
    colors = colors or ['blue', 'red', 'green', 'orange'][:N]
    Fmts = Fmts or ['o', 's', 'v', '^', 'p', '*', 'D', 'x'][:N] 
        #linestyles;: ['--', '-.', ':','-']
    
    # === 3. Resolve Accessors ===
    '''
    in case std x has no error (None), I will asign a low error, but not 0, so
    the loops can be executed. If I state 0, they would not work
    '''
    
    if std_x_list == None:  #setting a error close to 0, but not being 0
        std_x_list = [x/9999999 for x in x_list]
    
    is_series = isinstance(x_list[0], pd.Series)
    def get_x(x, i): return x if is_series else x.loc[i]
    def get_std(std, i): return std if is_series else std.loc[i]

    # === 4. Plotting Loop ===
    t_start = tr.time()

    for i in element_index:
        element_name = i#[:-4]
        relevant = i in Nucl_rel if Nucl_rel is not None else False

        if relevant or plot_everything:
            plt.figure(figsize=(11, 8))
            plt.title(pre_title_plt + element_name, fontsize=22, wrap=True)

            for k in range(N):      #loop plotting, for all the datasets
                x = get_x(x_list[k], i)
                y = y_list[k].loc[i]
                sx = get_std(std_x_list[k], i)
                sy = std_y_list[k].loc[i]

                Color = colors[k % len(colors)]
                Fmt = Fmts[k % len(Fmts)]
                Label = labels[k]

                if Blank_here:
                    plt.hlines(y[0], min(x), max(x), color=Color, linestyle='-', 
                               label=Label + ' MS')
                    plt.errorbar(x[1:], y[1:], yerr=sy[1:], xerr=sx[1:], 
                        color=Color, label=Label, fmt = Fmt, markersize=7)
                else:
                    plt.errorbar(x, y, yerr=sy, xerr=sx, 
                                 color=Color, label=Label, fmt = Fmt, markersize=7)

            plt.xlabel(x_label, fontsize=font_size)
            plt.ylabel(y_label, fontsize=font_size)
            plt.tick_params(axis='both', labelsize=font_size)
            if Logscale:
                plt.xscale('log')
                plt.yscale('log')
            plt.minorticks_on()
            plt.grid(which='minor', linestyle=':', linewidth=0.5)
            plt.grid(which='major')
            plt.legend(fontsize=font_size)

            # Save
            safe_name = element_name.replace("(LR)", "_LR").replace("(MR)", "_MR")
                    #replacing (MR) with _MR, not to save ( in the name
            save_path = os.path.join(
                path_relevant if relevant else path_main,
                f"{pre_save_name}_{safe_name}.png"
            )
            plt.savefig(save_path, format='png', bbox_inches='tight')
            plt.close()

    # === 5. Report Timing ===
    t_run = tr.time() - t_start
    print('###############################################')
    print(f'Plotting completed for {N} datasets in {t_run:.2f} seconds.')
    print('###############################################')
    
    

##############################################################################
#%%### 1.20) ICPMS multibar plotter ############
##############################################################################

def ICPMS_MultiBar_plotter(df, df_std, Elements, b = 0.2,
    Xlabel = 'X axis', Ylabel = 'Y axis', Title= 'PLot', Savename = 'Plot1'):
    '''
    Function that will do a multibar plot, with the number of abrs you give (Nbars),
    of a given df, chosing the elements you desire. The df_std shoudl also be included!
    
    
    *Inputs:
        .df/df_std: df with the icpms data
        .Elements: Array containing the elemtns you want to include. Eg:
                ['Sr(LR)', 'Si(MR)']
        .b: blank space between the bar groups. b < 1
        .Xlabel/Ylabel: string with the label for the x/y axis. Default: 'X/Y axis'
        .Title: string with the title of the plot. Defalt: 'Plot'
        .Savename: string with the name of the file to save (png). Eg: 'Plot1'
    
    *Outputs:
        No outputs, jsut generate a plot!
    
    '''

    #Parameters for the multibar plot
    Nbars = len(Elements)       # number of bars to plot
    w = (1 - b) / Nbars             #Width of each bar
    offsets = (np.arange(Nbars) - (Nbars - 1) / 2) * w  #displacement between the bar
            # Si N es par: [-3.5w, -2.5w, ..., +3.5w]
            # Si N es impar: [-3w, -2w, ..., +3w]
    X_axis =  np.linspace(1, df.shape[1], num =  df.shape[1]) #plotting all columns!
    
    
    # === Plot ===
    plt.figure(figsize=(12, 8))
    plt.title(Title, fontsize=22, wrap=True)

    for i, elem in enumerate(Elements):
        plt.bar(X_axis + offsets[i], 
            df.loc[elem],
            yerr=df_std.loc[elem],
            width=w, edgecolor="black",
            label=elem,
            align='center')
    
    plt.ylabel(Ylabel, fontsize= Font)              #ylabel
    plt.xlabel(Xlabel, fontsize= Font)   
    #plt.xticks(X_axis, Dict_el_MS['dat'].columns, rotation=90)
    #plt.yscale('log') 
    plt.legend(fontsize = Font)
    plt.tick_params(axis='both', labelsize=Font)              #size of axis
    plt.minorticks_on()             #enabling minor grid lines
    plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        
                                #which both to plot major and minor grid lines
    plt.grid(which = 'major')
    plt.savefig(Savename + '.png', format='png', bbox_inches='tight')
    plt.show() 
    
    
    
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
            y_label = 'y', Color = 'b', save_name = '', post_title = ' ',
            Fit_type = 1):    
    '''
    Function to do and compute some variables relevant to the PSO (pseudo second 
    order) kinetic model. THis model comes from
    d(Q(t))/dt = K * (Q_e - Q(t))**2,
    whose solution is:
        Q(t) = Qe**2 * K * t / (1+Qe*K*t)
    where Q_e = Q(t ==> \infty).
    This funciton could do either the non-linear fit (that eq), or a linear one.
    THe solution of that eq can be casted in linear form:
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
        .Fit type: number to indnicate fit type. 1 = linear (default). 0 for non linear
         else will return NaN
    
    *Outputs
        .df series with all the relevant info, from the fit and computed 
        quantities, errors (quadratic propagation) included The input units 
        define those units!! Remember saltpepper!
    
    
    #------------- To Do -----------
        -Non linear fit version, to be developed still, not work well
    
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


    if Fit_type==1:                     #Linear fit type
        #
        t__Q = t / Q          #t/Q(t) for S
        Delta_t__Q = np.abs(t__Q) * np.sqrt((delta_Q / Q )**2 + 
                                        (delta_t /t )**2 )  #error, unused!!
        #Fit
        fit = Fits.LinearRegression(t, t__Q, delta_t, Delta_t__Q,
                                   x_label = x_label, y_label = y_label, 
                                   Color = Color, 
                save_name = folder_name +'/' + save_name, post_title = post_title)       
                            #Fit (i dont use npo variable, fit variable)
        ################ 3) Model parameters ################
        '''
        From that I can also get Qe and K easy: y = ax + b;
         a = 1/Qe ==> Qe = 1/a
         b = 1/KQe**2 = a**2/K == > K = a**2 /b
         '''
        fit['Q_e'] = 1 / fit['a']         #Qe = 1/a, y= ax + b
        fit['\Delta(Q_e)'] = fit['\Delta(a)'] /fit['a']**2     #Delta(Qe)
        fit['%rsd Q_e'] = fit['\Delta(Q_e)']/np.abs(fit['Q_e']) * 100
        fit['K'] = fit['a']**2 /fit['b']         #K = 1/b * Qe**2 = a**2/b
        fit['\Delta(K)'] = np.abs(fit['K']) * np.sqrt( 2*(fit['\Delta(a)'] / fit[
        'a'] )**2 + (fit['\Delta(b)'] / fit['b'])**2 )  
                            #Delta(K) np.abs() so its always >0
        fit['%rsd K'] = fit['\Delta(K)']/fit['K']*100
        
    elif Fit_type ==0:                                    #Non linear tyoe

        def PSO_eq(t, Qe, K):       #t the independent variable, other parameterets
            return Qe**2*K*t/(1+ Qe*K*t)
    
        #Initial paratermers for the guess
        Qe0 = np.nanmax(Q)
        K0 = 1 / (Qe0 * (t[np.argmax(Q > Qe0/2)] if np.any(Q > Qe0/2) else t.mean()))
        p0 = [Qe0, K0]  
        popt, pcov = curve_fit(PSO_eq, t, Q, p0=p0, maxfev=10000)
        Qe,K = popt
        perr = np.sqrt(np.diag(pcov))  # errors
        dQe, dK = perr
        
        # --- Predicted values
        Qe_pred = PSO_eq(t, *popt)

        # --- Goodness of fit (R^2)
        ss_res = np.sum((Q - Qe_pred) ** 2)
        ss_tot = np.sum((Q - np.mean(Q)) ** 2)
        R2 = 1 - (ss_res / ss_tot)
        
        #Debug
        print(f"ss_res={ss_res}, ss_tot={ss_tot}, mean(Q)={np.mean(Q)}")
        #If ss_tot is ~0 or ss_res negative, you’ve found the culprit.
        
        #Storage
        fit = pd.Series({ 'Q_e': Qe, '\Delta(Q_e)' : dQe, '%rsd Qe': dQe/np.abs(Qe)*100,
                         'K': K, '\Delta(K)': dK, '%rsd K': dK/K * 100,
                         'r': R2})
        
    else:                                                        #Wong case
        print('Wrong fit type given. Accepted 1(Linear) and 0 (non linear)')
        fit = np.nan
        print('Returning nan!')
        
    ########## 4) Return ###########
    
    return fit
    

    #---------- Non linear fit
    # from scipy.optimize import curve_fit

    # def PSO_model(t, Qe, K):
    #     return (K * Qe**2 * t) / (1 + K * Qe * t)

    # popt, pcov = curve_fit(PSO_model, Dict_Time_S['< >'].drop(['Data 1', 'Data 13']).values,
    #         Dict_Qe_S['< >'].loc['U(LR)' ].drop('Data 12').values, 
    #         p0=[max(Dict_Qe_S['< >'].loc['U(LR)' ].drop('Data 12')), 0.01])
    # Qe, K = popt
    # perr = np.sqrt(np.diag(pcov))   # uncertainties


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
    ax.errorbar(t, Q, delta_Q, delta_t, 'o', color = Color, markersize = Markersize, label = 'Data')
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
def Fre_fit(Ce, Qe, delta_Ce=0, delta_Qe =0, folder_name = 'Fits',
            x_label = 'log($C_e [ng/g]$)', y_label = 'log($Q_e [ng/g_{be}]$)',
            Color = 'b', save_name = '', post_title = ' ', npo=100,
            Fit_type = 1):   
    '''
    Function to do and compute some variables relevant to the Freundlich fit
    of an adsorption isotherm Q_e = f (C_e), the sorbed quantity as a function
    of the equilibrium concentration. Its equation is
    
        Q_e = K_F * C_e**(1/n),           K_F, n constants
    Being 1/n the heterogeneity factor. n= 1 ==> linear sorption. That can also be
    linearized (log 10 used normally):
                    loq Q_e = log K_F + 1/n log C_e
    
    *n>1 ==> 1/n <1 ==> 
        .The adsorbent becomes gradually saturated, but adsorption still occurs 
        efficiently even as concentration rises.
        .It reflects a heterogeneous surface with high-affinity sites that get 
                occupied first.
    *n<1 ==> 1/n >1 ==>
        .Adsorption becomes more difficult as concentration rises.
        .It suggests possible cooperative or multilayer effects — not a simple 
        surface sorption process.
        
    The lineal fit of that is trivial:  y= ax + b
                y= log Q_e
                x= Log C_e
                a = 1/n
                b= log K_F
    (easier than if I use 1/n as constant!)
    
    n is adimensional, but K_F has unit, which depends on Ce; Kf = Qe/Ce**n
        i) Ce in ppb, ng/g. Then K_F is in ng^(1-n)*g_tot^n/(g_be)
        ii) Ce in M = mol/L. Then K_F is in L^n/(mol^(n-1)*Kg) = 
            L^n mol^(1-n)/kg
    
    Note I add + (LR) to the column name in the fit serie!! Watch out, maybe 
    you need to modify it in the future??????
    
    *Inputs
        .Ce, Qe: df series containing C_e and Q_e data. Expected the averaged values
            . Must have same index as the df columns in order to plot them!!
        .delta_Ce, delta_Qe: df with their uncertainties. Default value = 0, 
        since I do not use  them!
        .x_label, y_label= x and y label, for the plot. Default value: 
            'log($C_e [ng/g]$)' and 
                    log($Q_e [ng/g_{be}]$)' respectively!
        .post_title = '' : title to add after 'Linear fit '
        .save_name = filename of the fit plot, if it wants to be save. 
        Default value = '' ==> no saving.
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
    First the folder to store the plots will be created. IN the main folder a
    subfolder with the relevant elements, to be given, will be created
    '''
    
    path_bar_pl = os.getcwd() + '/' + folder_name + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)


    
    if Fit_type ==1:                #linear type
    #
        logCe = np.log10(Ce)                        #cal of the log
        logQe = np.log10(Qe)
        delta_logCe = delta_Ce / (np.abs(Ce) * np.log(10))   #error of the log10!!
        delta_logQe = delta_Qe / (np.abs(Qe) * np.log(10))
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
    From that I can also get the constants; applying log10 == log, no ln!!!, 
    since log properties are applicable regardless of the base:
        y= ax + b : loq Qe = 1/n log Ce + log KF
                    y= log Q_e
                    x= Log C_e
                    a = 1/n
                    b= log K_F ==> K_F = 10**b
                    delta_K_F = delta_b * ln(10) * K_F
                    delta_n = delta_a/a**2
                    (partial log(Kf) / partial Kf = 1/ (KF ln(10)) )
    
    And the units?
    n is adimensional, ofc
    K_F not, since Q_e = K_F * C_e**1/n ==> K_F = Q_e * C_e*n ==>
    K_F is ng/g_be * (ng/g_tot)**n = ng**n+1/(g_be * g_tot**n)
        '''
        fit['n'] = 1/fit['a']         
        fit['\Delta(n)'] = fit['\Delta(a)']/ fit['a'] **2
        fit['K[L^n/(kg*mol^{n-1})]'] = 10**fit['b']        
        fit['\Delta(K[L^n/(kg*mol^{n-1})])'] = fit['K[L^n/(kg*mol^{n-1})]'] *fit[
            '\Delta(b)'] * np.log(10)
    #Return the %rsd will be useful to know how relevants are the aprameters!
        fit['%rsd K'] = fit['\Delta(K[L^n/(kg*mol^{n-1})])'] / fit['K[L^n/(kg*mol^{n-1})]'
                                                                   ] * 100
        fit['%rsd n'] = fit['\Delta(n)'] / fit['n'] * 100
        '''
        K_F error may be higher than you expected, but since you work with lgoaritm
        this could happen!
        '''
    elif Fit_type == 0:                         #non linear type
        #
        def Fre_fit_eq(C,K, n):         #fre fit eq, non linear (original)
            return K*C**(1/n)
        # --- Initial guess: K ~ max (Qe)/max(Ce), n  ~ 1
        p0 = [np.max(Qe)/np.max(Ce), 1]  
        popt, pcov = curve_fit(Fre_fit_eq, Ce, Qe, p0=p0, maxfev=10000)
        K,n = popt
        perr = np.sqrt(np.diag(pcov))  # errors
        dK, dn = perr

        # --- Predicted values
        Qe_pred = Fre_fit_eq(Ce, *popt)

        # --- Goodness of fit (R^2)
        ss_res = np.sum((Qe - Qe_pred) ** 2)
        ss_tot = np.sum((Qe - np.mean(Qe)) ** 2)
        R2 = 1 - (ss_res / ss_tot)
        
        fit = pd.Series({
        'K[L^n/(kg*mol^{n-1})]': K,
        '\Delta(K[L^n/(kg*mol^{n-1})])': dK,
        'n' : n,
        '\Delta(n)' : dn,
        'r': R2,
        '%rsd K' : dK/K*100, '%rsd n': dn/n*100})
        
        # --- Plot
        Ce_range = np.linspace(min(Ce), max(Ce), 200)
        plt.figure(figsize=(11,8))
        plt.errorbar(Ce, Qe, xerr=delta_Ce, yerr=delta_Qe,
               fmt='o', color=Color, label='Data', markersize = Markersize)
        plt.plot(Ce_range, Fre_fit_eq(Ce_range, *popt),
           '--', color=Color, label=f'Fit (R²={R2:.3f})')
        plt.xlabel('$C_e$ [M]', fontsize = Font)
        plt.ylabel('$Q_e$ [mol/kg$_{be}$]', fontsize = Font)
        plt.title(f'{post_title}', fontsize=22)
        plt.grid(True)
        plt.tick_params(axis='both', labelsize=Font)
        plt.legend(fontsize = Font)
        #plt.tight_layout()
        plt.savefig(folder_name +'/' + save_name + '.png', format = 'png',
                    bbox_inches='tight')
        plt.show()
        #
        #
    else:       #wrong type given
        print('Wrong fit type given. Accepted 1(Linear) and 0 (non linear)')
        fit = np.nan
        print('Returning nan!')
        
    ########## 4) Return ###########
    
    return fit



#%% ######### 2.4) Langmuir isot fit #############################
###################################################

def Lang_fit(Ce, Qe, delta_Ce=0, delta_Qe =0, Fit_type = 1,
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
                                   uncertainties)
        
        Here you plot Ce (x) vs Ce/Qe (y)
        2) Q_e/C_e = -K_L* Q_e + Q_max * K_L
        Here you plot Qe (x) vs Qe/Ce (y)
    I will do these function so that you can do one of the 2 linearizations, or
    the normal version.
    
    *Inputs
        .Ce, Qe: df series containing C_e and Q_e data. Expected the averaged
                values. Must have same index as the df columns in order to 
                plot them!!
        .delta_Ce, delta_Qe: df with their uncertainties. Default value = 0, 
        since I do not use them!
        .x_label, y_label= x and y label, for the plot. Default value: 
            'log($C_e [ng/g]$)' and log($Q_e [ng/g_{be}]$)' respectively!
        .post_title = '' : title to add after 'Linear fit '
        .save_name = filename of the fit plot, if it wants to be save. 
            Default value = '' ==> no saving. this variable is followed 
            by .png for savinf
        .Color = 'b': color for the plot
        .Folder_name: folder name, where to store the fit plots
        .npo=100: number of points for the fit plot
        .Fit_type: linearization type. Default: 1, meaning linearization 1, 2 for
        lin 2. other that 1 or 2 for non-linear lin!
    
    
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
    delta_Ce_Qe = np.abs(Ce_Qe) * np.sqrt((delta_Ce / Ce)**2 + 
                                  (delta_Qe / Qe)**2 )
    Qe_Ce = Qe / Ce
    delta_Qe_Ce = np.abs(Qe_Ce) * np.sqrt((delta_Ce / Ce)**2 + 
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
        print("#------------------------------- ")
        print("Linearization 1 done, Ce/Qe vs Ce! ")
        print("#------------------------------- \n")
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
        fit["Q_max[mol/kg_be]"] = 1/ fit["a"]
        fit["\Delta(Q_max[mol/kg_be])"] = fit['\Delta(a)']/fit['a']**2
        fit["K_L[L/kg]"] = 1/ (fit["b"] * fit["Q_max[mol/kg_be]"] )
        fit["/Delta(K_L[L/kg])"] = fit["K_L[L/kg]"] * np.sqrt(
             (fit["\Delta(Q_max[mol/kg_be])"]/fit["Q_max[mol/kg_be]"])**2 + 
             (fit['\Delta(b)']/fit["b"])**2 )
        fit['%rsd Q_max'] = fit["\Delta(Q_max[mol/kg_be])"]/fit["Q_max[mol/kg_be]"] * 100
        fit['%rsd K_L'] = 100 * fit["/Delta(K_L[L/kg])"] / fit["K_L[L/kg]"]
    elif Fit_type ==2:       #Linearization 2!
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
        print("#-------------------------------")
        print("Linearization 2 done, Qe/Ce vs Qe! ")
        print("#----------------------------------\n")
    
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
        fit["K_L[L/kg]"] = -fit["a"]
        fit["Q_max[mol/kg_be]"] = fit["b"] / fit["K_L[L/kg]"]
        
        fit["/Delta(K_L[L/kg])"] = fit['\Delta(a)']        
        fit["\Delta(Q_max[mol/kg_be])"] = fit["Q_max[mol/kg_be]"] * np.sqrt(
             (fit["/Delta(K_L[L/kg])"]/fit["K_L[L/kg]"])**2 + 
             (fit['\Delta(b)']/fit["b"])**2 )

        fit['%rsd Q_max'] = fit["\Delta(Q_max[mol/kg_be])"]/fit["Q_max[mol/kg_be]"] * 100
        fit['%rsd K_L'] = 100 * fit["/Delta(K_L[L/kg])"] / fit["K_L[L/kg]"]   

    elif Fit_type == 0:           #no linear fit, normal fit!
        print("#------------------------------- ")
        print("Non linear fit will be done! (or tried at least xD ")
        print('Beware with r, might be misleading. Parameters must be checked!')
        print("#----------------------------------\n")
        #
        def Lang_fit_eq(Ce, Qmax, K):      #equaition fo langmuir fit
            #Q_e = Q_max * K_L/(1+K_L* C_e) * C_e
            return Qmax * K * Ce /(1+K*Ce) 
        # --- Initial guess: Qmax ~ max(Qe), KL ~ 1/mean(Ce)
        p0 = [np.max(Qe), 1/np.mean(Ce)]  
        popt, pcov = curve_fit(Lang_fit_eq, Ce, Qe, p0=p0, maxfev=5000)
        Qmax, KL = popt
        perr = np.sqrt(np.diag(pcov))  # errors
        dQmax, dKL = perr

        # --- Predicted values
        Qe_pred = Lang_fit_eq(Ce, *popt)

        # --- Goodness of fit (R^2)
        ss_res = np.sum((Qe - Qe_pred) ** 2)
        ss_tot = np.sum((Qe - np.mean(Qe)) ** 2)
        R2 = 1 - (ss_res / ss_tot)
        
        
        fit = pd.Series({
        'Q_max[mol/kg_be]': Qmax,
        '\Delta(Q_max[mol/kg_be])': dQmax,
        'K_L[L/kg]': KL,
        '/Delta(K_L[L/kg])': dKL,
        'r': R2,
        '%rsd K_L' : dKL/KL*100, '%rsd Q_max': dQmax/Qmax*100})
        
        # --- Plot
        Ce_range = np.linspace(min(Ce), max(Ce), 200)
        plt.figure(figsize=(11,8))
        plt.errorbar(Ce, Qe, xerr=delta_Ce, yerr=delta_Qe,
               fmt='o', color=Color, label='Data', markersize = Markersize)
        plt.plot(Ce_range, Lang_fit_eq(Ce_range, *popt),
           '--', color=Color, label=f'Fit (R²={R2:.3f})')
        plt.xlabel('$C_e$ [M]', fontsize = Font)
        plt.ylabel('$Q_e$ [mol/kg$_{be}$]', fontsize = Font)
        plt.title(f'{post_title}', fontsize=22)
        plt.grid(True)
        plt.tick_params(axis='both', labelsize=Font)
        plt.legend(fontsize = Font)
        #plt.tight_layout()
        plt.savefig(folder_name +'/' + save_name + '.png', format = 'png',
                    bbox_inches='tight')
        plt.show()
        
    else:
        print('Wrong fit type introduced! Returning fit = np.nan ')
        print('Values: 1 for lin 1, 2 for lin 2, 0 for non linear version')
        print('-------------------------------\n')
        fit = np.nan
    
    ########## 5) Return ###########
    
    return fit



#%% ######### 2.5) D-R iso fit #############################
###################################################
def D_R_fit(Ce, Qe, delta_Ce=0, delta_Qe =0, T = 293.15, delta_T = .1, 
            folder_name = 'Fits', x_label = 'log($C_e [ng/g]$)', 
            y_label = 'log($Q_e [ng/g_{be}]$)',
            Color = 'b', save_name = '', post_title = ' '):   
    '''
    Function to do and compute the D-R fit model. Its equation is:
        
    Qe = Qs * exp(-beta eps**2),
    eps = RT*ln(Ct/Ce),
    
    Ct = constant, usually taken as 1M, or 1g/L, depending on your data
    [Wang2020] said that it could be the solubility limit. Watch out!

    Linearizing it:
        log(Qe) = log(Qs) - beta*eps**2
    
    Note that lin fit we will do. The units of Qe will define the units of Qs!

    Note that you could compute the mean free energy F as
    F = 1/sqrt(-2beta)
    And from that now the nature of the process:
        F<8KJ/mol ==> physisorption (van der Waals, weak electrostatic)
        8<F<16 ==> ion exchangee / weak chemisorption
        F> 16kJ/mol ==> strong chemisorption
    [Chabani2006]
    Note I add + (LR) to the column name in the fit serie!! Watch out, maybe 
    you need to modify it in the future??????
    
    *Inputs
        .Ce, Qe: df series containing C_e and Q_e data. Expected the averaged values
            . Must have same index as the df columns in order to plot them!!
        .delta_Ce, delta_Qe: df with their uncertainties. Default value = 0, 
        since I do not use  them!
        .T: temperature of the experiment [K], or of each sample. Default: 
                293.15 (20°)
        .delta_T: uncertainty of the temperature [K]. Default: 0.1
        .x_label, y_label= x and y label, for the plot. Default value: 
            'log($C_e [ng/g]$)' and 
                    log($Q_e [ng/g_{be}]$)' respectively!
        .post_title = '' : title to add after 'Linear fit '
        .save_name = filename of the fit plot, if it wants to be save. 
        Default value = '' ==> no saving.
                    this variable is followed by .png for saving
        .Color = 'b': color for the plot
        .Folder_name: folder name, where to store the fit plots
    
    
    *Outputs
        .df series with all the relevant info, from the fit and computed quantities, 
        errors (quadratic propagation) included
            The input units define those units!! Remember saltpepper!
    
    
    '''    
    ############# 0.1) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder a
    subfolder with the relevant elements, to be given, will be created
    '''
    
    path_bar_pl = os.getcwd() + '/' + folder_name + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)

    
    ############## 1) Calcs #################
    R = 8.31446261815324    #J/(K*mol) molar gas constant
    
    #I need to compute the logarithms!
    logQe = np.log(Qe)    
    delta_logQe = delta_Qe / np.abs(Qe)       #error of the natural log!
    
    eps = R*T*np.log(1/Ce)                      #J/mol
    delta_eps = eps* np.sqrt((delta_T/T)**2 + (Ce**4))
    
    #Elevating it by square:
    eps2 = eps*eps          #J2/mol2
    delta_eps2 = np.sqrt(2) * eps * delta_eps
    
        #delta(log(1/Ce))/(1/Ce) = 1/(1/Ce**2)=Ce**2
    ############# 2)Fit ######################
    
    fit = Fits.LinearRegression(eps2, logQe, delta_eps2, delta_logQe,
                                   x_label = x_label, y_label = y_label, 
                                   x_legend = '$\epsilon$^2', y_legend = 'log($Q_e$)',
                                   Color = Color, 
                                   save_name = folder_name +'/' + save_name, 
                                   post_title = post_title)       
                            #Fit (i dont use npo variable, fit variable)
                #note that for the legnd I delete the units!!
    
    
    ################ 3) Model parameters ################
    '''
    From that I can also get the constants; applying log10 == log, no ln!!!, 
    since log properties are applicable regardless of the base:
        y= ax + b : loq Qe = log Qs - beta*eps**2
                    y= log Q_e
                    x= eps**2
                    a = -beta
                    b= log Qs ==> Q_s = e**b
                    delta_beta = delta_a
                    delta_b = delta_Qs/Q_s ==Y delta_Qs = Q_s * delta_b
    '''
    fit['beta[mol2/J2]'] = -fit['a']         
    fit['\Delta(beta[mol2/J2])'] = fit['\Delta(a)'] 
    fit['Q_s'] = 10**fit['b']        
    fit['\Delta(Q_s)'] = fit['Q_s'] *fit['\Delta(b)']
    
    #Lets finally compute the mean free energy F, to see the nature of the process
    F = 1/np.sqrt(2*fit['beta[mol2/J2]'])   #J/mol
    delta_F = F *   fit['\Delta(beta[mol2/J2])'] / (2*fit['beta[mol2/J2]'])      
        #After some uncertainty calcs I derived that formula
    
    #Now we can print which type of sorption is, based on F:
    print('##################################################')
    print(f' F = {F:.1e} +- {delta_F:.1e}' + 'J/mol')
    if F*10**-3 < 8: #physi
        print('F < 8KJ/mol, indicating it is physisorption ')
    elif (F*10**-3 < 8 and F*10**-3 > 16): #chemi
        print('F> 16kJ/mol, indicating it is chemisorption')
    else : #F*10**-3 > 16     #WTF is that?
        print('8KJ/mol < F < 16KJ/mol, Ion exchange/soft chemisorption')
    
    print('##################################################')
    
    #Finally we store them
    fit['F[J/mol]'] = F    
    fit['\Delta(F[J/mol])'] = delta_F
    
    ########## 4) Return ###########
    
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

def Read_XRD_F130 (name, t_paso = 10, Skip_rows = 266, Compute_d = 0, 
                   DosTheta_inter = [4,6], XRD_Type = 'Cold'):
    '''
    Function that reads the .ras file from the XRD for both:
        . F130, Cold (Olaf Walter)
        . A232, Active, (Olaf also xD)
    returning a df with the relevant info. It will also plot it and save it, 
    in the same folder, and if desired, also compute the interlaminar space from
    the 001 peak (aorund 6°)
    
    The 2 XRD gives differnt, output files, so will be treated differently:
        .Cold: .ras file contain lot of text at the beginning (1st 266 lines), then 
            3 variables:
            2Theta, Counts, Unknown
            The unknown seem to be always 1, so I delete it.
        .Active: .uxd file, containing text at ebginning and after 2 columns
                2Theta PSD 
    
    Note that some of the text I skip contian relevant info. Eg:
            *MEAS_SCAN_START_TIME "04/10/24 15:21:09" (near the last lines of text)
    
    *Inputs:
        .name: name of the file. Ej: 'file.ras'
        .t_paso: time per step, in seconds. Ej: 10 [s]. Thi assumes the operation 
        mode in step, the usual
        .Skip_rows: number of rows of the .ras file to skip. Default: 266. Spotted 
        from opening the file. For Hot, 62
        .Compute_d : boolean to indicate if you can to compute the interlaminar
        space distance from the 001 peka or not. Default: 0 (no)
        .2Theta_inter: array indicating the interval of the 001 peak, in case that
        Compute_d = 1. Default: [4,6]
        .XRD_type: string indicating the XRD device: 'Cold' or 'Hot'. Default:
            'Cold' (F130)
        
    *Output
        .df with the 2Theta, Counts and cps
    
    '''
    
    ##### 1. Reading #########
    #Reading was not trivial, but I came up with the way of doing it:
        

    if XRD_Type== 'Cold':                           #Cold XRD
        aux = pd.read_csv(name, skiprows = Skip_rows, sep = ' ', 
                      names = ['2Theta[°]', 'Counts', 'npi'], encoding='latin')
        '''That worked, but 2 last columns are text, that I can delete easily. 
        Also I have 3 column,  being the 3rd always 1, so I will also delte it:'''
    
        df = aux.iloc[:-2,:-1]        #Removing 3rd column, and last 2 rows
    
    elif XRD_Type == 'Hot':                 #hot XRD
        df = pd.read_csv('XRD/DL_Bentonite_raw.uxd', skiprows = Skip_rows, sep = '      ', 
                          names = ['2Theta[°]', 'Counts'], encoding='latin')
                    #separato spotted from the .udx file
    else:       
        print('Wrong XRD type, accepting only Hot or Cold !')
        print('Retunrin nan, pssiibly error in the function !')
        print('------------------------------------\n')
        df = np.nan
        
    #Since the 2theta data is not numeric because of these 2 rows, I need to make
    #them numeric:
    df = df.apply(pd.to_numeric)            #conversion to numeric in case its not
    
    #The cps are easy to get, since I define the time measuring each step:
        
    df['CPS'] = df['Counts']/t_paso         #Computing Counts Per Secons (CPS)


    '''
    I will also return the cps nromalized, min value 0 and max 1. This can be 
    accomplished by appling the operation:
            (x-min(x)) /(max(x)-min(x))         para todo x
    '''
    
    df['CPS_norm'] = ( df['CPS'] -df['CPS'].min() ) / (df['CPS'].max() - df['CPS'].min())
                        #normalized cps, min value 0, and max 1
    
    
    ##Prints useful for command line sytling (also for interlaminar space)
    print('\n##############')
    print('Currently working with:\n') 
    print( name[:-4] + '\n')
    
    
    
    ########### 2. interl space calc ##############
    if Compute_d:       #if true, compute it
        print('\n##############')    
        print('Beware woth the initial interval, you might need to modify it'+
        'since the first guess may not be okay/good enough and could give error \n')    
        try:    #Try the fit
            Fit = XRD_Get_interl_sp (df, DosTheta_inter)
        #
        #Now saving the results in the previous df
            df[Fit.columns] = Fit
        except RuntimeError:        
                #IF error because Fit does not get a fit xD (is RuntimeError)
            print('\n##############') 
            print('Impossible to do the fit with the given 2Theta interval, give'+
                  'a better one, puto!!!!!!!')
            print('The legend of the plot will be remarking that xD')
            print('\n##############')
        
        #
        #Now in the plot we will include a textbox with the d value:
        plt.figure(figsize=(11, 8))  # width, heigh 6.4*4.8 inches by default (11,8)
        # I need to enlarge since because of the title is big, I think
        plt.title('XRD: ' + name[:-4], fontsize=22,
                  wrap=True, loc='center')  # title
        try:
            plt.plot(df['2Theta[°]'], df['CPS_norm'], 
                 label='$d_{001}$ =(' + f'{df['d[A]'].iloc[0]:.2f}' + '+-' +
                 f'{df['\Delta(d)[A]'].iloc[0]:.2f}' + ')Angs')
        except KeyError:    
            #If fit did not work, d(A) does not exist, to previous label will not work
            plt.plot(df['2Theta[°]'], df['CPS_norm'], 
                 label='No fit possible bro xD')
        plt.xlabel("2"r"$\theta $ (°) ", fontsize= Font)  # xlabel
        plt.ylabel('cps', fontsize= Font)
        plt.tick_params(axis='both', labelsize= Font)  # size of axis
        plt.minorticks_on()             #enabling minor grid lines
        plt.grid(which = 'minor', linestyle=':', linewidth=0.5)  
                    #which both to plot major and minor grid lines
        plt.grid(which = 'major')
        plt.legend(fontsize= Font)
        plt.savefig( name[:-4]+ '.png', format='png',
                    bbox_inches='tight')  
                    # To save plot, same name as file, change extensio
        plt.show()    
        #
    else:                           #False, so I did not compute the distance
        'False, no distance computed. Then I only do the plot'
    #
        plt.figure(figsize=(11, 8))  # width, heigh 6.4*4.8 inches by default (11,8)
        # I need to enlarge since because of the title is big, I think
        plt.title('XRD: ' + name[:-4], fontsize=22,
              wrap=True, loc='center')  # title
        plt.plot(df['2Theta[°]'], df['CPS_norm'], label='data')
        plt.xlabel("2"r"$\theta $ (°) ", fontsize= Font)  # xlabel
        plt.ylabel('cps', fontsize= Font)
        plt.tick_params(axis='both', labelsize= Font)  # size of axis
        plt.minorticks_on()             #enabling minor grid lines
        plt.grid(which = 'minor', linestyle=':', linewidth=0.5)  
                #which both to plot major and minor grid lines
        plt.grid(which = 'major')
                #plt.legend(fontsize= Font)
        plt.savefig( name[:-4]+ '.png', format='png',
                bbox_inches='tight')  # To save plot, same name as file, change extensio
        plt.show()
    
    
    #Finally, return the df:
    return df


#%% ######### 5.3) XRD, Get interlaminar space bentonites #################### 
######################################################

def XRD_Get_interl_sp (XRD_df, DosTheta_inter, Kalpha = 1.5401):
    '''
    Function that will get the interlaminar space of the bentonites, based on the
    first basal refelction 001, which will satisfy:
            lambda = 2d sin (theta)    n = 1
            
    For this, we need:
        1. Perform the gaussian fit of the peak (Fits module used!)
        2. Compute the interlaminar space d
    
    The .ras file contain lot of text at the beginning (1st 266 lines), then 3 
    variables:  2Theta, Counts, Unknown
        The unknown seem to be always 1, so I delete it.
        
    Some of the text I skip contian relevant info. Eg:
            *MEAS_SCAN_START_TIME "04/10/24 15:21:09" (near the last lines of text)
    
    *Inputs:
        .XRD_df: df containing the info from the XRD. Important the names, must have
            '2Theta[°]' and 'CPS_norm' (I refer to those
        .DosTheta_inter: array with min and maximum values of the 2Theta variable 
        from the XRD df, to do the fit. Ej: [5,30]
        .kalpha: Wavelength of the Kalpha1 radiation of the Cu, the element of
            the XRD = 1.5401Angs
        .name: name of the file. Ej: 'file.ras'

    *Output
        .df with the interlaminar space!
    
    '''
    
    
    ####### 1. Gaussian fit #########
    #We need to get the desired interval to perform the fit, eye spotted from 
    #the XRD diagram
    
    inter_x =XRD_df['2Theta[°]'].loc[ (XRD_df['2Theta[°]'] < DosTheta_inter[1]) & 
                                        (XRD_df['2Theta[°]'] > DosTheta_inter[0]) ]
        #getting desired interval
    #Now I need the cps values (y) of those x values. Since the index is preserve, 
    #I could use it
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
    print('Interlaminar space of  (see above):' + f'{aux["d[A]"] : .3f}' + ' +-' + f'{aux["\Delta(d)[A]"] : .3f}')
    print('#######################\n ')
    
    #Finally, return d and its error in a df, recycling the one for the fit
    Fit['d[A]'] = d
    Fit['\Delta(d)[A]'] = Delta_d
    
    return Fit


#%% ######### 6) FTIR, read and plot #################### 
######################################################

def Read_FTIR (name, Type = 'A', Plot = 'A', Sep = ','):
    '''
    Function that reads the .dpt file from the FTIR in F130 (Olaf Walter), returning a df
    with the relevant info. It will also plot it, and save it!
    In case the data its in absorbance, it will convert it to transmitance!
    

    *Inputs:
        .name: name of the file. Ej: 'file.dpt'. If in a folder: 'Folder/name.dpt'
        .Type: string indicating whether the data its Absorbance (A) or transmitance (T).
            data. Default: 'A'
        .Sep: string indicating the separator. Default: a comma, ','. Could also be
            '/t'  (a tab)
        .Plot: string to define what shall we plot, the transmitance (T) or 
            absorbance (A). Default: 'A'
        
    *Output
        .df with the 1/lambda and the transmitance!
    
    '''
    
    ##### 1. Reading #########
    '''
    Reading the data is trivial. But of course it depends on whether the file
    contain the absorbance or the transmitance. The relation between them are:
        A = -log (T)  (log10) ==> T = 10**-A
    '''
    
    if Type == 'A': #It is absorbance
        df = pd.read_csv(name, sep = Sep, names = ['1/lambda[cm-1]','Absorbance'])
        df['Transmitance'] = 10**(-df['Absorbance'])
        #Now we put them in % scale
        df['Absorbance[%]'] =  df['Absorbance']*100
        df['Transmitance[%]'] = df['Transmitance']*100
    else: #It is transmitance
        df = pd.read_csv(name, sep = Sep, names = ['1/lambda[cm-1]','Transmitance'])
        df['Absorbance'] = -np.log10(df['Transmitance'])
        df['Absorbance[%]'] =  df['Absorbance']*100
        df['Transmitance[%]'] =  df['Transmitance']*100
        

    ##Prints useful for command line sytling (also for interlaminar space)
    print('\n##############')
    print('Currently working with:\n') 
    print( name[:-4])
    
    
    ########### 2. Plotting ##
    plt.figure(figsize=(11, 8))  # width, heigh 6.4*4.8 inches by default (11,8)
    # I need to enlarge since because of the title is big, I think
    plt.title('FTIR: ' + name[:-4], fontsize=22,
              wrap=True, loc='center')  # title
    #plt.plot(df['1/lambda[cm-1]'], df['Transmitance[%]'], label='data')
    if Plot== 'A':           #plotting Abs
        plt.plot(df['1/lambda[cm-1]'], df['Absorbance[%]'], label='data')
        plt.ylabel('Absorbance [%]', fontsize= Font)
    else:   #plotting transm
        plt.plot(df['1/lambda[cm-1]'], df['Transmitance[%]'], label='data')
        plt.ylabel('Transmitance [%]', fontsize= Font)
    #
    plt.xlabel(" $1/\lambda[cm^{-1}]$", fontsize= Font)  # ylabe
    plt.gca().invert_xaxis()            #to invert x axis (1st higher values, then lower)
    plt.tick_params(axis='both', labelsize= Font)  # size of axis
    plt.minorticks_on()             #enabling minor grid lines
    plt.grid(which = 'minor', linestyle=':', linewidth=0.5)       
                #which both to plot major and minor grid lines
    plt.grid(which = 'major')
    #plt.legend(fontsize= Font)
    plt.savefig( name[:-4]+ '.png', format='png',
                bbox_inches='tight')  # To save plot, same name as file, change extensio
    plt.show()
    
    
    #Finally, return the df:
    return df


