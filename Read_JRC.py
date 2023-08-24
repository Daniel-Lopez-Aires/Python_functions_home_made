# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 10:29:58 2023

@author: lopedan

This script will contain all the functions I create to read files from experimental measurements,
say TGA, XRD, etc (the ones that needed ofc)
"""

#%%######################################
############## 0) General packages ###########
######################################


import matplotlib.pyplot as plt  #for simplicity, to not write matplotlib.pyplot
        #everytime we want to plot something
#from scipy.stats import norm               ##norm.fit() fit to gaussian
import numpy as np
    #np contain linspaces as np.linspace(a,b,N)
import pandas as pd
import os, sys                   #to import functions from other folders!!
sys.path.insert(0, '//net1.cec.eu.int/jrc-services/KRU-Users/lopedan/Desktop/PhD_Residuos_nucleares/Python/Functions')   
                                    #path where I have the functions
import Fits, Peak_analyis_spectra
import time as tr                                #to measure the running time

#############################################################



#%%######################################
########### 1.1) ICPMS excel reader #############
#####################################

def Read_ICPMS_excel (excel_name, cps_sheet_name = 'To_read', return_debug = False):
    '''
    Function that will read the excel file from ICPMS and will return a df with the relevant
    information, for easier handling /plotting. Note the excel should be a bit preprocessed:
            
        1) Clean sheet where only the relevant data(cps, ppb, whatever) is, to load it. 
                .You can remove ICPMS blanks (std 0, etc). THe Isotopes column in COlumn A in excel
                THe 1st isotope, Co59(LR) in row 7 in excel). Sample names in row 2 in excel
            
    
    You could use this function to get the raw data (output from ICPMS) or to correct 
    them for the ICPMS dilution factor.     

    Note sometimes some random NaN data from excel can be added. Easy solution, go to the excel sheet,
    and delete those rows/columns. No clue why this happens, but that solves it (:   
                                                                                 
    Recommended to delete the "wash" sample, will only bring problems in analyis xD

    
    *Inputs:
        .excel_name: string with the name of the excel, with the .xlsx. note if you select the file
        on the FIle viewer and do copy paste, the name will be there. But dont forget the '', it must
        be an string! Eg: 'Excel.xlsx'
        .cps_sheet_name: string with the name of the sheet with the data to read 
            future maybe also concentration values?). Default value: 'To_read' 
        (from acid vs no acid test)
        .return_debug: if you want to get some extra df for debug 
        (raw data, without cleaning, so like the excel). Default value = False

        
    *Outputs:
        .several df with the cps/%rsd or whatver it is reading. Depending whether you want the debug you may 
            have 1 or 2 outputs. THe isotopes are the index of the df, so the columns are the sample data!
            the column names are the sample names
            
            
    Note that if you have N outputs, if you want to obtain a variable per output, 
    in the script I should call X variables, like:
            a, b, ..n = Read_ICPMS_excel(name)
    If you write less, say 1, 2, some variable will contain more data, in a dictionary
            
            
    ######## TO DO ########
        . Arbitrary inputs (type *arb_arguments, with other stuff like, variable1, *aribtrary_arg)
            so you choose if you want rsd or no, to save some time?
            
    #######################
        '''
    
    
    ########### 1) Raw load ###########
    '''
    Can be done easily with pandas. Since th excel sheet containing the cps and the excel file only
    differs in the .xlsx we can define the excel sheet name with the name given as input:
    '''    
    #Load
    Dat = pd.read_excel(excel_name, cps_sheet_name, header = [1], index_col=0)
        #header 1 means take row 1 to give names to the columns
        #That contains the cps and cps*dil factor
        #index col = 0 to use first column as index!
    '''
    Note once I suffered that the dimesions of those were not similar, and in one sheet they
    were loadingn NaN values. I just erase those empty stuff in excel (selecting and delete)
    and after it worked!
    '''
    
    ############### 2) Clean df ############
    
    '''
    After the raw load, we can clean that a bit, creating a handful df, not the preovious, which
    are literally the excel in a df. That is, only collecting the relevant columns and putting them
    in a df, for further analysis (plotting, etc).
    
    Nevertheless, here could be relevant the fact of removing the blank or not, and for that you should
    say if you have a blank and if you want to delete it.
    
    The 1st cleaning is removing the 1st 4 rows, which contain bullshit, so we could do it with the 
    .drop method. Note the index are no longer used, simply erased.
    '''
    df_cps = Dat.drop(index = [Dat.index[0], Dat.index[1], Dat.index[2],Dat.index[3] ], axis = 0)   
                                                    #Removing rows 0, 1,2,3 (their index)
                
    #Another cleaning will be putting the df in numeric format. It is in object format, which gives problems
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



#%%######################################
########### 1.2) Future std computer!!!!!!!!!!!!!!!!!!!!!! #############
#####################################

def ICPMS_std_calculator (df_cps, df_rsd):
    '''
    Function that will compute the std from the %rsd data and the cps data. The data
    have already been read with the reader function.
    
    *Inputs:
        .df_cps: df containing the cps/ppb data. Isotopes are indexes, each row a measurement, 1st 1st repl, then
		2nd replicate, etc
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
        Not tested, but was copied literally from a function that worked (the reader was bigger), 
        so it should work, right? not even changed the names xD
    
    ##################
            
        '''
    
    
    ########### 1) Calc ###########
    '''
    Note the sig is RSD = relative standard deviation. I should compute sigma (std dev),
    which could be plotted in the temporal plots for ex as error bar. Thats easy:
            RSD = sigma / <x> * 100 ==> sigma = RSD * <x> /100
    
    Note RSD I have in 1 df (for each measurement (column) I have lot of elements (rows) ),
    and 1 df is for RSD; the other is for the mean values, so I need to do operations between 
    them. 
        
    Still, note the cps are multiplied by the dilution factor Df and applied some corrections
    (blanks, IS). i will forget about the 2nd things, more or less minor, and will consider the
    Df correction, which is essentially multiplying by it. So, I should multiply the RSD by the Df,
    and then apply that
    
    '''
    df_std= pd.DataFrame(df_cps * df_rsd / 100, 
                  columns=df_cps.columns, index=df_cps.index)   #std df    
    
    '''
    TO apply the Df corrections or not we set a variable to choose, true or false
    '''    
        
    ############## 2) Ouptuts ######################

    return df_std

    


#%%######################################
########### 1.3) ICPMS Dilution factor finder #############
#####################################

def ICPMS_Df_finder (excel_name, D_f_data, samp_prep_sheet_name = 'Sample_prep'):
    '''
    Function that will find the ICPMS dilution factor from the excel containing the ICPMS sample
    preparation info.
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
            Df_data = [1, 3] means names in row 1, Df in row 3 (in the excel)
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
    Since its object type, I will make it numeric, since everything will be easier with it (and strictly its true)
    '''
    D_f = D_f.apply(pd.to_numeric)                    #conversion to numeric data type
    
    ########### 4) Return #############
    '''
    Finally we return Df, which is a pandas serie
    '''
    
    return D_f            #return




#%%######################################
########### 1.4) ICPMS Dilution factor corrector #############
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
    
    
    ########### 2) Return #############
    return df_corrected             #return



#%%######################################
########### 1.5) ICPMS Sample blank substraction #############
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
    
    Only for 2 replicates, but could be generalized for N replicates easily.


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



#########################################################################
#%% ############### 1.6) Get isotope number from name ########################
########################################################################

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




#%%######################################
########### 1.7) ICPMS IS sens calculation #############
#####################################

def IS_sens_calculator_plotter(df_cps_ppb_dat, 
                               IS_meas = ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)'],
                               name_IS_sens_LR_plot = 'IS_sensLR_plot', 
                               name_IS_sens_MR_plot = 'IS_sensMR_plot'):
    '''
    Function part of the ICPMS Data processing!
    
    Function th<t will compute the IS sens (cps/ppb) for a df containing the cps and ppb. The format is like
    Stefaans excel, his first sheet. It IS needed the ppb table below the cps data. How below? Do not care, I
    find the data using locate functions ;)
    
    Note the plots appear in the plots pannel in Spyder, but are also saved ;)

    
    * Input
        .df_cps_ppb_dat: df containing the cps and ppb data. THe fashion is like Stefaans raw sheet. Isotopes are index.
            Take care of the names of the columns (like in IS conc table or so, the values you find), if they dont exist will
            give error!
        .IS_meas: array containing in a list the measured Internal Standards, containing its resolution, like 
            how they appear in the isotopes column. Default value:
                ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)']
            That mean those isotopes were measured. If Ho165(MR) also measured, just included it, and fine ;)
        .name_IS_sens_LR_plot: name for the plot of the IS sens for LR case. Similar for MR. Default values:
            'IS_sensLR_plot' and 'IS_sensMR_plot'. To that is added the .png to create the file.
            
        
    #Output
        .df with the IS sens data
        .Plot with the IS sens plot, 2 plots, 1 for LR and other for MR, in .png file
        
        
    ###### TO Do
        .Function to compute the ppb data, from Sum isobars? May be complex, possibly too time consuming ?
    '''
    
    ############ Data finder ################
    '''
    1st thing to find the ppb data. Should be somehow below the cps data, which has Stefaans format.
    To do the find in a loop way I define arrays which the elements to find. 
    Given as an input now, easier, KISS!

    Note they are associated, you take i element of both arrays, they go in pairs. I ahve seen a way, is to create in the loop
    the lines, and add them separately. The loop find the elements and perform the division cps/ppb
    '''

    df_IS_sens = pd.DataFrame()         #Empty df to store the values


    for i in range(0, len(IS_meas)):
        value_to_find = IS_meas[i]
        
        #####ppb value #########
        '''
        This will be made from the IS_meas variable. Since each element has a different letter, I can just read the 1st letter and
        say:
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
            
        matching_rows = df_cps_ppb_dat.loc[df_cps_ppb_dat.index == value_to_find]  #give the full row!
        matching_rows_ppb = df_cps_ppb_dat.loc[df_cps_ppb_dat.index == value_to_find_ppb]  #give the full row!

        #Now the operations, element wise, in array mode to avoid the index problems
        aux = matching_rows.iloc[:,:].values / matching_rows_ppb.iloc[:,:].values
            #eleement wise operation, like this you dont care about having diferent indexes, since you are
            #multiplying arrays. I erased the 1st value (Co name), so I needit to add it again
        
        #TO store temporarily those values I create an auxiliary df
        df_aux = pd.DataFrame(data = aux)

        #And I add that to the storing df
        df_IS_sens = df_IS_sens.append(df_aux, ignore_index= True)

    '''
    Now, to add the isotopes as index, we first add them as a column, and then we set
    the column to index:we need to insert the isotopes column, giving as values the isotopes names or so:
    '''
    #df_IS_sens['Isotopes'] = ['Co LR', 'In LR', 'Ho LR', 'Th LR', 'Co MR', 'In MR']
    df_IS_sens['Isotopes'] = IS_meas        #setting isotopes name from the loop variable, better
    
    df_IS_sens.set_index('Isotopes', inplace = True)
        
    
    '''
    As a final stylish, I can put the same column names as the df used to create the IS df:
    '''
    df_IS_sens.columns = df_cps_ppb_dat.columns
    
    
    #################################################
    ############# IS sens plotter ####################
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
            axL.plot(list(range(0, df_IS_sens.shape[1])), df_IS_sens.iloc[i,:], label = df_IS_sens.index[i]) 
            
        else:   #MR
            axM.plot(list(range(0, df_IS_sens.shape[1])), df_IS_sens.iloc[i,:], label = df_IS_sens.index[i])  
    
    #Final styling of the plot
    axL.legend()
    axM.legend()
    axL.grid(True)
    axM.grid(True)
    axL.set_xlabel('Sample number', size = 14)
    axM.set_xlabel('Sample number', size = 14)
    axL.set_ylabel("cps/ppb", size = 14)
    axM.set_ylabel("cps/ppb", size = 14)   
    axL.tick_params(axis='both', labelsize=14)              #size of axis
    axM.tick_params(axis='both', labelsize=14)              #size of axis
    pltL.savefig(name_IS_sens_LR_plot + '.png', format='png', bbox_inches='tight')    #note I call plt, not ax!
    pltM.savefig(name_IS_sens_MR_plot + '.png', format='png', bbox_inches='tight')
    
         
    ############## Return of values ##########################
    return df_IS_sens            #mass is an string!



#%%######################################
########### 1.8) ICPMS IS sens correction #############
#####################################

def IS_sens_correction(df_raw, df_IS_sens, 
                       IS_meas = ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)']):
    '''
    PART OF THE ICPMS Data processing!
    
    Function that will apply the IS sens correction to the cps data, from the raw ICPMS data.
    This is aprt of the ICPMS data analysis. You need to correct for IS sens, then apply ICPMS
    blanks, and then calibrate with 2+4 and 3+5 IS.
    
    'IS conc [ppb]' should be the box that indicates where the ppb chart start. That is, in A column in excel,
    that should be written, since that is used to find the ppb data!!!
    
    Note the plots appear in the plots pannel in Spyder

    
    * Input
        .df_raw: df containing the raw cps and ppb data. THe fashion is like Stefaans raw sheet
        .df_IS_sens df containing the IS sens data (cps/ppb). Ensure the size andindexing are
        appropiates!
        .IS_meas: array containing in a list the measured Internal Standards, containing its resolution, like 
            how they appear in the isotopes column. Default value:
                ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)']
        
    #Output
        .df with the IS corrected data, containing the cps and ppb data
        
        
    ###### TO Do #############
        .Switch cases if other corrections needede (avoid one IS because the measuring was wrong?)
    #########################
    '''
    
    
    ############### 0) Import neccessary functions#################
    '''
    Well, not needed to import anything, since it is also here the function to get the mass number A!
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
    
    Now I can try to do a loop. I could give where the IS sens start from the df inspection (240), and I know the order:
        .Co LR
        .In LR
        .Ho R
        . Th LR
        .CO MR
        . In MR
        
    Next step is to create a function to get the mass and apply one or other correction.

    We already have the function, so now we need to create a loop (after a function) that get the mass, and
    apply one or other correction
    '''   
    
    df_IS_correct = df_raw.copy()       #Thats the proper way to create it, copy so it is not asigned to it!
    
    	#The loop should go until the alst isotope, which I can find by finding IS conc ppb, the first thing for the ppb chart!
        #so, THIS data is mandatory that exist like that!!
    '''
    To generalize, I will do if statement for the different IS measured cases. By far only 2, the sequence done in the past, the
    4 in LR and 2 in MR; and now (8/23) 4 in MR also, the most detailed case.
    '''
    
    if IS_meas == ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)']:   
                                                                    ########case 1, the old scenario (exp bef 8/23)
     
        for i in range(np.where(df_raw.index == 'IS conc [ppb]')[0][0]): #loop through all isotopes
                #df_raw.loc[229,'Isotopes'] = df_raw.iloc[226,0]#relation iloc, loc
        
            mass, resol = Get_A_Resol(df_raw.index[i])
        #No we neeed to see the resolution, since for LR we have 4, and for MR only 2
            if resol == 'LR':           #low resolution
            #Now, which mass?
                if int(mass) in range(59,81): #mass from 59 to 80 included both
                    df_IS_correct.loc[df_raw.index[i],df_IS_correct.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[0],df_IS_sens.columns[:] ] * df_IS_sens.loc[df_IS_sens.index[0],df_IS_sens.columns[:]].mean() 
        
                elif int(mass) in range (81, 139):
                    df_IS_correct.loc[df_raw.index[i],df_IS_correct.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[1],df_IS_sens.columns[:] ] * df_IS_sens.loc[df_IS_sens.index[1],df_IS_sens.columns[:]].mean() 
            
                elif int(mass) in range (139, 210):
                    df_IS_correct.loc[df_raw.index[i],df_IS_correct.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[2],df_IS_sens.columns[:] ] * df_IS_sens.loc[df_IS_sens.index[2],df_IS_sens.columns[:]].mean() 
            
                else:   #if mass in range(210, 249):    rest of mass range
                    df_IS_correct.loc[df_raw.index[i],df_IS_correct.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[3],df_IS_sens.columns[:] ] * df_IS_sens.loc[df_IS_sens.index[3],df_IS_sens.columns[:]].mean() 

            else:       #Medium resolution,  MR
                if int(mass) in range(59,81): #mass from 59 to 80 included both
                    df_IS_correct.loc[df_raw.index[i],df_IS_correct.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[4],df_IS_sens.columns[:] ] * df_IS_sens.loc[df_IS_sens.index[4],df_IS_sens.columns[:]].mean() 
            
                else:       #rest of mass ranges, only 2 IS here xD
                    df_IS_correct.loc[df_raw.index[i],df_IS_correct.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[5],df_IS_sens.columns[:] ] * df_IS_sens.loc[df_IS_sens.index[5],df_IS_sens.columns[:]].mean() 
     
    elif  IS_meas == ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 
               'Co59(MR)', 'In115(MR)', 'Ho165(MR)', 'Th232(MR)']:   #######case 2, all IS measured in LR and MR! 

        for i in range(np.where(df_raw.index == 'IS conc [ppb]')[0][0]): #loop through all isotopes
                #df_raw.loc[229,'Isotopes'] = df_raw.iloc[226,0]#relation iloc, loc
        
            mass, resol = Get_A_Resol(df_raw.index[i])
        #No we neeed to see the resolution, since for LR we have 4, and for MR only 2
            if resol == 'LR':           #low resolution
            #Now, which mass?
                if int(mass) in range(59,81): #mass from 59 to 80 included both
                    df_IS_correct.loc[df_raw.index[i],df_IS_correct.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[0],df_IS_sens.columns[:] ] * df_IS_sens.loc[df_IS_sens.index[0],df_IS_sens.columns[:]].mean() 
        
                elif int(mass) in range (81, 139):
                    df_IS_correct.loc[df_raw.index[i],df_IS_correct.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[1],df_IS_sens.columns[:] ] * df_IS_sens.loc[df_IS_sens.index[1],df_IS_sens.columns[:]].mean() 
            
                elif int(mass) in range (139, 210):
                    df_IS_correct.loc[df_raw.index[i],df_IS_correct.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[2],df_IS_sens.columns[:] ] * df_IS_sens.loc[df_IS_sens.index[2],df_IS_sens.columns[:]].mean() 
            
                else:   #if mass in range(210, 249):    rest of mass range
                    df_IS_correct.loc[df_raw.index[i],df_IS_correct.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[3],df_IS_sens.columns[:] ] * df_IS_sens.loc[df_IS_sens.index[3],df_IS_sens.columns[:]].mean() 

            else:       #Medium resolution,  MR
                if int(mass) in range(59,81): #mass from 59 to 80 included both
                    df_IS_correct.loc[df_raw.index[i],df_IS_correct.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[4],df_IS_sens.columns[:] ] * df_IS_sens.loc[df_IS_sens.index[4],df_IS_sens.columns[:]].mean() 
            
                elif int(mass) in range (81, 139):
                    df_IS_correct.loc[df_raw.index[i],df_IS_correct.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[5],df_IS_sens.columns[:] ] * df_IS_sens.loc[df_IS_sens.index[5],df_IS_sens.columns[:]].mean() 
            
                elif int(mass) in range (139, 210):
                    df_IS_correct.loc[df_raw.index[i],df_IS_correct.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[6],df_IS_sens.columns[:] ] * df_IS_sens.loc[df_IS_sens.index[6],df_IS_sens.columns[:]].mean() 
            
                else:   #if mass in range(210, 249):    rest of mass range
                    df_IS_correct.loc[df_raw.index[i],df_IS_correct.columns[:]] = df_raw.loc[df_raw.index[i],df_raw.columns[:]] / df_IS_sens.loc[df_IS_sens.index[7],df_IS_sens.columns[:] ] * df_IS_sens.loc[df_IS_sens.index[7],df_IS_sens.columns[:]].mean() 


    ################## Return ###############
    '''
    After th eloop ends, we can return it
    '''
    return df_IS_correct



#%%######################################
########### 1.9) ICPMS Blank correction #############
#####################################

def ICPMS_ICPMSBlanks_corrector(df_IS_co, columns_blks):
    
    '''
    PART OF THE ICPMS Data processing!
    
    Function to be run after the IS sensitivity correction, to apply the next step in the ICPMS data
    analysis process, the ICPMS blanks correction. 
    
    Note that in the excel, int he part of the ppb table, you should delete the names, that stefaan write twice,
    for this to work!!!! 
    

    
    
    *Inputs:
        .columns_blks: np.array([]) indicating the number of the columns contaning blanks (std 0ppt and blank std).
            Ex: columns_blanks = np.array([1, 2]). Those columns number are from excel!
        .df_IS_co: df containing the cps (also ppb, not needed but there it is), output from the IS sens correction funciton.
            The formatting is like the excel. 
    *Outputs:
        .df containinng the cps data corrected for the ICPMS blanks. Note it also contains the ppd data, but now modified so they
            are random numbers. 
            
    ##### To DO ###########
        *Improve style and remove ppb data (would need modifications for the IS sens function)
    '''
    
    ###################################################################
    '''
    So, for the Blank correction, I need:
        1) To put ina  df all of the blanks. I could indicate column numbers, KISS
        2) Compute average from the cps of the blank
        3) Substract either the average, or the minimum cps value (from cps IS corrected), to
            ensure no negative values. Well, since the values are always positive, the mean will
            always be superior thatn the average value! So, we substract the minimum value!

    I would really need the IS corrected data without the sens also bro, but I could get it easily I think
    '''
    
    columns_blks = columns_blks - 1 #to adapt to the system in the df, isotopes are index, not columns!
    
    #1) is trivial:
    df_blanks = df_IS_co.iloc[:, columns_blks ]          #df with blanks
                    #no isotope label contained!!

    #2) could be done with df
    blk_means = df_blanks.mean(axis = 1)

    #3) The substraction
    '''
    THe first step is getting the minimun values. Those values are for the columns of the samples except the blanks!
    '''
    df_Is_co_no_blk = df_IS_co.drop(df_IS_co.columns[columns_blks ], axis = 1)
                        #df initial without the blanks!! no isotope also!

    #I could do the min to that, but since also contain the cps data there are some NaN, which make things not work. If I do
    #fill nan with the mean values, that could fix it. Lets try bro! 

    df_Is_co_no_blk.fillna(99999, inplace = True)           #filling NaN values with 99999 (easily recognizable)

    #Now the min values are_

    min_values = df_Is_co_no_blk.min(axis = 1) 
                #axis = 1 for columns!
    '''
    For this I would need a loop to check for each row which value to substract. I would say we always substract
    the mean, and Stefaan also think that. adn I probe that matematically, so it is like that xD

    So, now the problem is perform that operation here on python, since the df are different.

    I needed to delete some column nsames in teh ppb data!!!!!!!!!!
    
    The loop is in the same fashion as the one for IS sens correction :D
    '''

    df_IS_blk_co = df_IS_co.copy(  )        #initilization df for the loop

    for i in range(np.where(df_IS_co.index == 'IS conc [ppb]')[0][0]): #loop through all isotopes
                    #[0]s we get the desired value!
        if min_values[i] <= blk_means[i]:     #min low, so substracting min
            
            df_IS_blk_co.loc[df_IS_co.index[i],df_IS_blk_co.columns[:]]= df_IS_co.loc[df_IS_co.index[i],df_IS_co.columns[:]] - min_values[i]
                            #substraction of the min value        
            
        else:           #mean lower, so substracting the mean
            df_IS_blk_co.loc[df_IS_co.index[i],df_IS_blk_co.columns[:]]= df_IS_co.loc[df_IS_co.index[i],df_IS_co.columns[:]] - blk_means[i]

    
    #FInally we delete the blanks before returning the data: 
        
    df_IS_blk_co_final = df_IS_blk_co.drop(df_IS_blk_co.columns[columns_blks], axis = 1)
                #for some reason I can not do the drop and inplace = True, dont work


    ########### Return ###########
    return df_IS_blk_co_final           #return of the data



#%%######################################
########### 1.10) ICPMS data processing automatized #############
#####################################

def ICPMS_data_process(df_cps, ICPblk_columns, 
                       name_plot_LR_bef = 'IS_sensLR_plotBEF', name_plot_MR_bef = 'IS_sensMR_plotBEF',
                       name_plot_LR_aft = 'IS_sensLR_plot', name_plot_MR_aft = 'IS_sensMR_plot',
                       IS_meas = ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)'],
                       excel_name = 'df_IS_Blks_corr.xlsx'):
    '''
    SUITE of ICPMS Data Processing!
    
    Function that will apply all the steps for the automatization of the ICPMS data processing:
        1) Read cps amd ppb data 
        2) Compute sens = cps/ppb and plot it
        3) Apply the sensitivity correction
        4) plot the new sensitivtiy plots (should be straight lines)
        5) Substract ICPMS blanks
        6) Save that data (to calibrate after, by hand, for the moment..)
    
    To do 2), its needed that the ppb data table begins with IS conc [ppb]!
    
    *Inputs:
        .df_cps: df containing the cps data and also the ppb data, in the classic format. Must not contain the wash, will give 
            errors (divide by zero). Take care of the names (like for the cps table), they are crutial for the sens calc and 
                correction! Also about the format, not rows with 0s etc. Take a lot of care!!!
        .IS_meas: array containing in a list the IS and its resolution, like how they appear in the isotopes column. Default value:
                ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)', 'Co59(MR)', 'In115(MR)']
            That mean those isotopes were measured. If Ho165(MR) also measured, just included it, and fine ;)
        .ICPblk_columns: np array containing the columns numbers where the ICPMS blanks are (blank and std 0ppt). Numbers from
            the excel (A = 0, B = 1, etc)
        .name_plot_L(M)R_bef(aft): name of the plots of the IS sensitivity before and after doing the IS correction. 
        .excel_name: string that will be the name of the file containing the output. Default: 'df_IS_Blks_corr.xlsx'
        
    *Output:
        .df containing the data applying the ICPMS blanks and sensitivity corrections :)
    
    
    To Do:
            .Automatize more stuff? 
            .Optimize excel saving!
    '''
    
    ###### 1) IS sens calc ######
    df_IS_sens = IS_sens_calculator_plotter(df_cps, IS_meas, name_IS_sens_LR_plot = name_plot_LR_bef, 
    name_IS_sens_MR_plot = name_plot_MR_bef)         #calculation the IS correction
                #I define names of the plot, so the other ones have the default name
    
    
    ###### 2) IS sens correction and new sens calc ##########
    df_IS_corrected = IS_sens_correction(df_cps, df_IS_sens, IS_meas)       #applying the IS correction

    df_IS_sens_co = IS_sens_calculator_plotter(df_IS_corrected, IS_meas, name_IS_sens_LR_plot= name_plot_LR_aft,
                                               name_IS_sens_MR_plot= name_plot_MR_aft)    #getting and plotting new IS sens
    
    
    ##### 3)ICPMS Blk correction #########
    df_IS_Blks_co = ICPMS_ICPMSBlanks_corrector(df_IS_corrected, ICPblk_columns)    #correcting for ICPMS blanks


    ##### 4) Saving and Output #########
    '''
    Here I want to save the df after IS correction, and after the Blk correction. Both steps would be nice, for debugging!
    
    '''
    writer = pd.ExcelWriter(excel_name, engine = 'xlsxwriter')      #excel writer

    df_IS_Blks_co.to_excel(writer, sheet_name = 'Blk_correction', startrow = 5, freeze_panes = (6, 1))            
                #saving to excel. i freeze row 6 and column 1, so I can see all the data in a good way :)
                
    df_IS_corrected.to_excel(writer, sheet_name = 'IS_correction', startrow = 5, freeze_panes = (6, 1))        #saving to excel in another sheet
            #Note it does not have perfect format, to optimize it!!!
            #THe start trow make that Co59 is on row 7, as it should be!
    df_IS_sens_co.to_excel(writer, sheet_name = 'IS_correction', startrow = 5 + df_IS_corrected.shape[0] + 2)
                        #putting the new IS sensitivity below the IS corrected data!!
    
    writer.save()                                           #critical step, save the excel xD 


    return df_IS_Blks_co




#%%######################################
########### 1.10) ICPMS Isotope selector #############
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





#%%######################################
########### 1.11) Kd calculaor #############
#####################################

def ICPMS_KdQe_calc (df_data, df_VoM_disol, df_m_be, Nrepl = 2):
    '''
    Function that will compute the distribution constant Kd and the adsorption quantity
    q_e from the ppb data obtained with ICPMS. Note that data must be corrected
    for the dilutions factors. 
    
    Based on the blk corrector function.
    
    THe distribution constant Kd is:
        K_d = (C_0 - C_eq)/C_eq * V/m = Q_e / C_eq;
    being            
        C_0 = initial concentration [M]
        C_eq = concentration at equilibrium [M]
        m = mass of dry bentonite [g]
        V = volume of the solution [L]
        Q_e = absorbed quantity in equil [g soluto / g bent]
        
    In our case, that we measure C in ppb = ng soluto /g disol, we need to mutiply by g tot / m bent, to
    achieve the same units in Q_e!!

    Necessary that the data
    contain no Div0, ensure in the excel by using the iferror(operation, 0) function!
    Only for 2 replicates, but could be generalized for N replicates easily.
    
    Note this requires a df series with the volume, that you were not measuring in the first exp
    (up to 8/23). Note ten that units are involved!. If measuring mass ing and volumes in L, Q


    *Inputs:
        .df_data: dataframe containing the data, the full data, with the 2 replicates. Should be Dfs corrected
            Format: isotopes as index, columns the samples, 1st 1st replicate, then 2nd replicate. 2 replicates assume
            this function!!!!
        .df_VoM_disol: pd series containing the volume [mL] added to the bottle of the solution, BIC, 
        or whatever. normally 50ml OR the total mass of the solution [g]. If df_data in ppb, this must be the total mass
            so that Q_e is in g/g !
        .df_m_bent: pd series contaning the mass of bentonite [g] in the bottle (normally 250mg)
        .Nrepl: number of replicates. Default value = 2. 3 also accepted
    
    *Outputs (in that order):
        .df with the Kd data
        .df with q_e data
        
        '''
    
    
    ########## 0) Precalcs ##########
    '''
    To avoid the div0 error, which occurs when I have 0 as a values in the df, which I have for all the elements
    that were not found by ICPMS, I can just put NaN instead, since that will not give the Div0 error when computing Kd
    '''
    
    df_data.replace(0, np.nan, inplace=True)                    #replace 0 values with NaN, to avoid Div0 error!
    
    
    ########### 1) Calcs ###########
    '''
    The operations to perform are:
        1) C_0 - C_eq (>0)
        2) 1) * V/m = q_e
        3) 2)	1/C_eq = Kd
    
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
    
    #So, lets split into the 2 replicates!
    '''
    For 2 replicates its easy, for 3 it could be more tricky. Beware! TO create a function you should say
    the number of replicates and so!
    '''
    
    if Nrepl == 2:          #Standard case, 2 replicates
        df_1 = df_data.iloc[ :, 0: round( ( df_data.shape[1] ) / 2 ) ]      #1st replicate
        df_2 = df_data.iloc[ :, round( ( df_data.shape[1] ) / 2 ) :  ]       #replicate 2
    
        df_VoM_1 = df_VoM_disol.iloc[ 0: round( ( df_VoM_disol.shape[0] ) / 2 ) ]      #1st replicate
        df_VoM_2 = df_VoM_disol.iloc[ round( ( df_VoM_disol.shape[0] ) / 2 ) :  ]       #replicate 2
            #Achtung! In shape I put 0, because they are series, so 1D!!!!
            #VoM = Volume or Mass!
        df_m_1 = df_m_be.iloc[ 0: round( ( df_m_be.shape[0] ) / 2 ) ]      #1st replicat
        df_m_2 = df_m_be.iloc[ round( ( df_m_be.shape[0] ) / 2 ) :  ]       #replicate 2    
    
    #######Future note: here you see the automatization to N-replicates, doing this with a function.
        #Then the operations you can done them, grouping the df in an array, and for element in array, perform
        #them!
    
    ###### 1) C_0 - C_eq = - (C_eq - C0)
    #I will do C_eq - C0, and then invert that, since its easier. C0 is the blank data, 
    #thats why is easier, so I can copy paste the blank substraction
    
        dfCeq_C0_1 = df_1.subtract(df_1.iloc[:,0], axis = 0 )       #doing the substraction
        dfCeq_C0_1.drop( [df_1.iloc[:,0].name], axis = 1, inplace = True)   #drop blank column
        #
        dfCeq_C0_2 = df_2.subtract(df_2.iloc[:,0], axis = 0 )               #Replicate 2
        dfCeq_C0_2.drop( [df_2.iloc[:,0].name], axis = 1, inplace = True)
    
        #Now lets invert the sign:
        dfC0_Ceq_1 = - dfCeq_C0_1
        dfC0_Ceq_2 = - dfCeq_C0_2

    ######## 2) Apply the V/ m giving q_e (from Df_exp)
    #For this I ned to remove the blank columns to both m and V, since from C0-Ceq they are removed!

        df_m_1 = df_m_1[1:]         #fast way to delete 1st elemen (blank) in a series
                        #new_series = data.drop(data.index[0]) also work, from Chatgpt
        df_m_2 = df_m_2[1:]
        df_VoM_1 = df_VoM_1[1:]
        df_VoM_2 = df_VoM_2[1:]
    
    #And now I can operate:
        
        df_Qe_1 = dfC0_Ceq_1 * df_VoM_1 / df_m_1
        df_Qe_2 = dfC0_Ceq_2 * df_VoM_2 / df_m_2 
    
    ######## 3) Apply 1/C_eq = Kd
        df_Kd_1 = df_Qe_1 / df_1.drop( [df_1.iloc[:,0].name], axis = 1)   
                        #Not df_1 contains blk (1st column), so I remove it for the operation!    
                        #This also works: df_Qe_1.div(df_1.drop( [df_1.iloc[:,0].name], axis = 1) )
        df_Kd_2 = df_Qe_2 / df_2.drop( [df_2.iloc[:,0].name], axis = 1)
    
    #Now lets add them together

        df_Kd = pd.concat( [df_Kd_1, df_Kd_2], axis = 1)         #mergind the 2 little df ina  huge one
        df_Qe = pd.concat( [df_Qe_1, df_Qe_2 ] , axis = 1)

    elif Nrepl == 3:            #3 replicates
    #Gathering the replicates sepparately    
        df_1 = df_data.iloc[:, : round(df_data.shape[1] / 3)]
        df_2 = df_data.iloc[:, round(df_data.shape[1] / 3): 2*round(df_data.shape[1] / 3)]
        df_3 = df_data.iloc[:, 2*round(df_data.shape[1] / 3) :]
        
        df_VoM_1 = df_VoM_disol.iloc[ 0: round( ( df_VoM_disol.shape[0] ) / 3 ) ]      #1st replicate
        df_VoM_2 = df_VoM_disol.iloc[ round( ( df_VoM_disol.shape[0] ) / 3 ): 2* round( ( df_VoM_disol.shape[0] ) / 3 ) ]      
        df_VoM_3 = df_VoM_disol.iloc[ 2 *round( ( df_VoM_disol.shape[0] ) / 3 ) : ]      #3rd replicate
        
        df_m_1 = df_m_be.iloc[ : round( ( df_m_be.shape[0] ) / 3 ) ]      #1st replicat
        df_m_2 = df_m_be.iloc[ round( ( df_m_be.shape[0] ) / 3 ) : 2* round( ( df_m_be.shape[0] ) / 3 )]      #2nd replicat
        df_m_3 = df_m_be.iloc[ 2 *round( ( df_m_be.shape[0] ) / 3 ) : ]      #3rd replicat
    
        #1) 
        dfCeq_C0_1 = df_1.subtract(df_1.iloc[:,0], axis = 0 )       #doing the substraction
        dfCeq_C0_2 = df_2.subtract(df_2.iloc[:,0], axis = 0 )       #doing the substraction        
        dfCeq_C0_3 = df_3.subtract(df_3.iloc[:,0], axis = 0 )       #doing the substraction
        
        dfCeq_C0_1.drop( [df_1.iloc[:,0].name], axis = 1, inplace = True)   #drop blank column
        dfCeq_C0_2.drop( [df_2.iloc[:,0].name], axis = 1, inplace = True)   #drop blank column
        dfCeq_C0_3.drop( [df_3.iloc[:,0].name], axis = 1, inplace = True)   #drop blank column        
    
        dfC0_Ceq_1 = - dfCeq_C0_1
        dfC0_Ceq_2 = - dfCeq_C0_2
        dfC0_Ceq_3 = - dfCeq_C0_3
    
    ######## 2) Apply the V/ m giving q_e (from Df_exp)
    #For this I ned to remove the blank columns to both m and V, since from C0-Ceq they are removed!

        df_m_1 = df_m_1[1:]         #fast way to delete 1st elemen (blank) in a series
        df_m_2 = df_m_2[1:]
        df_m_3 = df_m_3[1:]
        df_VoM_1 = df_VoM_1[1:]
        df_VoM_2 = df_VoM_2[1:]
        df_VoM_3 = df_VoM_3[1:]
        
    #And now I can operate:
        
        df_Qe_1 = dfC0_Ceq_1 * df_VoM_1 / df_m_1
        df_Qe_2 = dfC0_Ceq_2 * df_VoM_2 / df_m_2 
        df_Qe_3 = dfC0_Ceq_3 * df_VoM_3 / df_m_3
        
        
    ######## 3) Apply 1/C_eq = Kd
        df_Kd_1 = df_Qe_1 / df_1.drop( [df_1.iloc[:,0].name], axis = 1)   
        df_Kd_2 = df_Qe_2 / df_2.drop( [df_2.iloc[:,0].name], axis = 1)
        df_Kd_3 = df_Qe_3 / df_3.drop( [df_3.iloc[:,0].name], axis = 1)
    
    #Now lets add them together
        df_Kd = pd.concat( [df_Kd_1, df_Kd_2, df_Kd_3], axis = 1)         
        df_Qe = pd.concat( [df_Qe_1, df_Qe_2, df_Qe_3 ] , axis = 1)  
    
    else:               #Error case
        print('Erro case, wrong Nrepl introduced!')    
        df_Kd = 0
        df_Qe =0
                
        
    ########### 2) Return #############
    return df_Kd, df_Qe             #return



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
####################### PLOTTERS ####################################
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#%%######################################
########### 1.12) ICPMS Bar plotter #############
#####################################

def ICPMS_Barplotter (df_cps, df_rstd, folder_name = 'Bar_plots',
                      pre_title_plt = "Concentration of ", pre_save_name = 'Conc',
                      Elem_rel = ['Si28', 'Si29', 'Si30',
                              'Al27',
                              'Mg24', 'Mg25', 'Mg26',
                              'Mn55',
                              'Fe56', 'Fe57',
                              'Ca42', 'Ca43', 'Ca44', 
                              'Na23', 
                              'K', 
                              'Ti46', 'Ti47', 'Ti48', 'Ti49', 'Ti50',
                              'Sr86', 'Sr87', 'Sr88', 'Cs133', 'U238', 'U235', 'U234', 'Eu151', 'Eu153', 'La139'] ):
    '''
    Function that will do bar plots of the raw data from the ICPMS, the cps and the rstd. By raw
    I mean withoutany corrections/calibrations (neither the dilution factor correction). This is a
    preliminary plot, to check if everything is allright, with the rstd (rstd should be < 20%).
    
    *Inputs:
        .df_cps, df_rstd: dataframes containing the cps and the relative standard deviation. Those are
        outputs for the Read_ICPMS_excel function. Note the Isotopes column are the inde!
        .folder_name: folder name. Default value: 'Bar_plots'
        .pre_title_plt : title of the graph, part that appears before the name of the elements (thats why pre title).
            Detault value: "Concentration of " (note the space after of, so the element is not together with that!)
        . pre_save_name: name of the graph files to save. Default: 'Conc', giving Conc_Mg24.png for ex    
        
        .Elem_rel: array containing the name of the relevant elemtns, which are the elements that will be saved
            in a specific folder. Default value: 
                ['Si28', 'Si29', 'Si30',
                        'Al27',
                        'Mg24', 'Mg25', 'Mg26',
                        'Mn55',
                        'Fe56', 'Fe57',
                        'Ca42', 'Ca43', 'Ca44', 
                        'Na23', 
                        'K', 
                        'Ti46', 'Ti47', 'Ti48', 'Ti49', 'Ti50',
                        'Sr86', 'Sr87', 'Sr88', 'Cs133', 'U238', 'U235', 'U234', 'Eu151', 'Eu153', 'La139']      
                #List of relevant elements. Note T sheet made of Si and Al,
                                            #Oct sheet by Al, Mg, Mn, Fe, and rest in interlaminar spaces.
                                            #commoninterlaminars are Ca, Na, K, Ti, li, Sr.
                                    #Previous elements that were deleted: S32,33,34, P31 (impurities)
        
    *Outputs:
        .Plots (saving them) of the raw data
    
    ######### To DO #########
    .TInclude option to plot all the elements or only the relevants, since for all (250) take 140s!
    
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
    This is a loop plot, so beware, will take long (2-3mins!).
    
    
    ##### Mutiple bar plot; number of bars = 2 ############

The relations that must be satisfied are, if the width of the bar is w and the 
blank space between values (value is x value) is b<1, is:

        centroid = value - w/2
        centroid (right of value) = value + w/2
        2w + b = 1 ==> w = (1-b)/2
        
Setting b gives w. In fact the general equations for 2n bars per X tick (n = 1,2,..) are:
    w = (1-b)/2n
    value - w/2*n <= centroid <= value + w/2*n, in steps of w (of course)
    
    ##############################
    '''
    t_start = tr.time()       #[s] start time of the plot execution
    
    ###Plot
    #Some parameters for the plot
    X_axis = np.arange( len(df_cps.axes[1]))                 #To do the 2 bar plot
            #choosing the number of columns. [0] for rows
    b = .4                              #[au] blank space between values <1
    w = (1-b)/2          #bar width

    for i in list( range(df_cps.shape[0] ) ):     #Loop thorugh all rows
                    #df_cps.index give the index values, low and high
                    #
                    ########### Bar plot ###########################
        plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
        plt.title(pre_title_plt + df_cps.index[i][:-4], fontsize=22, wrap=True)           #title
        a = plt.bar(X_axis - w/2, df_cps.loc[df_cps.index[i]], width = w, edgecolor="black", 
            label = "I ", align='center') 
            #-2 not to plot the blank!! Remove it to plot it!
        plt.ylabel("I [cps]", fontsize=14)              #ylabel
        plt.xlabel('Sample', fontsize = 14)
        plt.tick_params(axis='both', labelsize=14)              #size of axis
        plt.yscale('log') 
        plt.grid(True)
        plt.xticks(X_axis, [ df_cps.columns[:][j] for j in range(len(df_cps.axes[1]) ) ]
            , rotation = 90)
        plt.twinx()             #For setting 2 axes
        aa = plt.bar( X_axis + w/2, df_rstd.loc[df_cps.index[i]], width = w, edgecolor="black", 
             label = '$\sigma_{rel}$', align='center', color = 'red') 
        plt.ylabel("$\sigma_{rel}$ [%]", fontsize=14)              #ylabel
        #
        aaa = [a, aa]
        plt.legend(aaa, [p_.get_label() for p_ in aaa])
        
        #Saving in the folder
        if df_cps.index[i][:-4] in Elem_rel:  #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
            plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        'Conc_rsd_' + df_cps.index[i][:-4] + '.png', format='png', bbox_inches='tight')
            #
        else:        #if the element is not relevant
            plt.savefig(folder_name +'/' +  
                        'Conc_rsd_' + df_cps.index[i][:-4] +'.png', format='png', bbox_inches='tight')
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
    
    

#%%######################################
########### 1.12) ICPMS plotter #############
#####################################

def ICPMS_Plotter (x, df_cps, x_label, y_label, folder_name = 'Plots', 
                   pre_title_plt = "Concentration of ", pre_save_name = 'Conc',
                   Elem_rel =['Si28', 'Si29', 'Si30',
                           'Al27',
                           'Mg24', 'Mg25', 'Mg26',
                           'Mn55',
                           'Fe56', 'Fe57',
                           'Ca42', 'Ca43', 'Ca44', 
                           'Na23', 
                           'K', 
                           'Ti46', 'Ti47', 'Ti48', 'Ti49', 'Ti50',
                           'Sr86', 'Sr87', 'Sr88', 'Cs133', 'U238', 'U235', 'U234', 'Eu151', 'Eu153', 'La139'], 
                   plot_everything = False, Nrepl = 2 ):
    '''
    Function that will plots of the data from the ICPMS (cps) vs another variable, initially
    time, the cps and the rstd. This assume we have 2 replicates, 1 series after the other.
    Stimated running time around 80s.
    
    *Inputs:
        .x: x axis variable in the plot. This should be a df series
        .df_cps: dataframes containing the cps. Those are
        outputs for the Read_ICPMS_excel function. Note the 1st column must be the one
        with the isotopes (like in the excel)
        .x_label: string that will be the x label for the plot (for math stuff, 
                                    use $$. eg: '$\Delta t[h]$')
        .y_label: string that will be the y label for the plot
        .folder_name: string defining the name of the folder to create to store the plots
            default value: 'Plots'
        . plot_everything: string defining if you want to plot all the elements or only the
            relevant ones. Default value: False (only plot relevants)
        .pre_title_plt : title of the graph, part that appears before the name of the elements (thats why pre title).
            Detault value: "Concentration of " (note the space after of, so the element is not together with that!)
        . pre_save_name: name of the graph files to save. Default: 'Conc', giving Conc_Mg24.png for ex    
        . Nrepl: number of replicates. Default value: 2. 3 also accepted
        .Elem_rel: array containing the name of the relevant elemtns, which are the elements that will be saved
            in a specific folder. Default value: 
                ['Si28', 'Si29', 'Si30',
                        'Al27',
                        'Mg24', 'Mg25', 'Mg26',
                        'Mn55',
                        'Fe56', 'Fe57',
                        'Ca42', 'Ca43', 'Ca44', 
                        'Na23', 
                        'K', 
                        'Ti46', 'Ti47', 'Ti48', 'Ti49', 'Ti50',
                        'Sr86', 'Sr87', 'Sr88', 'Cs133', 'U238', 'U235', 'U234', 'Eu151', 'Eu153', 'La139']      
                #List of relevant elements. Note T sheet made of Si and Al,
                                            #Oct sheet by Al, Mg, Mn, Fe, and rest in interlaminar spaces.
                                            #commoninterlaminars are Ca, Na, K, Ti, li, Sr.
                                    #Previous elements that were deleted: S32,33,34, P31 (impurities)
                                    
    *Outputs:
        .Plots (saving them) of the x and df_cps data, cps vs x!
    
    
    ### TO DO: ####
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	.Implement error plotting (in an errorbar pyplot)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
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
    This is a loop plot, so beware, will take long if you plot all the elements (280) (2-3mins!).
    THe best way is loop with python style
    '''
    t_start = tr.time()       #[s] start time of the plot execution
    
    '''
    Before plotting, I need to see whether I have 2 or 3 replicates, to divide the df in one or another way.
    That will be with an if loop, as usual
    '''
    
    if Nrepl == 2:          #2 replicates, standard case
    
    ###Plot
        for index, row in df_cps.iterrows():   #Loop throught all the rows, the isotopes
                    #df_cps.index give the index values, low and high
		   # 4 because of the way the df is created (and hence the excel tabelle)

            #Saving in the folder
            if index[:-4] in Elem_rel:  #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
                plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                plt.title(pre_title_plt + index[:-4], fontsize=22, wrap=True)           #title
                plt.plot(x[:int(len(x)/2)], row[:int(len(x)/2)], 'bo--', 
                     MarkerSize = 5, label = 'Repl_1') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x[int(len(x)/2):], row[int(len(x)/2) :], 'ro--', 
                     MarkerSize = 5, label = 'Repl_2') 
                plt.ylabel(y_label, fontsize=14)              #ylabel
                plt.xlabel(x_label, fontsize = 14)
                plt.tick_params(axis='both', labelsize=14)              #size of axis
                #plt.yscale('log') 
                plt.grid(True)
                plt.legend()
                plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_' + index[:-4] + '.png', format='png', bbox_inches='tight')
            #
            else:        #if the element is not relevant
                if plot_everything == True :        #if you want to plot all the elements (may be desired?)
                #    
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + index[:-4], fontsize=22, wrap=True)     #title
                    plt.plot(x[:int(len(x)/2)], row[:int(len(x)/2)], 'bo--', 
                         MarkerSize = 5, label = 'Repl_1') 
                    plt.plot(x[int(len(x)/2):], row[int(len(x)/2):], 'ro--', 
                         MarkerSize = 5, label = 'Repl_2') 
                    plt.ylabel(y_label, fontsize=14)              #ylabel
                    plt.xlabel(x_label, fontsize = 14)
                    plt.tick_params(axis='both', labelsize=14)              #size of axis
                    #plt.yscale('log') 
                    plt.grid(True)
                    plt.legend()            
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
            if index[:-4] in Elem_rel:  #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
                plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                plt.title(pre_title_plt + index[:-4], fontsize=22, wrap=True)           #title
                plt.plot(x[:int(len(x)/3)], row[:int(len(x)/3)], 'bo--', 
                     MarkerSize = 5, label = 'Repl_1') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x[int(len(x)/3): 2* int(len(x)/3)], row[int(len(x)/3) :2* int(len(x)/3)], 'ro--', 
                     MarkerSize = 5, label = 'Repl_2') 
                plt.plot(x[2* int(len(x)/3):], row[2* int(len(x)/3) :], 'go--', 
                     MarkerSize = 5, label = 'Repl_3') 
                plt.ylabel(y_label, fontsize=14)              #ylabel
                plt.xlabel(x_label, fontsize = 14)
                plt.tick_params(axis='both', labelsize=14)              #size of axis
                #plt.yscale('log') 
                plt.grid(True)
                plt.legend()
                plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_' + index[:-4] + '.png', format='png', bbox_inches='tight')
            #
            else:        #if the element is not relevant
                if plot_everything == True :        #if you want to plot all the elements (may be desired?)
                #    
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + index[:-4], fontsize=22, wrap=True)     #title
                    plt.plot(x[:int(len(x)/3)], row[:int(len(x)/3)], 'bo--', 
                     MarkerSize = 5, label = 'Repl_1') 
                    #+1 needed since the df contain a row with the column names!
                    plt.plot(x[int(len(x)/3): 2* int(len(x)/3)], row[int(len(x)/3) :2* int(len(x)/3)], 'ro--', 
                     MarkerSize = 5, label = 'Repl_2') 
                    plt.plot(x[2* int(len(x)/3):], row[2* int(len(x)/3) :], 'go--', 
                     MarkerSize = 5, label = 'Repl_3') 
                    plt.ylabel(y_label, fontsize=14)              #ylabel
                    plt.xlabel(x_label, fontsize = 14)
                    plt.tick_params(axis='both', labelsize=14)              #size of axis
                    #plt.yscale('log') 
                    plt.grid(True)
                    plt.legend()            
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
    
    
    
 

#%%######################################
########### 1.13) ICPMS plotter 3 bentonites #############
#####################################

def ICPMS_Plotter3 (x, df_cps, x_label, y_label, folder_name = 'Plots', plot_everything = False,
                    pre_title_plt = "Concentration of ", pre_save_name = 'Conc',
                    Elem_rel = ['Si28', 'Si29', 'Si30',
                            'Al27',
                            'Mg24', 'Mg25', 'Mg26',
                            'Mn55',
                            'Fe56', 'Fe57',
                            'Ca42', 'Ca43', 'Ca44', 
                            'Na23', 
                            'K', 
                            'Ti46', 'Ti47', 'Ti48', 'Ti49', 'Ti50',
                            'Sr86', 'Sr87', 'Sr88', 'Cs133', 'U238', 'U235', 'U234', 'Eu151', 'Eu153', 'La139'] ):
    '''
    Function that will plots of the data from the ICPMS (cps) vs another variable, initially
    time, for the 3 bentonites. This assume we have 2 replicates, 1 series after the other.
    Stimated running time around 80s.
    
    *Inputs:
        .x: x axis variable in the plot. dictionary, contaning the 3 df series with the keys: Sard, Tur, BK
        .df_cps: dataframes containing the cps. dictionary containing the 3 df series, same order as x. Those are
        outputs for the Read_ICPMS_excel function. Note the 1st column must be the one
        with the isotopes (like in the excel)
        .x_label: string that will be the x label for the plot (for math stuff, 
                                    use $$. eg: '$\Delta t[h]$')
        .y_label: string that will be the y label for the plot
        .folder_name: string defining the name of the folder to create to store the plots
            default value: 'Plots'
        . plot_everything: string defining if you want to plot all the elements or only the
            relevant ones. Default value: False (only plot relevants)
        .pre_title_plt : title of the graph, part that appears before the name of the elements (thats why pre title).
            Detault value: "Concentration of " (note the space after of, so the element is not together with that!)
        . pre_save_name: name of the graph files to save. Default: 'Conc', giving Conc_Mg24.png for ex    
        .Elem_rel: array containing the name of the relevant elemtns, which are the elements that will be saved
            in a specific folder. Default value: 
                ['Si28', 'Si29', 'Si30',
                        'Al27',
                        'Mg24', 'Mg25', 'Mg26',
                        'Mn55',
                        'Fe56', 'Fe57',
                        'Ca42', 'Ca43', 'Ca44', 
                        'Na23', 
                        'K', 
                        'Ti46', 'Ti47', 'Ti48', 'Ti49', 'Ti50',
                        'Sr86', 'Sr87', 'Sr88', 'Cs133', 'U238', 'U235', 'U234', 'Eu151', 'Eu153', 'La139']      
                #List of relevant elements. Note T sheet made of Si and Al,
                                            #Oct sheet by Al, Mg, Mn, Fe, and rest in interlaminar spaces.
                                            #commoninterlaminars are Ca, Na, K, Ti, li, Sr.
                                    #Previous elements that were deleted: S32,33,34, P31 (impurities)
                                    
    *Outputs:
        .Plots (saving them) of the x and df_cps data, cps vs x!
    
    
    ### TO DO: ####
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	.Implement error plotting (in an errorbar pyplot)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    .Improve it, not soo god, averaging replicates and plotting average and its error would be better!
    . TO include 3 replicates version!
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
    This is a loop plot, so beware, will take long if you plot all the elements (280) (2-3mins!).
    
    '''
    Bent_color = {'Sard' : (.68,.24,.31), 'Tur' :  '#EEE8AA', 'BK' : 'grey'} 
                #Color for the ploting, color amtch the actual bentonite color
                
    t_start = tr.time()       #[s] start time of the plot execution
    
    ###Plot

    for i in list( range(df_cps['Sard'].shape[0] ) ):       #Loop thorugh all rows
		   # 4 because of the way the df is created (and hence the excel tabelle)
        #
        
        #Saving in the folder
        if df_cps['Sard'].index[i][:-4] in Elem_rel:        #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
            plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
            plt.title(pre_title_plt + df_cps['Sard'].index[i][:-4], fontsize=22, wrap=True)           #title
            #PLot bentonite 1, Sard
            plt.plot(x['Sard'][:int(len(x['Sard'])/2)], df_cps['Sard'].loc[df_cps['Sard'].index[i] ][:int(len(x['Sard'])/2)], 'o--', color = Bent_color['Sard'],
                     MarkerSize = 5, label = 'Repl_1 S') 
                    #+1 needed since the df contain a row with the column names!
            plt.plot(x['Sard'][int(len(x['Sard'])/2):], df_cps['Sard'].loc[df_cps['Sard'].index[i] ][int(len(x['Sard'])/2):], 'o--', color = Bent_color['Sard'],
                     MarkerSize = 5, label = 'Repl_2 S') 
            #PLot bentonite 2, T
            plt.plot(x['Tur'][:int(len(x['Tur'])/2)], df_cps['Tur'].loc[df_cps['Sard'].index[i] ][:int(len(x['Tur'])/2 )], 'o--', color = Bent_color['Tur'],
                     MarkerSize = 5, label = 'Repl_1 T') 
                    #+1 needed since the df contain a row with the column names!
            plt.plot(x['Tur'][int(len(x['Tur'])/2):], df_cps['Tur'].loc[df_cps['Sard'].index[i] ][int(len(x['Tur'])/2 ):], 'o--', color = Bent_color['Tur'],
                     MarkerSize = 5, label = 'Repl_2 T') 
            #PLot bentonite 3, BK
            plt.plot(x['BK'][:int(len(x['BK'])/2)], df_cps['BK'].loc[df_cps['Sard'].index[i] ][:int(len(x['BK'])/2 )], 'o--', color = Bent_color['BK'],
                     MarkerSize = 5, label = 'Repl_1 BK') 
                    #+1 needed since the df contain a row with the column names!
            plt.plot(x['BK'][int(len(x['BK'])/2):], df_cps['BK'].loc[df_cps['Sard'].index[i] ][int(len(x['BK'])/2 ):], 'o--', color = Bent_color['BK'],
                     MarkerSize = 5, label = 'Repl_2 BK') 
            plt.ylabel(y_label, fontsize=14)              #ylabel
            plt.xlabel(x_label, fontsize = 14)
            plt.tick_params(axis='both', labelsize=14)              #size of axis
            #plt.yscale('log') 
            plt.grid(True)
            plt.legend()
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
                     MarkerSize = 5, label = 'Repl_1 S') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x['Sard'][int(len(x['Sard'])/2):], df_cps['Sard'].loc[df_cps['Sard'].index[i] ][int(len(x['Sard'])/2):], 'o--', color = Bent_color['Sard'],
                     MarkerSize = 5, label = 'Repl_2 S') 
                    #PLot bentonite 2, T
                plt.plot(x['Tur'][:int(len(x['Tur'])/2)], df_cps['Tur'].loc[df_cps['Sard'].index[i] ][:int(len(x['Tur'])/2 )], 'o--', color = Bent_color['Tur'],
                     MarkerSize = 5, label = 'Repl_1 T') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x['Tur'][int(len(x['Tur'])/2):], df_cps['Tur'].loc[df_cps['Sard'].index[i] ][int(len(x['Tur'])/2 ):], 'o--', color = Bent_color['Tur'],
                     MarkerSize = 5, label = 'Repl_2 T') 
                    #PLot bentonite 3, BK
                plt.plot(x['BK'][:int(len(x['BK'])/2)], df_cps['BK'].loc[df_cps['Sard'].index[i] ][:int(len(x['BK'])/2 )], 'o--', color = Bent_color['BK'],
                     MarkerSize = 5, label = 'Repl_1 BK') 
                    #+1 needed since the df contain a row with the column names!
                plt.plot(x['BK'][int(len(x['BK'])/2):], df_cps['BK'].loc[df_cps['Sard'].index[i] ][int(len(x['BK'])/2 ):], 'o--', color = Bent_color['BK'],
                     MarkerSize = 5, label = 'Repl_2 BK') 
                plt.ylabel(y_label, fontsize=14)              #ylabel
                plt.xlabel(x_label, fontsize = 14)
                plt.tick_params(axis='both', labelsize=14)              #size of axis
                #plt.yscale('log') 
                plt.grid(True)
                plt.legend()
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

    
    
#%%######################################
########### 1.14) ICPMS plotter blank appart #############
#####################################

def ICPMS_Plotter_blk (x, df_cps, x_label, y_label, folder_name = 'Plots', plot_everything = False, 
                       pre_title_plt = "Concentration of ", pre_save_name = 'Conc', Nrepl = 2,
                       Elem_rel = ['Si28', 'Si29', 'Si30',
                               'Al27',
                               'Mg24', 'Mg25', 'Mg26',
                               'Mn55',
                               'Fe56', 'Fe57',
                               'Ca42', 'Ca43', 'Ca44', 
                               'Na23', 
                               'K', 
                               'Ti46', 'Ti47', 'Ti48', 'Ti49', 'Ti50',
                               'Sr86', 'Sr87', 'Sr88', 'Cs133', 'U238', 'U235', 'U234', 'Eu151', 'Eu153', 'La139'] ):
    '''
    Function that will plots of the data from the ICPMS (cps) vs another variable, initially
    time, the cps and the rstd. This assume we have 2 replicates, 1 series after the other.
    Blank is plotted sepparately, so the data must include a blank, which should be 1st, the number 1
    
    *Inputs:
        .x: x axis variable in the plot. This should be a df series
        .df_cps: dataframes containing the cps. Those are
        outputs for the Read_ICPMS_excel function. Note the isotopes are the index, so 1st column is 1_1!
        .x_label: string that will be the x label for the plot (for math stuff, 
                                    use $$. eg: '$\Delta t[h]$')
        .y_label: string that will be the y label for the plot
        .folder_name: string defining the name of the folder to create to store the plots
            default value: 'Plots'
        . plot_everything: string defining if you want to plot all the elements or only the
            relevant ones. Default value: False (only plot relevants)
        .pre_title_plt : title of the graph, part that appears before the name of the elements (thats why pre title).
                Detault value: "Concentration of " (note the space after of, so the element is not together with that!)
        . pre_save_name: name of the graph files to save. Default: 'Conc', giving Conc_Mg24.png for ex    
        .Nrepl : number of replicates. Default value : 2. 3 value also accepted
        .Elem_rel: array containing the name of the relevant elemtns, which are the elements that will be saved
            in a specific folder. Default value: 
                ['Si28', 'Si29', 'Si30',
                        'Al27',
                        'Mg24', 'Mg25', 'Mg26',
                        'Mn55',
                        'Fe56', 'Fe57',
                        'Ca42', 'Ca43', 'Ca44', 
                        'Na23', 
                        'K', 
                        'Ti46', 'Ti47', 'Ti48', 'Ti49', 'Ti50',
                        'Sr86', 'Sr87', 'Sr88', 'Cs133', 'U238', 'U235', 'U234', 'Eu151', 'Eu153', 'La139']      
                #List of relevant elements. Note T sheet made of Si and Al,
                                            #Oct sheet by Al, Mg, Mn, Fe, and rest in interlaminar spaces.
                                            #commoninterlaminars are Ca, Na, K, Ti, li, Sr.
                                    #Previous elements that were deleted: S32,33,34, P31 (impurities)
                                    
    *Outputs:
        .Plots (saving them) of the x and df_cps data, cps vs x!
    
    
    ### TO DO: ####
	.Implement error plotting (in an errorbar pyplot)
    .Implement variable for modifying tittle (eg adding bentonite name, as a input)
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
    This is a loop plot, so beware, will take long if you plot all the elements (280) (2-3mins!).
    I inlcude in if statement the numer of replicates, currently only 2 and 3!

    '''
    t_start = tr.time()       #[s] start time of the plot execution
    
    if Nrepl == 2:               #2 replicates, standard case!
    ###Plot
        for index,row in df_cps.iterrows(): 
                    #df_cps.index give the index values, low and high
		   # 4 because of the way the df is created (and hence the excel tabelle)

            #Saving in the folder
            if index[:-4] in Elem_rel:                      #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
                plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                plt.title(pre_title_plt + index[:-4], fontsize=22, wrap=True)     #title
                plt.plot(x[:int(len(x)/2) ], row[:int(len(x)/2)], 'bo--', 
                     MarkerSize = 5, label = 'Repl_1') 
                plt.hlines(row[0], min(x[:int(len(x)/2)]), max(x[:int(len(x)/2)] ), label = 'Blk_1' )
                    #row[0] is 1_1, 1st sample of 1st replicate!
            #Now replicate 2
                plt.plot(x[int(len(x)/2):], row[int(len(x)/2 ):], 'ro--', 
                     MarkerSize = 5, label = 'Repl_2') 
            
                plt.hlines(row[int(len(x)/2)], min(x[int(len(x)/2):] ), max(x[int(len(x)/2):] ), label = 'Blk_2', colors = 'r' )
                plt.ylabel(y_label, fontsize=14)              #ylabel
                plt.xlabel(x_label, fontsize = 14)
                plt.tick_params(axis='both', labelsize=14)              #size of axis
            #plt.yscale('log') 
                plt.grid(True)
                plt.legend()
                plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_'  + index[:-4] + '.png', format='png', bbox_inches='tight')
            #
            else:        #if the element is not relevant
                if plot_everything == True :     #if you want to plot all the elements (may be desired?)
                #    
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + index[:-4], fontsize=22, wrap=True)     #title
                    plt.plot(x[:int(len(x)/2)], row[:int(len(x)/2)], 'bo--', 
                         MarkerSize = 5, label = 'Repl_1') 
                    plt.plot(x[int(len(x)/2):], row[int(len(x)/2):], 'ro--', 
                         MarkerSize = 5, label = 'Repl_2') 
                    plt.ylabel(y_label, fontsize=14)                #ylabel
                    plt.xlabel(x_label, fontsize = 14)
                    plt.tick_params(axis='both', labelsize=14)              #size of axis
                #plt.yscale('log') 
                    plt.grid(True)
                    plt.legend()            
                    plt.savefig(folder_name +'/' +  
                        pre_save_name + '_'  + index[:-4] +'.png', format='png', bbox_inches='tight')
                    #To save plot in folder
        
            plt.close()             #to clsoe the plot not to consume too much resources
            
    elif Nrepl == 3:                        #3 replicates case
        for index,row in df_cps.iterrows(): 
                    #df_cps.index give the index values, low and high
		   # 4 because of the way the df is created (and hence the excel tabelle)

            #Saving in the folder
            if index[:-4] in Elem_rel:                      #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
                plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                plt.title(pre_title_plt + index[:-4], fontsize=22, wrap=True)     #title
                plt.plot(x[:int(len(x)/3) ], row[:int(len(x)/3)], 'bo--', 
                     MarkerSize = 5, label = 'Repl_1') 
                plt.plot(x[int(len(x)/3): 2*int(len(x)/3) ], row[ int(len(x)/3): 2*int(len(x)/3) ], 'ro--', 
                     MarkerSize = 5, label = 'Repl_2') 
                plt.plot(x[ 2*int(len(x)/3): ], row[ 2*int(len(x)/3): ], 'go--', 
                     MarkerSize = 5, label = 'Repl_3') 
                plt.hlines(row[0], min(x[:int(len(x)/3)]), max(x[:int(len(x)/3)] ), label = 'Blk_1', colors = 'b' )
                    #row[0] is 1_1, 1st sample of 1st replicate!
                plt.hlines(row[int(len(x)/3)], min(x[int(len(x)/3): 2*int(len(x)/3)] ), max(x[int(len(x)/3): 2* int(len(x)/3)] ), 
                           label = 'Blk_2', colors = 'r' )
                plt.hlines(row[2* int(len(x)/3)], min(x[2*int(len(x)/3):] ), max(x[2* int(len(x)/3):] ), 
                           label = 'Blk_3', colors = 'g' )
                plt.ylabel(y_label, fontsize=14)              #ylabel
                plt.xlabel(x_label, fontsize = 14)
                plt.tick_params(axis='both', labelsize=14)              #size of axis
            #plt.yscale('log') 
                plt.grid(True)
                plt.legend()
                plt.savefig(folder_name + '/' + 'Relevants' + '/' +
                        pre_save_name + '_'  + index[:-4] + '.png', format='png', bbox_inches='tight')
            #
            else:        #if the element is not relevant
                if plot_everything == True :     #if you want to plot all the elements (may be desired?)
                #    
                    plt.figure(figsize=(11,8))          #width, heigh 6.4*4.8 inches by default
                    plt.title(pre_title_plt + index[:-4], fontsize=22, wrap=True)     #title
                    plt.plot(x[:int(len(x)/2)], row[:int(len(x)/2)], 'bo--', 
                         MarkerSize = 5, label = 'Repl_1') 
                    plt.plot(x[int(len(x)/2):], row[int(len(x)/2):], 'ro--', 
                         MarkerSize = 5, label = 'Repl_2') 
                    plt.ylabel(y_label, fontsize=14)                #ylabel
                    plt.xlabel(x_label, fontsize = 14)
                    plt.tick_params(axis='both', labelsize=14)              #size of axis
                #plt.yscale('log') 
                    plt.grid(True)
                    plt.legend()            
                    plt.savefig(folder_name +'/' +  
                        pre_save_name + '_'  + index[:-4] +'.png', format='png', bbox_inches='tight')
                    #To save plot in folder
        
            plt.close()             #to clsoe the plot not to consume too much resources    
    
    else:           #Wrong case
        print('\n Wrong Nrepl introduced, nothing being done (:')
    
    
    
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
#@@@@@@@@@@@@@@@@@@@@@@@@@@ End ICPMS shit xD @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@q    
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
    
    
    
    
    
    
#%% ###############################################
################### 2) TGA reader ##################### 
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
    
    
 
#%%##################################################
##################### 8) XRD reader #################### 
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
