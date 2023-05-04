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
########### ICPMS excel reader #############
#####################################

def Read_ICPMS_excel (exc_name, D_f_data, cps_sheet_name = 'To_read', 
                      rsd_sheet_name = '%rsd_to_read', Df_correction = False,
                      return_debug = False):
    '''
    Function that will read the excel file from ICPMS and will return df with the relevant
    information, for easier handling /plotting. Note the excel should be a bit preprocessed:
            1) Including sample preparation (neccesary to get the Dilution factor). In this excel,
               the sample names (1st colum where the Df is) must be the same as the names in the top
               of the ICPMS data (do it manually, lazy spaniard, less siesta and more work!)
            2) Clean sheet where only the intensity data is, to load it. THe name must be:
                "To_read". Remove ICPMS blanks (std 0, etc)
            3) Prepare a similar sheet for the %rsd, called "%rsd_to_read" (same dimension as the 
                                                        other sheet, ofc!)
    
    You could use this to get the raw data (output from ICPMS) or to correct them for the
    ICPMS dilution factor.                                                                        
    Maybe that could be automatized? note that requires computing stuff from different sheets,
    and the D_f position could differ from file to file, so maybe more challenging that simply
    working with the excels a bit (Eww)
    
    *Inputs:
        .exc_name: string with the name of the excel, without the .xlsx
        .D_f_data: array with the column number and the row interval in which that data (D_f) is found.
            Note that info is also important for getting the index labels!
            Df_data = [1, 3, 5] means from row 1 to 3 (included) and column 5 (E in letters)
        .cps_sheet_name: string with the name of the sheet with the data to read 
            future maybe also concentration values?). Default value: 'To_read' 
        (from acid vs no acid test)
        .rsd_sheet_name: string with the name of the sheet with the rsd data to read. 
            Default value: '%rsd_to_read'
        .Df_correction: whether you can to correct for the dilution factor in the ICPMS
            correction or not. This correction is compute cps*Df, %rsd*Df, but if you
            only want the raw data, you do not need it. I need to do this also for the 
            return of outputs. Default value = False
        .return_debug: if you want to get some extra df for debug 
        (raw data, without cleaning, so like the excel). Default value = False

        
    *Outputs:
        .several df with the cps, %rsd, std, and cps and std corrected for the 
            dilution factor in the ICPMS sample preparation. Depending whether you want Df
            corrections and/or debug you may have 3,2 or 1 output (always the raw returned)
            
    Note that if you have N outputs, if you want to obtain a variable per output, 
    in the script I should call X variables, like:
            a, b, ..n = Read_ICPMS_excel(name)
    If you write less, say 1, 2, some variable will contain more data, in a dictionary
            
            
    ######## TO DO ########
        . Try-excepts blocks for error in spotting D_f? Maybe relevant when I implement plotting in
            another function
            
    #######################
        '''
    
    
    ########### 1) Raw load ###########
    '''
    Can be done easily with pandas. Since th excel sheet containing the cps and the excel file only
    differs in the .xlsx we can define the excel sheet name with the name given as input:
    '''
    excel_name = exc_name + '.xlsx'
    
    #Load
    Dat_cps = pd.read_excel(excel_name, cps_sheet_name, header = [1])
        #header 1 means take row 1 to give names to the columns
        #That contains the cps and cps*dil factor
        
    Dat_rsd = pd.read_excel(excel_name, rsd_sheet_name, header = [1])
        #This contains the rsd values (measured 3 times, automatically computed average)
    Dat_sa_prep = pd.read_excel(excel_name, 'Sampl_prep', header = None)
                #Sample prep sheet. Ensure it has that name!!!!    
    
    '''
    Note once I suffered that the dimesions of those were not similar, and in one sheet they
    were loadingn NaN values. I just erase those empty stuff in excel (selecting and delete)
    and after it worked!
    '''
    
    
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
    TO apply the Df corrections or not we set a variable to choose, true or false
    '''    
    
    if Df_correction == True:   #compute the Df corrections

        #
        '''
    Now, I can apply the Dilution factor to both the std and the cps, should be straightforward, same
    fashion than above. Well, not as simple, since for the multiplication the indexes should be the same
    so, I redefined (above) the Df indexes so they matched the Df ones, and then that calc is straightforward
        '''
        #
        df_stdDf = pd.DataFrame(df_std.iloc[:,1:] * D_f, 
                  columns=df_std.columns, index=df_std.index)           # std*Df df
        df_cpsDf = pd.DataFrame(df_cps.iloc[:,1:] * D_f, 
                  columns=df_cps.columns, index=df_std.index)           # std*Df df    
        #
        '''
    TO end that, note the isotopes column non is NaN since the multiplications, so lets redefine it again:
        '''
        df_stdDf['Isotopes'] = Dat_cps['Isotopes']              #redefining the isotopes column
        df_cpsDf['Isotopes'] = Dat_cps['Isotopes']
        #
        
        
    
    ############## 4) Ouptuts ######################
    '''
    To make it more compact, let´s group the outputs. Depending on whether you want
    the Df or the debug we return more or less things
    '''

    dfs_raw = {'cps' : df_cps, 'rsd' : df_rsd, 'std' : df_std}
        #dictionary with the raw data (no Df correction applied)
    dfs_debug = {'raw_cps' : Dat_cps, 'raw_rsd' : Dat_rsd }   #Dataframes for debug  ,
        #whihc consist of raw read of excel, without any cleaning         
        
    if Df_correction == True:    #Apply the correction
         Df_dfs = {'cpsDf' : df_cpsDf, 'stdDf' : df_stdDf, 'Df': D_f} #Dictionray containing
        #the df corrected by the dilution factor of ICPMS sample prep   
        
         if return_debug == True:    #return the debug
        #
            return dfs_raw, Df_dfs, dfs_debug               #return of values
    
         else: #no return debug, but yes corrections
            return dfs_raw, Df_dfs 
    else:  #No corrections
        if return_debug == True: #return debug
            return dfs_raw, dfs_debug
        else: #no return debug, neither correectinos
            return dfs_raw
    
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
########### ICPMS Bar plotter #############
#####################################

def ICPMS_Barplotter (df_cps, df_rstd):
    '''
    Function that will do bar plots of the raw data from the ICPMS, the cps and the rstd. By raw
    I mean withoutany corrections/calibrations (neither the dilution factor correction). This is a
    preliminary plot, to check if everything is allright, with the rstd (rstd should be < 20%).
    
    *Inputs:
        .df_cps, df_rstd: dataframes containing the cps and the relative standard deviation. Those are
        outputs for the Read_ICPMS_excel function.
        
    *Outputs:
        .Plots (saving them) of the raw data
    
    
    '''
    
    
    ############# 1) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder a subfolder
    with the relevant elements, to be given, will be created
    '''
    Elem_rel = ['Si28', 'Si29', 'Si30',
            'Al27',
            'Mg24', 'Mg25', 'Mg26',
            'Mn55',
            'Fe56', 'Fe57',
            'Ca42', 'Ca43', 'Ca44', 
            'Na23', 
            'K', 
            'Ti46', 'Ti47', 'Ti48', 'Ti49', 'Ti50',
            'P31', 
            'S32', 'S33', 'S34']      #List of relevant elements
    
    fold_name_Bar = 'Bar_plots'
    path_bar_pl = os.getcwd() + '/' + fold_name_Bar + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)

    #Subfolder with relevant plots:
    path_bar_pl_rel = os.getcwd() + '/' + fold_name_Bar + '/' + 'Relevants' + '/' 
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
    X_axis = np.arange( len(df_cps.axes[1]) -1)                 #To do the 2 bar plot
            #choosing the number of columns. [0] for rows
    b = .4                              #[au] blank space between values <1
    w = (1-b)/2          #bar width

    for i in list( range(4, df_cps.index[-1]) ):     #Loop for the 250 graph plotting
                    #df_cps.index give the index values, low and high
                    ########### Bar plot ###########################
        plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
        plt.title("Concentration of " + df_cps['Isotopes'][i], fontsize=22, wrap=True)           #title
        a = plt.bar(X_axis - w/2, df_cps.loc[i][1:], width = w, edgecolor="black", 
            label = "I ", align='center') 
            #-2 not to plot the blank!! Remove it to plot it!
        plt.ylabel("I [cps]", fontsize=14)              #ylabel
        plt.xlabel('Sample', fontsize = 14)
        plt.tick_params(axis='both', labelsize=14)              #size of axis
        plt.yscale('log') 
        plt.grid(True)
        plt.xticks(X_axis, [ df_cps.columns[1:][j] for j in range(0,len(df_cps.axes[1]) -1 ) ]
            , rotation = 90)
        plt.twinx()             #For setting 2 axes
        aa = plt.bar( X_axis + w/2, df_rstd.loc[i][1:], width = w, edgecolor="black", 
             label = '$\sigma_{rel}$', align='center', color = 'red') 
        plt.ylabel("$\sigma_{rel}$ [%]", fontsize=14)              #ylabel
        #
        aaa = [a, aa]
        plt.legend(aaa, [p_.get_label() for p_ in aaa])
        
        #Saving in the folder
        if df_cps['Isotopes'][i][:-4] in Elem_rel:  #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
            plt.savefig(fold_name_Bar + '/' + 'Relevants' + '/' +
                        'Conc_rsd_' + df_cps['Isotopes'][i] + '.png', format='png', bbox_inches='tight')
            #
        else:        #if the element is not relevant
            plt.savefig(fold_name_Bar +'/' +  
                        'Conc_rsd_' + df_cps['Isotopes'][i] +'.png', format='png', bbox_inches='tight')
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
########### ICPMS plotter #############
#####################################

def ICPMS_Plotter (x, df_cps, x_label, y_label, folder_name = 'Plots' ):
    '''
    Function that will plots of the data from the ICPMS (cps) vs another variable, initially
    time, the cps and the rstd. This assume we have 2 replicates, 1 series after the other
    *Inputs:
        .x: x axis variable in the plot. This should be a df series
        .df_cps, df_rstd: dataframes containing the cps and the relative standard deviation. Those are
        outputs for the Read_ICPMS_excel function.
        .x_label: string that will be the x label for the plot (for math stuff, 
                                    use $$. eg: '$\Delta t[h]$')
        .y_label: string that will be the y label for the plot
        .folder_name: string defining the name of the folder to create to store the plots
            default value: 'Plots'
        
    *Outputs:
        .Plots (saving them) of the x and df_cps data, cps vs x!
    
    ### TO DO: ####
	.Implement error plotting (in an errorbar pyplot)
    '''
    
    
    ############# 1) Folder creation ###############
    '''
    First the folder to store the plots will be created. IN the main folder a subfolder
    with the relevant elements, to be given, will be created
    '''
    Elem_rel = ['Si28', 'Si29', 'Si30',
            'Al27',
            'Mg24', 'Mg25', 'Mg26',
            'Mn55',
            'Fe56', 'Fe57',
            'Ca42', 'Ca43', 'Ca44', 
            'Na23', 
            'K', 
            'Ti46', 'Ti47', 'Ti48', 'Ti49', 'Ti50',
            'P31', 
            'S32', 'S33', 'S34']      #List of relevant elements
    
    fold_name_Bar = 'Plots'
    path_bar_pl = os.getcwd() + '/' + fold_name_Bar + '/'
        #Note os.getcwd() give current directory. With that structure we are able
        #to automatize the plotting!!!
        
    if not os.path.exists(path_bar_pl):
        os.makedirs(path_bar_pl)

    #Subfolder with relevant plots:
    path_bar_pl_rel = os.getcwd() + '/' + fold_name_Bar + '/' + 'Relevants' + '/' 
        #folder path for the relevant plots
    if not os.path.exists(path_bar_pl_rel):
        os.makedirs(path_bar_pl_rel)   
    
    
    ######### 2) plotting ###############
    '''
    This is a loop plot, so beware, will take long (2-3mins!).
    

    '''
    t_start = tr.time()       #[s] start time of the plot execution
    
    ###Plot

    for i in list( range(4, df_cps.index[-1]) ):     #Loop for the 250 graph plotting
                    #df_cps.index give the index values, low and high
		   # 4 because of the way the df is created (and hence the excel tabelle)
        #
        plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
        plt.title("Concentration of " + df_cps['Isotopes'][i], fontsize=22, wrap=True)           #title
        plt.plot(x[:int(len(x)/2)], df_cps.loc[i][1:int(len(x)/2 +1)], 'bo--', MarkerSize = 5, label = 'Repl_1') 
                #+1 needed since the df contain a row with the column names!
        plt.plot(x[int(len(x)/2):], df_cps.loc[i][int(len(x)/2 +1):], 'ro--', MarkerSize = 5, label = 'Repl_2') 
        plt.ylabel(y_label, fontsize=14)              #ylabel
        plt.xlabel(x_label, fontsize = 14)
        plt.tick_params(axis='both', labelsize=14)              #size of axis
        plt.yscale('log') 
        plt.grid(True)
        plt.legend()
        
        #Saving in the folder
        if df_cps['Isotopes'][i][:-4] in Elem_rel:  #if the element is relevant
            #note the -4 is so that that element contain only name and number, like Mg26, not Mg26 (MR),
            #in order to check with the list!
            plt.savefig(fold_name_Bar + '/' + 'Relevants' + '/' +
                        'Conc_' + df_cps['Isotopes'][i] + '.png', format='png', bbox_inches='tight')
            #
        else:        #if the element is not relevant
            plt.savefig(fold_name_Bar +'/' +  
                        'Conc_' + df_cps['Isotopes'][i] +'.png', format='png', bbox_inches='tight')
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
    
    


    
    
#%% ###############################################
################### TGA reader ##################### 
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
##################### XRD reader #################### 
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
