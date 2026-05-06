# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 11:53:09 2025

@author: Administrator

Okay bro, the time has come, time to analyze the hot adsorption experiment!

Note that for the moment I do not do Cleach correction ==> I am not substracting
the blank (sample 1, only BIC and bentonite)! So, for the moment I do not do
anything with that data!!!


KLIAO! No plt title if caption si included (reporting!)



#Tail patttern in Qe for La, Cs, Np, Am, Pu. Hence, removing first 2 samples==>

sample 2,3 gone! This is done by doing [:,3:]. BEfore I was doing [:,1:], only
to remove the procedural blank


La:
    Ce has for 1,2,3,4 less than 0.2ppb!
    
Cs:
    Ce has less than 0.2ppb!
    
    
    
Colours (Gemini defined):
U (Uranium)	Steel Blue	#4682B4	Professional and standard
Pu (Plutonium)	Crimson	#DC143C	Strong contrast against blue
Np (Neptunium)	Forest Green	#228B22	Deep and easy to track
Am (Americium)	Dark Orange	#FF8C00	Bright but readable
Cm (Curium)	Purple	#9932CC	Distinct from the reds and blues

Cs	Slate Gray	#708090	Neutral, distinct from the "colorful" lanthanides.
La	Goldenrod	#DAA520	Warm yellow/gold.
Sm	Lime Green	#32CD32	High visibility; contrasts sharply with Gray and Gold.
Eu	Saddle Brown	 #8B4513 	Dark, earthy, and solid.
Nd	Hot Pink #FF69B4	 Bright and high-energy.
"""



#----------------------------
#%%############################## -1) Cleaning ######################
#----------------------------
from IPython import get_ipython;   
get_ipython().magic('reset -sf')     #Delete all variables


#----------------------------
#%%## ##### 0) General packages ################################
#----------------------------

import time as tr                                #to measure the running time
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
                #Path for personal pc, Ubuntu!
sys.path.insert(0,'C:/Users/Administrator/Desktop/Python/Functions_homemade')
            #Path for guest laptop jrc
import Read_JRC, Fits                                     #my functions
import joblib
from os import path
from datetime import timedelta as td 
        #TO work with time differences in format hours:minute:seconds
from datetime import time   #to work with time data (no differences since hour <24h)

from scipy.stats import pearsonr, spearmanr #statistical correlation tests

#------- Useful variables ----------

Bent_color = {'Sar' : (.68,.24,.31), 'Tur' :  '#F6BE00', 'BK' : 'grey'} 
Font = 18
Markersize = 7
Linewidth = 3

Elem_rel_ben =  ['Si(MR)', 'Al(MR)', 'Mg(MR)', 'Mn(MR)', 'Fe(MR)', 
            'P(MR)', 'Ti(MR)','S(MR)','Sr(LR)', 'Ca(MR)']
                #bentonite relevant elements. No Na since trikcy for ICPMS!
                #Ca better to measure
                #Ojo, Sr here since bentonite leached it!
Elem_rel_CL = ['Cs(LR)', 'Eu(LR)','La(LR)', 'U(LR)', 'Mo(LR)' ] 
                    #no Sr here ni Na!
Elem_rel = Elem_rel_ben + Elem_rel_CL + ['Na(MR)', 'Pu(LR)', 'Np(LR)','Am(LR)', 
            'Cm(LR)', 'Zr(LR)', 'Tc(LR)', 'Nb(LR)', 'Sn(LR)', 'Sb(LR)', 'Nd(LR)']
                #Sn126 half life of 1.9e5y (hides Sb126)
Iso_rel = ['Sr86(LR)','Sr88(LR)', 'Sr90(LR)',  'Cs133(LR)' ,
              #Sr84 0.54wt%, low, Sr87 with inter, so discarded both xD
           'Eu151(LR)', 'Eu153(LR)', 'Sm154(LR)',
'La139(LR)',
 'U234(LR)', 'U235(LR)', 'U236(LR)','Np(237(LR)', 'U238(LR)',  
'Pu239(LR)', 'Pu240(LR)', 'Pu242(LR)',
'Am241(LR)', 'Am243(LR)',
'Cm244(LR)', 'Cm245(LR)', 'Cm246(LR)', 'Cm247(LR)', 'Cm248(LR)']
#plt show after pl save, otherwise not saved properly!



#---------------------------------------
#%%         0.01) Alpha measurements
#----------------------------------------
'''
This was moved to a separate script!
'''





#-------------------------------------------------------
#%%          0.1) ICPMS data processing, running automatization 
#------------------------------------------------

'''
In the 1st part I need to process the ICPMS data. I have automatized part of the process, #
so just run it!

But, also needed to run once, so once its properly done, I just comment it, since the
 important data will be the
ppb data, that I will obtain after the calibration in excel (ewww)
'''

Excel_name = 'Hot_adsorption.xlsx'           #excel with ICPMS remeasurements#
        #I load the main excel, the icpms excel I will not use it from here on!

df_raw = Read_JRC.Read_ICPMS_excel(Excel_name,'250903_S_Ad_Leach_DLA', 
        Is_wash_inside= 1)
df_rstd = Read_JRC.Read_ICPMS_excel(Excel_name,'%rsd', return_debug = False,
             Is_wash_inside= 1)
columns_blks = np.array( [1,2, 9, 10, 17, 42, 67, 71])      #columns where ICPMS 
                            #blanks (and std0) are (excel values)
        #same as for BK Ad U exp bro, perfect copy paste xDDD   
        
# Sr_corr =Read_JRC.ICPMS_Sr_correction(df_raw, df_rstd,
#         Excel_name = 'Sr_correction.xlsx', Sa_start_column = 20, N_sa = 11*4+3)
'''
---Ratio fission produced Sr88/Sr90 --------
Theoretical fission ratio (ORIGEN):  1.04
S Ad Leach 0.1          -2.42
S Ad Leach 0.2           7.07
S Ad Leach 0.3           2.66
S Ad Leach 0.4           1.03
S Ad Leach 0.5           1.33
S Ad Leach 0.6           2.01
S Ad Leach 0.7           1.83
S Ad Leach 0.8           1.49
S Ad Leach 0.9           1.16
S Ad Leach 0.10          1.21
S Ad Leach 0.11          1.08
S Ad Leach 1.1        -113.41
S Ad Leach 1.2         -30.18
S Ad Leach 1.3         -16.33
S Ad Leach 1.4         818.30
S Ad Leach 1.5           1.96
S Ad Leach 1.6           0.51
S Ad Leach 1.7           8.66
S Ad Leach 1.8           2.12
S Ad Leach 1.9           0.26
S Ad Leach 1.10          2.00
S Ad Leach 1.11          1.80
Blank std d             -0.36
std check 2+4 1ppb#     20.51
std check 3+5 1ppb#     -4.05
S Ad Leach 2.1          42.91
S Ad Leach 2.2          27.65
S Ad Leach 2.3          47.63
S Ad Leach 2.4        -239.57
S Ad Leach 2.5         -82.13
S Ad Leach 2.6          14.55
S Ad Leach 2.7          -3.01
S Ad Leach 2.8          -6.19
S Ad Leach 2.9           3.51
S Ad Leach 2.10          1.27
S Ad Leach 2.11          2.96
S Ad Leach 3.1          -8.15
S Ad Leach 3.2         -32.66
S Ad Leach 3.3         -12.58
S Ad Leach 3.4        -444.02
S Ad Leach 3.5           1.44
S Ad Leach 3.6          24.69
S Ad Leach 3.7          -3.32
S Ad Leach 3.8          12.43
S Ad Leach 3.9           6.53
S Ad Leach 3.10          1.80
S Ad Leach 3.11          1.92

Ojo, for 0.1,2,3,6 big ratio Fission Sr88/Sr90, and also high %rsd obtained!

FOr the samples almost never is that ratio preserved, for last samples, 10, 11, almost.
Consider:
    .Sr88 fission possibly wrong, since bentontie also released Sr. Might be that for
    for most concnetrnated MS, bentontie contribution < MS contribution
    .Sr90 computed from Zr91, so not affected

Hence, the ratio should be off, but, when MS contribution > bent contrubution, should
become similar to the theorecial, 1, which is what we kinda see

Am I right?????????

'''
#!!!
#!!!!!!

# df_cps_IS_ICPBlk_corr= Read_JRC.ICPMS_data_process(df_raw, df_rstd, columns_blks, 
#         name_plot_LR_bef= 'IS_sens_LR_bef', name_plot_MR_bef ='IS_sens_MR_bef',
#         name_plot_LR_aft= 'IS_sens_LR_aft', name_plot_MR_aft ='IS_sens_MR_aft',
#         IS_meas = ['Co59(LR)', 'In115(LR)', 'Ho165(LR)', 'Th232(LR)',
#                     'Co59(MR)', 'In115(MR)', 'Ho165(MR)', 'Th232(MR)'],
#         excel_name = 'Corrections_with_Sr90.xlsx')
                #df with the cps corrected to IS and ICPMS Blks :)
'''
The last step is to copy those excel sheets into the main excel, and to calibrate it, 
and then analyze it!

Run it once, and then comment it, thats the way xD

Note really good sensitivities obtained! Maybe having offf the device for some time
helped it?
'''

del (df_raw, df_rstd, columns_blks)   #deleting, to save memory



#----------------------------
#%% ----------------1) Data loading ----------------------------
#----------------------------
'''
d
'''

#The exp data is:  
Exp_df =   pd.read_excel(Excel_name, 'Data_exp', index_col= 0, 
                          header = 48-1, nrows = 80-48) #Exp data
Exp_MS_df = pd.read_excel(Excel_name, 'Data_exp', index_col= 0, 
                          header = 11-1, nrows = 40-11)
                        #data for only the MS
df_ppb = Read_JRC.Read_ICPMS_excel(Excel_name, 'To_read_ppb',
                    return_debug = False)               #ppb data  
df_ppb_std = Read_JRC.Read_ICPMS_excel(Excel_name, 'To_read_ppb_std',
                                     return_debug = False) #[ppb] std data
               
df_rsd = df_ppb_std/df_ppb *100         #For debug! Something going on here!!

df_sens = Read_JRC.ICPMS_Sens_finder(
    Excel_name, Excel_sheet = 'Calib', Sens_column = 12)    #cps/ppb
                  
############ pH data
pH = Exp_df.loc['pH']
pH_MS = Exp_MS_df.loc['pH'][:-3]
        #[:-3] to remove the NF values, that I did not prepare at the end
Delta_pH = Exp_df.loc['Delta pH']
Delta_pH_MS = Exp_MS_df.loc['Delta pH'][:-3]


N_sa = round(pH.shape[0]/3)        #number f samples per repl


#--------- Printing replicates which had more than 10d of shaking
'''
Since this exp was peformed while hot cells special work was performed, I was 
not always allowed in the Lab. Hence, some samples shaking time are not 10d,
but more. It could be good to have them in mind at least
'''

# print('---------------------------------------------------------------------')  
# print('---- Here the list with samples whose shaking time was > 1.1*240h = 264h --')

# for i, value in enumerate(Exp_df.loc['Delta_t [h]']):
#     if value > 1.1 * 240:   #printing when time > 1.1*240 = 260h!
#         print(Exp_df.loc['Delta_t [h]'].index[i][10:], ', Delta_t[h]:', 
#               np.round(value) )     #print index and value, rounded
#                 #[10:] not to print the initial name, S Ad Leach
# print('---------------------------------------------------------------------')        
'''
Many samples:
 1.3 , Delta_t[h]: 309.0
 1.5 , Delta_t[h]: 313.0
 1.6 , Delta_t[h]: 309.0
 1.7 , Delta_t[h]: 264.0
 1.9 , Delta_t[h]: 307.0
 2.2 , Delta_t[h]: 335.0
 2.3 , Delta_t[h]: 335.0
 2.4 , Delta_t[h]: 330.0
 2.5 , Delta_t[h]: 313.0
 2.7 , Delta_t[h]: 264.0
 2.9 , Delta_t[h]: 288.0
 3.2 , Delta_t[h]: 335.0
 3.3 , Delta_t[h]: 271.0
 3.5 , Delta_t[h]: 309.0
 3.6 , Delta_t[h]: 306.0
 3.8 , Delta_t[h]: 307.0
 3.9 , Delta_t[h]: 306.0
 
Bear that in mind kliao!
'''

#----------------------------
#%% ---------------- 2) Data cleaning ----------------------------
#----------------------------
'''
Main things:
    -Sr90 issue
    -Df correction
    -to M
    -elemental
    -Cleach (otpional, discarded)

Note here for the 1st time, the MS has also 2 Dilution factors!

'''


#----------- 2.1) Sr90 issue ---------
'''
Well, since in the funciton I did replcae negative values with 0 (NaN, read as 0),
not to missunderstand the data (imagine, MS has Sr90, but sample har <0, I put 0, so
                                I could think the sorption was 100%). I will
replace here in python 0 with NaN!
'''

df_ppb.loc['Sr90(LR)'].replace(0, np.nan,  inplace = True)  #Replace 0 (neg values) by NaN



#%% --------- 2.1.5) SNF leach correct

'''
I could also apply the SNF leach corrections:
    .Cs
    .Pu
    
I will only do Pu so far, since i did not do the Xe corrections for the Cs, since
I saw in the SNF charact that after Cs removal little Cs137 present, and not trustworthy
with ICPMS!

This is better to do before Df corrections since nefore MS and samples are together, so
I would need to apply this only once, and not twice (to MS; and then to samples)

'''

Dict_SNFcorr = Read_JRC.ICPMS_SNF_Leach_correction(df_ppb, df_ppb_std, 
                    df_sens, Correction_to_do= 'Pu')




################## 2.2) Df spotting #################
Df_exp = (Exp_df.loc['Dil fact']).apply(pd.to_numeric) 
            #Dilution factor of the bottle to use to do the ICPMS samp prep. 
                #making that a numeric data, not string 
Df_MS_exp = (Exp_MS_df.loc['Dil fact'][:-3]).apply(pd.to_numeric) 

Df_ICPMS_1 = Read_JRC.ICPMS_Df_finder (Excel_name, [56,66],
                                     samp_prep_sheet_name = 'Sampl_prep') 
            #Df for the ICPMS sample prep (around 50)

Df_ICPMS_2 = Read_JRC.ICPMS_Df_finder (Excel_name, [81, 91],
                                     samp_prep_sheet_name = 'Sampl_prep') 
            #2nd dilution. Note I need all the names!!
            
# print('\n ### Dilution factors ######\n')
# print('\n 1st ICPMS dilution: ')
# print(Df_ICPMS_1)
# print('\n 2nd ICPMS dilution: (NaN = no 2nd dilution made) ')
# print(Df_ICPMS_2)


'''
So, Spotting them is easy. NOw I need to combine them somehow, like 
ignoring names, 
since they have different. Or just modify the name so they have 
the same name! I did that!
 Now, to be able to multiply them, since I have
 blank data, which is 1 in the case of the dil2 (not all the 
 samples have dilutions 2), 
 we will substitute nan for 1, and then multiplying!
'''

Df_ICPMS_2.fillna(1, inplace = True)            #removing nan for 1

Df_ICPMS = Df_ICPMS_1 * Df_ICPMS_2      #getting total Df

#That works, but that variable has nan columns and row, we need to 
#remove them. I can do it like:
Df_ICPMS = Df_ICPMS[pd.notnull(Df_ICPMS.index)]
#Df_ICPMS.dropna(inplace = True)   #this remove all rows where there
# is NaN, so remove too much!

del Df_ICPMS_1, Df_ICPMS_2                 #deleting to save some memory



#--------------- 2.3) Sepparating MS and samples ----------------------
#This is trivial:
    
#Using the Dict SNF corr better!
df_sa = Dict_SNFcorr['Data']['dat'].drop(['S Ad Leach 0.1', 'S Ad Leach 0.2',
           'S Ad Leach 0.3', 'S Ad Leach 0.4','S Ad Leach 0.5',
            'S Ad Leach 0.6','S Ad Leach 0.7','S Ad Leach 0.8',
            'S Ad Leach 0.9','S Ad Leach 0.10', 'S Ad Leach 0.11'], axis = 1)

df_sa_std = Dict_SNFcorr['Data']['std'].drop(['S Ad Leach 0.1', 'S Ad Leach 0.2',
           'S Ad Leach 0.3', 'S Ad Leach 0.4','S Ad Leach 0.5',
            'S Ad Leach 0.6','S Ad Leach 0.7','S Ad Leach 0.8',
            'S Ad Leach 0.9','S Ad Leach 0.10', 'S Ad Leach 0.11'], axis = 1)

df_MS = Dict_SNFcorr['Data']['dat'][['S Ad Leach 0.1', 'S Ad Leach 0.2',
           'S Ad Leach 0.3', 'S Ad Leach 0.4','S Ad Leach 0.5',
            'S Ad Leach 0.6','S Ad Leach 0.7','S Ad Leach 0.8',
            'S Ad Leach 0.9','S Ad Leach 0.10', 'S Ad Leach 0.11']]     #df with MS

df_MS_std = Dict_SNFcorr['Data']['std'][['S Ad Leach 0.1', 'S Ad Leach 0.2',
           'S Ad Leach 0.3', 'S Ad Leach 0.4','S Ad Leach 0.5',
            'S Ad Leach 0.6','S Ad Leach 0.7','S Ad Leach 0.8',
            'S Ad Leach 0.9','S Ad Leach 0.10', 'S Ad Leach 0.11']]


df_sa_rsd = df_sa_std / df_sa *100



############# 2.4) Applying DF corrections! #####################
"""
last step is applying the Df corrections! remember names are
 crutial here!!!'
"""
df_MS_Df = Read_JRC.ICPMS_Df_corrector(df_MS, 
            Df_ICPMS[ ['S Ad Leach 0.1', 'S Ad Leach 0.2',
                       'S Ad Leach 0.3', 'S Ad Leach 0.4','S Ad Leach 0.5',
                        'S Ad Leach 0.6','S Ad Leach 0.7','S Ad Leach 0.8',
                        'S Ad Leach 0.9','S Ad Leach 0.10', 'S Ad Leach 0.11'] ] )

df_MS_std_Df  = Read_JRC.ICPMS_Df_corrector(df_MS_std, 
            Df_ICPMS[['S Ad Leach 0.1', 'S Ad Leach 0.2',
                       'S Ad Leach 0.3', 'S Ad Leach 0.4','S Ad Leach 0.5',
                        'S Ad Leach 0.6','S Ad Leach 0.7','S Ad Leach 0.8',
                        'S Ad Leach 0.9','S Ad Leach 0.10', 'S Ad Leach 0.11']])


df_sa_Df = Read_JRC.ICPMS_Df_corrector(df_sa, 
            Df_ICPMS.drop(['S Ad Leach 0.1', 'S Ad Leach 0.2',
                       'S Ad Leach 0.3', 'S Ad Leach 0.4','S Ad Leach 0.5',
                        'S Ad Leach 0.6','S Ad Leach 0.7','S Ad Leach 0.8',
                        'S Ad Leach 0.9','S Ad Leach 0.10', 'S Ad Leach 0.11']) )

df_sa_std_Df = Read_JRC.ICPMS_Df_corrector(df_sa_std, 
            Df_ICPMS.drop(['S Ad Leach 0.1', 'S Ad Leach 0.2',
                       'S Ad Leach 0.3', 'S Ad Leach 0.4','S Ad Leach 0.5',
                        'S Ad Leach 0.6','S Ad Leach 0.7','S Ad Leach 0.8',
                        'S Ad Leach 0.9','S Ad Leach 0.10', 'S Ad Leach 0.11'] ) )


df_sa_Df2 = Read_JRC.ICPMS_Df_corrector(df_sa_Df, Df_exp)
df_sa_std_Df2 = Read_JRC.ICPMS_Df_corrector(df_sa_std_Df, Df_exp)

df_MS_Df2 = Read_JRC.ICPMS_Df_corrector(df_MS_Df, Df_MS_exp)

df_MS_std_Df2 = Read_JRC.ICPMS_Df_corrector(df_MS_std_Df, Df_MS_exp)


#rsd calc
df_sa_rsd =df_sa_std_Df2 / df_sa_Df2 *100

df_MS_rsd = df_MS_std_Df2 / df_MS_Df2 *100


del (df_MS, df_MS_std, df_sa, df_sa_std, 
     df_sa_Df, df_sa_std_Df, df_MS_Df, df_MS_std_Df ) #deleting to remove unneded var




#%% ---------- 2.4) Conversion to M ---------
#Conversion to M (SI unit) needed, in order to compare with literature!

'''
Ojo, since here I have radioactive elements, I need to expand the excel of M,
including them!

Ojo, for the MS, V I know since I prepare it accurately. Masses I measured mostly
for all. Note for 0.11 I did not. The function will do the average of the other
masses. So far I will consdier that value. I do have the density from the SNF
Charact 4, but should be really similar to the obtained here, so I keep this

Beware!!!!!! That is an ASSUMPTION!!!!!!!
!!!!!!!
'''

V_MS = 165          #[mL] theoretical volume of the MS


#Note 0s and NaN could be there. The function ppb to M will take care of that!


#Samples
Dict_M = Read_JRC.ICPMS_ppb_to_M(df_sa_Df2,df_sa_std_Df2, 
            m_s = Exp_df.loc['m_leach [g]'], Delta_m_s= 0.001, 
            V_s = Exp_df.loc['V_leach [mL]'] * 10**-3 ) #samples
        #delta m_s higher since less precision in a GB!
#MS
Dict_M_MS = Read_JRC.ICPMS_ppb_to_M(df_MS_Df2, df_MS_std_Df2, 
            m_s = Exp_MS_df.loc['Total mass of solution (g)'][:-3], Delta_m_s= 0.001, 
            V_s = V_MS * 10**-3)


"""
Ojo, for the MS I have no rho, so they are NaN!! I need to take care on that.
Actually, for no MS I have those data, so I must guess it, or assume it :/

Note that I could measure this data when I am back, but right now not :/
"""


#del (df_sa_Df2,df_sa_std_Df2, df_MS_Df2, df_MS_std_Df2)   
    #deleting unncesary variables
    


#%% ------------ 2.5) Elemental conversion ------------
''''I need to conver both MS data and sample data!

Doing like with a dictionary keep a compact structure, less mess of variables :D
I think I will try to do more of this!

'''

# Dict_el = Read_JRC.ICPMS_Isot_to_Elem_from_dict(
#     {'dat': Dict_M['dat'], 'dat_std':Dict_M['dat_std'], 'MS': Dict_M_MS['dat'], 
#      'MS_std': Dict_M_MS['dat_std']} )      #to element, 
Dict_el_sa = Read_JRC.ICPMS_Isot_to_Elem( Dict_M['dat'], Dict_M['dat_std'], 
                                         Radiaoct_here= True )
Dict_el_MS = Read_JRC.ICPMS_Isot_to_Elem( Dict_M_MS['dat'], Dict_M_MS['dat_std'], 
                                         Radiaoct_here= True )
            #.loc[['Cm244(LR)','Cm245(LR)','Cm248(LR)']]
#Here already for Cm in sa than in MS, but not in Dict_M==> wtf?? was Cm248, Th232 interf!!
#Still Cm247 is not soo okay neither, see the Iso Qe plots



#%% --------------- 2.6) MS Plotting ----------
'''
Okay bro, how do I do the plotting of the MS? COuld that be done?
Maybe a massive plot with all the reelvant elements?

From the df, if we click to sort by values, we see the order is, from inc to 
decre:
    Na, Al, Si, Mg, Mo, Zr, Ca, S, Zr, P, Ba, Zn, U, Fe, Mn, Cu, Ni, Cr, Ti,
    Br, Cs, Se, Co, Sr, Cd, Sn, Pd, In, Ho, Th, V, La, Ga, Ce, Pr, Ge, Sm, 
    Ag, Cm, Am, Te, Pu, Ru, W, Sn, Rb, Sc, Gd, Pt, Pu, Eu, Hf, Np, Bi, 
    Rh, Xe, Dy, Nb, Ta, Au, Tb, Yb, Er, TL, Re, Lu, Tm
    (for some, I have 2, MR and LR, writing only the first one xD)
    
Bro, what a mess.  Classifying them a bit:
Note that those are:
    
    ACtinides: U, Pu, Np, Am, Cm
    Fuel specific things: Zr we see!, Cs, Ba, La, Pb, Cm, (Am, Pu, Np negligible)
    REsin: Mo! P!
    Fission Products?
    

Okay, I will try to do a massive bar plot, at least for the MS 11, and maybe
later for MS 2, MS 1, or MS 6
'''

#Dict_el_MS['dat']
#Dict_el_MS['dat']['S Ad Leach 0.11']


X_axis = np.arange(Dict_el_MS['dat'].shape[0])

## We can plot all the mS doing a for loop iterating through the columms:

    
#------------- I comment this since it takes a lot!!!
#
# for i in list( range(0, Dict_el_MS['dat'].shape[1]) ):    #loop thorught the columns (MS)
# #-- Massive bar plot for MS 11
#     plt.figure(figsize=(25,10))  #width, heigh 6.4*4.8 inches by default
#     plt.title("Elemental concentration of MS" + str(i+1), fontsize=Font)    #title  
#             #note index 0 is for  sample 1!
#     plt.bar(X_axis, Dict_el_MS['dat'].iloc[:,i], edgecolor="black",
#         yerr = Dict_el_MS['std'].iloc[:,i], width = 0.9) 
#     plt.xlabel("Element", fontsize=Font)                              #xlabel
#     plt.ylabel("Conc [M]", fontsize=Font)                          #ylabel
#     plt.yscale('log')                                    #y axis in log scale
# # Set size of tick labels.
#     plt.tick_params(axis='both', labelsize=Font)              #size of axis
#     plt.xticks(X_axis, Dict_el_MS['dat'].index, rotation = 90)
#     plt.minorticks_on()             #enabling minor grid lines
#     plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        
#             #which both to plot major and minor grid lines
#     plt.grid(which = 'major')
#     plt.savefig('Conc_MS_' + str(i+1) + '.png', format='png', bbox_inches='tight') 
#     plt.show()

'''
Lot of pltos, with lots of bar, but the most complete plots you can get xD

Most abundant elements are Na (MR), followed by Si(MR). Then, I have elements
surrounding those masses in MR. After, LR:
    .Mo, Zr, Ba, U, Cs, 
    

Ojo! More U than Cs!! But still most abundant Mo! From teh filtration? I must
have a look at SNF leach script to find this, I will do those plot as well!

I will possibly average replicates!
'''


#[j[:-4] for j in Dict_el_MS['dat'].index]



#%% ############ 2.6) Removal of the nuclides leached by bentonite ############
#18/2/25
'''
Okay my dear friend, its time to removed the elements leached by bentonite,
since considedring the Sr case, is evident we need to do it in order to get
an adsorption isotherm. 
We need to do this in a replicate by replicate case, since like this Qe/Kd
are computed

I have 3 replicates: 1 with 1_1, 1_2, ... (1 is blank)
                    2 with 2_1, 2_2, 2_3,..
                    3 with 3_1, 3_2, 3_3,...
                    
MS several, 0_1, 0_2, ....
0_1 is blank, BIC water

I need, for each replicate, to substract the effect of the bentonite, which 
would be 1- 0_1 for each replicate, 
        1_1-0_1 for 1
        2_1-0_1 for 2
        3_1 - 0_1 for 3
        
And then, apply that correction to each replicate from sample 2 on:
    1_2 - (1_1-0_1), 1_3 - (''), ...
    2_2 - (2_1-0_2), 2_3 - (''),...
    ...

I did that with a function, have a look at it if u need :)
'''

###Substraction correction
# df_ppb_brA,df_ppb_std_brA = Read_JRC.ICPMS_Removal_Bent_leach(
#     df_ppb_sa_Df2, df_ppb_sa_std_Df2, df_MS_Df, df_MS_std_Df, 
#     return_leached= 0)


#Ratio based correction
Dict_bent_rem = Read_JRC.ICPMS_Removal_Bent_leach_ratio(
    {'ppb': Dict_el_sa['dat'], 'ppb_std': Dict_el_sa['std'],
     'MS': Dict_el_MS['dat'],'MS_std': Dict_el_MS['std']}, return_leached= 1, 
    Nucl_rel = ['U(LR)', 'Cs(LR)', 'Sr(LR)', 'S(MR)'] )
            #Shall I put other element??



#No Cleach correction!!
Dict_bent_rem['dat_br'] = Dict_el_sa['dat']    #manually overwrite the relevant variables
Dict_bent_rem['std_br'] = Dict_el_sa['std']

print('---------- Kliao, Achtung ---------------------')
print('--------------------------------------------------')
print('No Cleach removal applied, beware !!!!!!')
print('--------------------------------------------------')
print('---------- Kliao, Achtung ---------------------\n')



#%% ------------ 2.7) Outliers finding and removal ----------------------



# Read_JRC.ICPMS_Plotter_mean_blk_N (
#         x_list = [Dict_el_MS['dat'].iloc[:,3:] ,
#                   Dict_el_MS['dat'].iloc[:,3:] ,
#                   Dict_el_MS['dat'].iloc[:,3:]   ], 
#         std_x_list = [Dict_el_MS['std'].iloc[:,3:] ,
#         Dict_el_MS['std'].iloc[:,3:], Dict_el_MS['std'].iloc[:,3:] ],
#         y_list = [ Dict_el_sa['dat'].iloc[:,1:N_sa],
#                   Dict_el_sa['dat'].iloc[:,N_sa+1: 2*N_sa],
#                   Dict_el_sa['dat'].iloc[:,2*N_sa+1:]],
#         std_y_list =[ Dict_el_sa['std'].iloc[:,1:N_sa],
#         Dict_el_sa['std'].iloc[:,N_sa+1:2*N_sa],
#         Dict_el_sa['std'].iloc[:,2*N_sa+1:] ],
#     element_index = Dict_el_sa['std'].index,
#     x_label = '$C_0$ [M]', y_label ="$C_e [M]$",
#     labels = ['Repl 1', 'Repl 2', 'Repl 3'],
#     folder_name = 'Conc_i_vs_f_repl', pre_title_plt = "Initial vs final concentration ", 
#     pre_save_name = 'Con_i_vs_f', Logscale = 1, Nucl_rel= Elem_rel, plot_everything= 1)

'''
A lot of pltos, 70s to make them, but with them I will spot the outliers!

Spotting relevant folder:
    3_1 seem outlier for Cs, Mo, Np (U could be), Zr, Cd, Ce, Nd, Pr, Ru, Sm, Tc
    1_7 for Cs, Mo, Np, U, Zr, Am, Pu, Cd, Ce, Nd, Pd, Ru, Sm, Tc

Okay, consistent enought, I will remove those. Others sample number not so clear,
so I will do that!
'''


Dict_el_sa['std']['S Ad Leach 3.1'] = (Dict_el_sa['std']['S Ad Leach 1.1']+
                                      Dict_el_sa['std']['S Ad Leach 2.1'])/2
Dict_el_sa['std']['S Ad Leach 1.7'] = (Dict_el_sa['std']['S Ad Leach 2.7']+
                                      Dict_el_sa['std']['S Ad Leach 3.7'])/2
Dict_el_sa['dat']['S Ad Leach 3.1'] = (Dict_el_sa['dat']['S Ad Leach 1.1']+
                                      Dict_el_sa['dat']['S Ad Leach 2.1'])/2
Dict_el_sa['dat']['S Ad Leach 1.7'] = (Dict_el_sa['dat']['S Ad Leach 2.7']+
                                      Dict_el_sa['dat']['S Ad Leach 3.7'])/2
pH['S Ad Leach 3.1'] = (pH['S Ad Leach 1.1']+
                                      pH['S Ad Leach 2.1'])/2
pH['S Ad Leach 1.7'] = (pH['S Ad Leach 2.7']+
                                      pH['S Ad Leach 3.7'])/2
Delta_pH['S Ad Leach 3.1'] = (Delta_pH['S Ad Leach 1.1']+
                                      Delta_pH['S Ad Leach 2.1'])/2
Delta_pH['S Ad Leach 1.7'] = (Delta_pH['S Ad Leach 2.7']+
                                      Delta_pH['S Ad Leach 3.7'])/2

print('----------------------------------')
print('Outliers removed: 1_7, 3_1')
print('----------------------------------\n')



#----------------------------------------------------------------
#%% ------------------- 3) pH plot --------------
#----------------------------------------------------------------

'''
Okay, so we could plot the pH values.Getting them is easy. 
We could do as a 1st plot a plot with 4 series, the 3 replicates + mother 
sol!

Note I shuold use the main variables, pH and DeltapH, since there I corrected
for the outliers!
'''
    #13,10'pH for the T Ads U experiment. Assigned $\Delta(pH[mV]) = 5mV$'
plt.figure(figsize=(13,10))          #width, heigh 6.4*4.8 inches by default (11,8)
                #I need to enlarge since because of the title is big, I think
#plt.title('pH for the S Ad Cleach exp', fontsize=22,
#          wrap=True, loc = 'center')     #title
plt.errorbar( np.array([x for x in range(1,N_sa+1) ]), pH_MS.values,
             Delta_pH_MS.values, 0,
             'o--', markersize = Markersize, label = 'Mother sols', elinewidth = Linewidth)  
plt.errorbar( np.array([x for x in range(1,N_sa+1) ]), pH[:N_sa],
             Delta_pH[:N_sa] , 0,
             '*', markersize = Markersize, label = 'Repl 1', elinewidth = Linewidth) 
plt.errorbar( np.array([x for x in range(1,N_sa+1) ]), pH[N_sa:2* N_sa],
             Delta_pH[N_sa:2* N_sa] , 0,
             's', markersize = Markersize, label = 'Repl 2', elinewidth = Linewidth) 
plt.errorbar( np.array([x for x in range(1,N_sa+1) ]), pH[2*N_sa:] ,
             Delta_pH[2*N_sa:] , 0,'^', #color = Bent_color['Sar'], 
             markersize = Markersize, label = 'Repl 3', elinewidth = Linewidth) 
             #'^', markersize = Markersize, label = 'Repl 3',) 
                    #Like that you can plot the blank
plt.xlabel(" Sample / mother sol", fontsize=18)              #ylabel
plt.ylabel('pH', fontsize = Font)
plt.tick_params(axis='both', labelsize=18)              #size of axis
plt.minorticks_on()             #enabling minor grid lines
plt.grid(which = 'minor', linestyle=':', linewidth=0.5)
plt.grid(which = 'major')
plt.legend(fontsize = Font)
plt.savefig('pH_SAdLeach.png', format='png', bbox_inches='tight')   
            #To save plot in folder
plt.show()  
''' Analysis

.MS fairly constant pJ betwen 7.8 and 8. The pattern indicates that the Leachate
was a bit more acid than the BIC water

.Samples flcutate moll between 8.2 and 8.8. Still overall slight basification
obtained!
Why the fluctuations? No clue, but they are relatively small (but greater than 
        the assigned pH uncertainty)

'''

# --------------------------
# %% ----------- 3.5) Dose measurement dosimeter  ------
# ---------------
'''
Okay, I measured the dose with an albedo soimeters of:
    .MS
    .Solid
    .Superantant (liquid)
    
I measure it thourhg the GB, touching the container of those things with the detector.
Touching laterl for liquid, for solid bottom, where solid was.

Max values were written, went relatively constant.

Note that the volumes of liquids differ, I should normalize them to the volume
Also the solid masses could vary, I should normalize for the masses, which I do not
have tought xD

Example: solid 1.11 les dose than 2.11, 3.11, but I used 1.11 for XRD, FTIR, so I have
less solid there. POssibly when normalizing they shoudl be somehwat similar

I can compute the ratio Dosis solid/liquid to begin with:
'''

Dose_ratio = Exp_df.loc['Dose solid [uSv/h]'] / Exp_df.loc['Dose Liquid [uSv/h]']
            #Dose ratio, solid/liquid

'''
Thats always >1, from 3 to 8. Where not, and is <1, is because the dose levels
were comparable to background ones (or lower).

the dose emasurement was also flucutating, and tough to do, since depending on how
to put the bottle, the glove, the solid inside the bottle, etc. But with this we
get a sense that solid got more dose than the liquid. Sadly can not be compared to the
initial dose. The reason I do not truly understand, since I thought that it should be
related to the container, but I put MS 11 in the same containter as the solid 11, but
still different values. 

And the used bottles didnt have dose, so the wall sorption hypothesis should not
happen.

'''



#----------------------------------------------------------------
#%% ################# 4) Kd and Qe calc ################################
#----------------------------------------------------------------
'''
I created a function, so should be trivial xD

Ojo, I encountered here that the data from Exp_df, mass and so, were not
numeric. I incldued in the funciton the apply(pd.to_numeric), but this was not
the issue, the issue is having no numbers in the excel, but strings. This
is result from the copy only numbers, so take care, copy formulas also!!!!!!!
(Linux)
'''


#isotp version
Dict_ads_iso = Read_JRC.ICPMS_KdQe_calc_Ad (
    df_MS = Dict_M_MS['dat'], df_MS_std= Dict_M_MS['dat_std'],
    df_dat = Dict_M[ 'dat'], df_dat_std= Dict_M[ 'dat_std'],
    df_VoM_disol= Exp_df.loc['V_leach [mL]']*10**-3, 
    df_m_be= Exp_df.loc['mass of bent [g]']*10**-3, ret_Co__Ceq= True)
'''
Beware bro, Sr90 is tricky;

Since I am replacing negative values with 0, if MS has not 0, then the sorption
would be assumed to be 100%, so I must change that!!
'''

Dict_Qe_iso = Read_JRC.ICPMS_MeanStd_calculator(Dict_ads_iso['Qe'].loc['Sr90(LR)'], Dict_ads_iso['Qe_std'].loc['Sr90(LR)'],
                                            Nrepl = 3)

#### M version
Dict_ads = Read_JRC.ICPMS_KdQe_calc_Ad (
    df_MS = Dict_el_MS['dat'], df_MS_std= Dict_el_MS['std'],
    df_dat = Dict_bent_rem[ 'dat_br'], df_dat_std= Dict_bent_rem[ 'std_br'],
    df_VoM_disol= Exp_df.loc['V_leach [mL]']*10**-3, 
    df_m_be= Exp_df.loc['mass of bent [g]']*10**-3, ret_Co__Ceq= True)
'''
Bro, wtf, what a chaos

+*Opening rsq, and clicking on a random sample, to sort it from higher to lower
, we see:
    .highest sorption for Am (99%), Pu, nd, Pu, Cs, La, U, Bi (60%)
    .Clicking on 2_11: highest for Al, Ga, Mn, Pu (88%), Cs, Ni,Bi, Nb, Cr, Pb, 
    Am (58%)
    
Essentially, a mess. Lets see the replciates, and if they kinda agree, we do
the mean!
'''


#--------------------------------------------
#%%                     6) Activity stimation 
#--------------------------------------------
'''
Okay, in order to do some radioactivity measurement:
    .Gamma (Joseph)
    .Alpha (Ana Sanchez/Ale Kindabum)
    
I need to give them some stimats of the dose I should have. Note I would be intersted
in measuring:
    .MS
    .Superantants
    .Solid
    
In theorey., MS = Solid + supernatant

TO get a sense, I could measure the dose of
    .1, 2, 6, 11 (1 as blank)

So, measurenets:
    *Liquid: MS 1, MS 2, MS 6, MS 11 (4), 1.1, 2.1, 3.1, 1.2, 2.2, 3.2, 1.6, 2.6, 3.6, 
        1.11, 2.11, 3.11 (12 or 4 if 1 replicate only) ==> 16 samples 
        (if all replicates measured, or 6 if only 1 repl)
    *Solid: 1.1, 2.1, 3.1, 1.2, 2.2, 3.2, 1.6, 2.6, 3.6, 
        1.11, 2.11, 3.11 (12 samples if all replicates measured, or 4 if one replicate only)

That would need some time to measure touhght. Fuck, I could have order it already xD


#--- Copy from SNF cahract script ------
To copmute the expected dose of the samples, in order to predict for the
gamma measurements, the specific activities will be used.

.The Gamma in B123 could measure from 0.1Bq to 2000Bq. 1Bq = 1 desin/s. 
We will set as target 2k s-1.
.V of the vial = 1.8mL.

Half lifes from nndc3!

With the half life, the specific activty could be computed:
    a = N_A * ln (2) /(T1/2 * M)
    
And with the specific activity the toal activty (A) could eb computed:
    a = A/m ==> A = m*a
    
Note I do not have mass, rather mass/total mass. SO I can get activity/total mass!

I created a funciton to do this, so lets fucking go! Ojo, the function requires
the data in ppb mf! Or, you modify it!

'''

# ------------- 9.1) Alpha activity -------------

'I put all the relevant elements: Am mostly, Pu, Np, U. In the function, so lets go'


Act_alp_MS = Read_JRC.ICPMS_Get_Activity(Dict_M_MS['dat'], Dict_M_MS['dat_std'], 
                            Type= 'Alpha', Conc_units = 'M' )   #Bq/L of MS

Act_alp_S = Read_JRC.ICPMS_Get_Activity(Dict_M['dat'], 
    Dict_M['dat_std'],  Type= 'Alpha', Conc_units = 'M' )   #Bq/L of Samples
'''
The results are too low, the total sunm is mu Bq, from Am241, or 10-5 if Cm244 in,
which possibly not there. THen maybe will not
be measured. I will ask anyhow, about the max limit. Note Max conc, for MS 11, is:
    
    10-8 M Am241, Am243
    10-9 M Pu239,240,242
    10-9 to -11M Cm244,248
    10-8 to -6 M for U234 to 238
    
So low concentrations, and hence low activity.

'''

'''
Okay, Alpha measurements were performed, see the script. Obtaining 2 peaks:
    peak 1 at 5.450 MeV
    peak 2 at 5.75 MeV, more intense than peak 1

The alpha emitors in that region are:
    Nuclide   E (MeV)
    Pu238---> 5.495 
    Am241--- 5.485
    Am243 -- 5.275 
    Cm243 -- 5.785
    Cm244 -- 5.804
    
Hence, peak 2 should be Cm243, Cm244
peak 1 should be Pu238, Am243, Am241

If I compute the expected activty of the samples, do i get the measured ones,
considering that? Note that Cm243 I do not have, I measured both Am243-Cm243.

Okay, I computed the Bq/L from the measurement. Does it match the predicted
Bq/L from the ICPMS?

To compare them, since I have the efficiency factor, I should scale the ICPMS
ones down, since not all the counts were measured. 


'''
ef_det = 13/100             #ef detector alpha 2, the one used
A0 = Act_alp_MS['Act [Bq/L]'].loc['A_tot(LR)']*1e-3 * ef_det     
        #Total activity [kBq/L] corrected for the detectrs
A0_std = Act_alp_MS['Delta[Act[Bq/L]]'].loc['A_tot(LR)']*1e-3 * ef_det     
Af = Act_alp_S['Act [Bq/L]'].loc['A_tot(LR)']*1e-3 * ef_det
        #Total activity 
Af_std = Act_alp_S['Delta[Act[Bq/L]]'].loc['A_tot(LR)']*1e-3 * ef_det   

#Lets gather that in a df better:

A0Af_det_corr_df = pd.DataFrame({'A0': A0, 'std_A0': A0_std,
                        'Af': Af, 'std_Af': Af_std})
        #not the best df, but nice to have a look and write those

#I did a function to compute the Asq from the ICPMS:
    
Asq = Read_JRC.ICPMS_Asq_calc (Act_alp_MS['Act [Bq/L]'], 
    Act_alp_MS['Delta[Act[Bq/L]]'], Act_alp_S['Act [Bq/L]'], 
        Act_alp_S['Delta[Act[Bq/L]]'])

'''

Computed from measurements, total activity:
    Sample      kBq/L from meas         kBq/L from ICPMS pred (see calcs below)
    0.2       2.55                          3.52
    0.6         22                          21
    0.11        540                         580
    
\
Brooo, thats a crazy good agreement, wow!!!
!!!
!!!

Ojo, in A_tot I do an error, since I am adding both the (LR) and (MR). For Asq
is not an issue, but for absolutes magnitudes might be. Still, since Cm was
not measure for MR, the error is not fatal, but bear in mind!!! 
Correct it better mF !!!
!!!!
!!!

Ojo!!!

I obtained alpha activity removal of 40% for sapmle 11, while from measurements
I obtained 20%, factor 2 difference ==>
ICPMS predict to remove more dose than what the spectrometry said

Note that 40% was the rsq of Cm, major
alpha contributor. So, you need to check the alpha data, what happened there?
??????
??




'''

# ------------- 9.1) Gamma dose -------------

Act_gam_MS = Read_JRC.ICPMS_Get_Activity(Dict_M_MS['dat'], Dict_M_MS['dat_std'], 
                                     Type= 'Gamma', Conc_units = 'M')   #Bq/L


#More work to be done....


#qwdsafewsaacdsfedwS         #Work under progress...

#-----------------------------------------------
#%% --------- Plolt replicates -----------------
#-----------------


####### rsq plot
# Read_JRC.ICPMS_Plotter_mean_blk_N (
#         x_list = [np.log10(Dict_el_MS['dat'].iloc[:,3:] ),
#                   np.log10(Dict_el_MS['dat'].iloc[:,3:] ),
#                   np.log10(Dict_el_MS['dat'].iloc[:,3:] )  ], 
#         std_x_list = [Dict_el_MS['std'].iloc[:,3:] / (Dict_el_MS['dat'].iloc[:,3:]
#                                                          * np.log(10) ) ,
#         Dict_el_MS['std'].iloc[:,3:] / (Dict_el_MS['dat'].iloc[:,3:]
#                                                          * np.log(10) ),
#         Dict_el_MS['std'].iloc[:,3:] / (Dict_el_MS['dat'].iloc[:,3:]
#                                                         * np.log(10) )],
#         y_list = [ Dict_ads['rsq'].iloc[:,1:N_sa],
#                   Dict_ads['rsq'].iloc[:,N_sa+1: 2*N_sa],
#                   Dict_ads['rsq'].iloc[:,2*N_sa+1:]],
#         std_y_list =[ Dict_ads['rsq_std'].iloc[:,1:N_sa],
#         Dict_ads['rsq_std'].iloc[:,N_sa+1:2*N_sa],
#         Dict_ads['rsq_std'].iloc[:,2*N_sa+1:] ],
#     element_index = Dict_ads['Qe'].index,
#     x_label = 'log($C_0$ [M])', y_label ="$(C_0-C_e)/C_0 [\%]$",
#     labels = ['Repl 1', 'Repl 2', 'Repl 3'],
#     folder_name = 'Rsq_repl', pre_title_plt = "Relative Sorbed Amount ", 
#     pre_save_name = 'rsq', Logscale = 0, Nucl_rel= Elem_rel, plot_everything= 1)
''' Analysis

Nice plots bro!!! I have a lot of plot, lets trz to pseudo analye them

*Relevants
    .Cs: rsq decrease as C increase! from 90 to 84
    .Eu: released always except for last values!
    .Mo: constant, no sorption nor release
    .Sr: release, typical patern as before
    .U: S pattern, first increase, then drecrease, for high conc==< competition?
    .
    BEntonite elements
    .Released, but the release decrease as conc increase!
    
*Other folder
    .Am: sorption high, 95%, decrease as conc increase!
    .Np: sorbed like 50%
    .Pu: high sorption, decrease from 95 to 50% as conc increase!! (last samples)
    .S: release, but release decrease as conc increase
    .Sm: from released to sorbed for high conc!!
    .
'''

######Iso Qe repl
# Read_JRC.ICPMS_Plotter_mean_blk_N (
#         x_list = [np.log10(Dict_el_sa['dat'].iloc[:,1:N_sa] ),
#                   np.log10(Dict_el_sa['dat'].iloc[:,N_sa+1: 2*N_sa] ),
#                   np.log10(Dict_el_sa['dat'].iloc[:,2*N_sa+1:] )  ], 
#         std_x_list = [ Dict_el_sa['std'].iloc[:,1:N_sa] / (Dict_el_sa['dat'].iloc[:,1:N_sa]
#                                                          * np.log(10) ),
#         Dict_el_sa['std'].iloc[:,N_sa+1: 2*N_sa] / (Dict_el_sa['dat'].iloc[:,N_sa+1: 2*N_sa]
#                                                          * np.log(10) ),
#         Dict_el_sa['std'].iloc[:,2*N_sa+1:] / (Dict_el_sa['dat'].iloc[:,2*N_sa+1:]
#                                                         * np.log(10) )],
#         y_list = [ np.log10(Dict_ads['Qe'].iloc[:,1:N_sa]),
#                   np.log10(Dict_ads['Qe'].iloc[:,N_sa+1: 2*N_sa]),
#                   np.log10(Dict_ads['Qe'].iloc[:,2*N_sa+1:]) ],
#         std_y_list =[ Dict_ads['Qe_std'].iloc[:,1:N_sa] / (np.log(10) * 
#                             np.abs(Dict_ads['Qe'].iloc[:,1:N_sa]) ),
#         Dict_ads['Qe_std'].iloc[:,N_sa+1:2*N_sa] /  (np.log(10) * 
#                             np.abs(Dict_ads['Qe'].iloc[:,N_sa+1:2*N_sa]) ),
#         Dict_ads['Qe_std'].iloc[:,2*N_sa+1:] /  (np.log(10) * 
#                          np.abs(Dict_ads['Qe'].iloc[:,2*N_sa+1:] ) )  ],
#     element_index = Dict_ads['Qe'].index,
#     x_label = 'log($C_e$ [M])', y_label ="log$(Q_e[mol/kg_{be}])$",
#     labels = ['Repl 1', 'Repl 2', 'Repl 3'],
#     folder_name = 'Iso_Qe_repl', pre_title_plt = "Adsorption isotherm ", 
#     pre_save_name = 'Qe_loglin', Logscale = 0, Nucl_rel= Elem_rel, plot_everything= 1)
''' Analysis

Interesting plots, and chaotic bro, wow. Nevertheless, I see some elements display
similar behaviour to the others
'''






################################################################
#%% ############# 5) Mean calcs ##############
#################################################################

'''
Well, here the replicates do not differ that much actually, so I will do the
 mean plot witohut removing any!!!

'''


Dict_Qe = Read_JRC.ICPMS_MeanStd_calculator(Dict_ads['Qe'], Dict_ads['Qe_std'],
                                            Nrepl = 3)
Dict_Kd = Read_JRC.ICPMS_MeanStd_calculator(Dict_ads['Kd'], Dict_ads['Kd_std'],
                                            Nrepl = 3)
# df_ppb_sa_S_Df2_mean, df_ppb_sa_S_Df2_std = Read_JRC.ICPMS_MeanStd_calculator(
#     df_ppb_sa_S_Df2.drop(['S Ad U 1.1', 'S Ad U 2.1','S Ad U 3.1'], axis = 1), Nrepl = 3)
#Dict_C0__Ceq = Read_JRC.ICPMS_MeanStd_calculator( Dict_ads['C0-Ceq'], Nrepl = 3)
Dict_rsq = Read_JRC.ICPMS_MeanStd_calculator( Dict_ads['rsq'], Dict_ads['rsq_std'], 
                                             Nrepl = 3)
Dict_el_sa_avg = Read_JRC.ICPMS_MeanStd_calculator( Dict_bent_rem['dat_br'], 
                                          Dict_bent_rem['std_br'],Nrepl = 3)
            #M data
Dict_Cleach = Read_JRC.ICPMS_MeanStd_calculator( Dict_bent_rem['C_leach'], 
                                        Dict_bent_rem['std_C_leach'],Nrepl = 3)
Dict_C0__Ceq = Read_JRC.ICPMS_MeanStd_calculator(Dict_ads['C0-Ceq'], 
                                        Dict_ads['C0-Ceq_std'],Nrepl = 3)


#Isotop version
#Beware, for Sr90, where sometimes I have NaN, for the mean NaN are not taken into account,
#but for the std, they are taken into account, NaN is the std!
Dict_Qe_iso = Read_JRC.ICPMS_MeanStd_calculator(Dict_ads_iso['Qe'], Dict_ads_iso['Qe_std'],
                                            Nrepl = 3)
Dict_Kd_iso = Read_JRC.ICPMS_MeanStd_calculator(Dict_ads_iso['Kd'], Dict_ads_iso['Kd_std'],
                                            Nrepl = 3)
Dict_rsq_iso = Read_JRC.ICPMS_MeanStd_calculator( Dict_ads_iso['rsq'], Dict_ads_iso['rsq_std'], 
                                             Nrepl = 3)
Dict_el_sa_avg_iso = Read_JRC.ICPMS_MeanStd_calculator( Dict_M['dat'], 
                                          Dict_M['dat_std'],Nrepl = 3)
#rstd of Sr90 3000% xD, huge in rsd plots


Dict_Asq_avg = Read_JRC.ICPMS_MeanStd_calculator( Asq['Asq'], 
                                          Asq['Asq_std'],Nrepl = 3)




#-------------------
#%% -------- 8) Correlatoin analysis ----------
#---------------
'''
As Paul Carbol told me (great his Phd, ahve a look, also to his researchgate), 
also gpt, and Even Sonia, lets do some correlations tests!!

I could do with 
            .rsq
            .Qe
            .Kd

I will firstly try with rsq, and then I will do with the others if needed.
Since I saw that for rsq, bentonite elements rsq increase, and U,Pun,Pu, La
decrease, I will firstly explore those.

Bradbury2005 modelling do correlation simply by plotting 1variable agains other, 
and fitting it!!

'''


Rsq_corr= Read_JRC.ICPMS_Correlation_test(Dict_rsq['< >'], 
                                      Elem_rel, type = 'Both')

'''

    *Pea analysis (lineal)
        .Hihg r (>0.8) between all bentonite elements!
        .High between U and La, Mo, Pu, Np, Am, Cs, Sn, Nd
    .No clear relations between bentonite elements and relevant elements!
    
    *Spe
        .Hihgh r between all bentonite elemenets!
        .High between U and Cs, La, Pu, Np, Am (0.6 the min there)
        .NEgative relations U and bent elem, from -0.4 to -0.6
        .Similar for Pu, Am, 
        .For Np not so high r!
'''

# Qe_corr= Read_JRC.ICPMS_Correlation_test(Dict_Qe['< >'], 
#                                       Elem_rel, type = 'Both')
'''
Does that really say something? 

#-----Pea (lin)

*Focusing at U, oit is strongly related to all the elements that are sorbed,
    Pu, Np, Am, Cm, Cs, La, ...
    though, hichest correlation for Pu, Np, 0.98

#--- SPe (non lin)
    *U
    .Corr 1 for Cs, Pu, Np, Am, Sn, Nd, 0.99 for La, 0
    0.5 for bent el 0.76 Zr
    
'''


#-------- Rsq U vs rsq Si
plt.figure(figsize=(12,10))  #width, heigh 6.4*4.8 inches by default
#plt.title("rsq of U (x) vs rsq of Si (y)", fontsize=22, wrap=True) #title
plt.errorbar(x = Dict_rsq['< >'].loc['U(LR)'].iloc[1:], 
             xerr =Dict_rsq['std'].loc['U(LR)'].iloc[1:],
              y =  Dict_rsq['< >'].loc['Si(MR)'].iloc[1:],
              yerr = Dict_rsq['std'].loc['Si(MR)'].iloc[1:], #removing blank!
              fmt = 'o', markersize = Markersize )
plt.ylabel("$(C_{eq} - C_0)/C_0$ Si [%]", fontsize= Font)              #ylabel
plt.xlabel("$(C_{eq} - C_0)/C_0$ U [%]", fontsize= Font)   
#plt.xticks(X_axis, Dict_el_MS['dat'].columns, rotation=90)
#plt.yscale('log') 
plt.legend(fontsize = Font)
plt.tick_params(axis='both', labelsize=Font)              #size of axis
plt.minorticks_on()             #enabling minor grid lines
plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        
                                #which both to plot major and minor grid lines
plt.grid(which = 'major')
plt.savefig('rsq_U_Si.png', format='png', bbox_inches='tight')
plt.show()  
'''
We see some kind of curve! Some some correlation exist, but not linear. Hence,
Pearson did not work!

That correlation indicates that, for high rsq, we have low Si, and viceversa!

Not C0 of U changed, but of Si no!
'''

#-------- Rsq U vs rsq Pu
plt.figure(figsize=(12,10))  #width, heigh 6.4*4.8 inches by default
#plt.title("rsq of U (x) vs rsq of Pu (y)", fontsize=22, wrap=True) #title
plt.errorbar(x = Dict_rsq['< >'].loc['U(LR)'].iloc[1:], 
             xerr =Dict_rsq['std'].loc['U(LR)'].iloc[1:],
              y =  Dict_rsq['< >'].loc['Pu(LR)'].iloc[1:],
              yerr = Dict_rsq['std'].loc['Pu(LR)'].iloc[1:], #removing blank!
              fmt = 'o', markersize = Markersize ) 
plt.ylabel("$(C_{eq} - C_0)/C_0$ Pu [%]", fontsize= Font)              #ylabel
plt.xlabel("$(C_{eq} - C_0)/C_0$ U [%]", fontsize= Font)   
#plt.xticks(X_axis, Dict_el_MS['dat'].columns, rotation=90)
#plt.yscale('log') 
plt.legend(fontsize = Font)
plt.tick_params(axis='both', labelsize=Font)              #size of axis
plt.minorticks_on()             #enabling minor grid lines
plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        
                                #which both to plot major and minor grid lines
plt.grid(which = 'major')
plt.savefig('rsq_U_Pu.png', format='png', bbox_inches='tight')
plt.show()  
'Here we see a linear pattern, they are related!'




#-----------------------------------------------------------
#%% 				 10) Plots mean values 
#-----------------------------------------------------------
'''
Okay, I already did the plot all
.Q
.Kd
.rsq

For all the samples. 

*Multielemental plots:
    That is a nice idea. I could so for the CL eleemnts to compare with the
    S Ad CL
    .I could also gather between oxidation status and so:
        .U(VI), NP(IV), Pu(IV, V)
        .Cs (I)
        .LA (III), CM(III), Am(III), Eu(III)
		
		
The initial samples will be removed from the plots since those points are not so nice,
and create a tail-like pattern
'''

S_rem = ['Data 1','Data 2', 'Data 3'] 	# Samples to remove
MS_rem = ['S Ad Leach 0.1', 'S Ad Leach 0.2', 'S Ad Leach 0.3'] #MS samples to remove
Removed_name = '123'    #Samples numbers used above

print('Bro, did u ensure the samples removed for plotting are okay, and names also?')
print('The samples removed are:')
print(S_rem)
print(MS_rem)
print('And the name: ' + Removed_name + '\n')
print('Does it all make sense? if not, modify these things!')


#---------- Bent elements conc ------

Read_JRC.ICPMS_Barplotter (Dict_el_sa_avg['< >'], Dict_el_sa_avg['%rsd'], 
        folder_name = 'Bar_plot_conc', Nucl_rel= Elem_rel_ben + ['Na(MR)'], 
        ylabel_1= 'Conc [M]')

'''
This I wanted to save, to have the final conc of those elements. MAybe could indicate
sth?
'''

#---------------- Cleach ------------------------------------
#Lets plot this, might be interesting/reveal something?

# plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
# plt.title("Cleach, for the SAdLeach exp", fontsize=22, wrap=True) #title

# for el in ['Si(MR)', 'Al(MR)', 'Mg(MR)', 'Mn(MR)', 'Fe(MR)', 
#             'P(MR)', 'Ti(MR)','S(MR)','Sr(LR)']:       
#             #plotting each element of the Bent_elem array except Ca(MR)
#     plt.errorbar(x = Dict_el_MS['dat'].loc[el][1:], xerr =Dict_el_MS['std'].loc[el][1:],
#              y = -Dict_C0__Ceq['< >'].loc[el][1:],
#              yerr = Dict_C0__Ceq['std'].loc[el][1:], #removing blank!
#              fmt = 'o:', label = '<' + el + '>') 
# plt.ylabel("$C_{eq} - C_0$ [M]", fontsize=18)              #ylabel
# plt.xlabel('$C_0$ [M]', fontsize = Font)
# plt.tick_params(axis='both', labelsize=18)              #size of axis
# plt.yscale('log')
# plt.minorticks_on()
# plt.grid(which='minor', linewidth=0.5)
# plt.grid(which='major') 
# plt.legend(fontsize = Font)
# plt.savefig('Cleach_SAdLeach.png', format='png', bbox_inches='tight')
# plt.show()
'''
That is a bad plot, since the x axis is similar for almost all the samples,
but the y value vary. I do not know how could I improve it, other than
doing a bar plot
'''


# ----------- Cleach vs [U]_0, element by elements  .loc['U(LR)']
Read_JRC.ICPMS_Plotter_mean_blk_N(
    x_list = [ Dict_el_MS['dat'].loc['U(LR)'][3:] ],
    std_x_list = [ Dict_el_MS['std'].loc['U(LR)'][3:] ] ,
    y_list = [ -Dict_C0__Ceq['< >'].drop(S_rem, axis = 1)  ],
    std_y_list= [ Dict_C0__Ceq['std'].drop(S_rem, axis = 1) ],
    element_index = Dict_el_MS['dat'].index, labels= ['<S>'],
    x_label = '$[U]_0 [M]$', 
    y_label = "$ C_{eq}-C_0$ [M]", folder_name = 'Cleach_SAdLeach',
    pre_title_plt = "Cleach vs initial U concentration ", pre_save_name = 'Cleach_plot',
    Nucl_rel=Elem_rel_ben, Blank_here= 0,
    colors= [ Bent_color['Sar']])
'''
Really strange plots:
    .For Sr, Si, many elements not clear pattern. The inital release and later a
     sorption even suggest that precipitation could ahve occured!! 
    .For Mg, it gets to a Plateau

'''

#------- rsq vs [U]_0
Read_JRC.ICPMS_Plotter_mean_blk_N(
    x_list = [ Dict_el_MS['dat'].loc['U(LR)'][3:] ],
    std_x_list = [ Dict_el_MS['std'].loc['U(LR)'][3:] ] ,
    y_list = [ Dict_rsq['< >'].drop(S_rem, axis = 1)  ],
    std_y_list= [ Dict_rsq['std'].drop(S_rem, axis = 1) ],
    element_index = Dict_el_MS['dat'].index, labels= ['<S>'],
    x_label = '$[U]_0 [M]$', 
    y_label = "$ (C_0- C_{eq})/C_0$ [%]", folder_name = 'Cleach_rsq_SAdLeach',
    pre_title_plt = "rsq vs initial U concentration ", pre_save_name = 'Cleach_plot',
    Nucl_rel=Elem_rel_ben, Blank_here= 0,
    colors= [ Bent_color['Sar']])


# Multibar plot for Cleach
Read_JRC.ICPMS_MultiBar_plotter( -Dict_C0__Ceq['< >'],  Dict_C0__Ceq['std'], 
            ['Si(MR)', 'Sr(LR)','Al(MR)','Mg(MR)','Fe(MR)','S(MR)'], b = 0.1,
        Xlabel = '< Samples >', Ylabel = "$C_{eq} - C_0$ [M]", 
        Title = 'Cleach for the S Ad Cleach exp',
        Savename = 'Cleach_SAdCleach', Log_y= False)

#Si only
Read_JRC.ICPMS_MultiBar_plotter( -Dict_C0__Ceq['< >'],  Dict_C0__Ceq['std'], 
            ['Al(MR)'], b = 0.1,
        Xlabel = '< Samples >', Ylabel = "$C_{eq} - C_0$ [M]", 
        Title = 'Cleach for the S Ad Cleach exp',
        Savename = 'Cleach_Al_SAdCleach', Log_y= False)
        #Single element also strange for Si, Sr, Al
        #For Mg decreasing as N increase
'''
Note that due to the scale:
    .No Sr data can be seen! (it is freleased in lower amount!)
    
Well, a strange plot, since:
    .For 8,9,10, Ceq-Co for Si <0!
    .From 8 to 11, Ceq-Co for Al <0!!
    .Mg seems stable
    .Na not included since huge errorbars!!    
Na and Ca chaotic (remember difficult to measure!)

To get a better sense, lets do the same but for the rsq also!
    
'''



# -------- Cleach as function of C_e [U for ex]
'''
A interesting plot could be to plot Ceq-C0 of bentonite elements as a function
of initial concentration of sorbing elements (u, Cs, etc), not bentonite elements.
That might be more clear than the bar plots. LEts give it a try!
'''

# plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
# plt.title("Cleach as function of [U]", fontsize=22, wrap=True) #title

# Elem_rel_ben = ['Sr(LR)', 'Si(MR)']
# for el in Elem_rel_ben:        #plotting each element of the Bent_elem array
#     plt.errorbar(x = Dict_el_MS['dat'].loc['U(LR)'][1:], 
#         xerr = Dict_el_MS['std'].loc['U(LR)'][1:],
#         y = -Dict_C0__Ceq['< >'].loc[el][1:], yerr = Dict_C0__Ceq['std'].loc[el][1:], 
#              fmt = 'o:', label = '<' + el + '>') 
# plt.ylabel("$C_{eq} - C_0$ [M]", fontsize= Font)              #ylabel
# plt.xlabel('$[U]_0$ [M]', fontsize = Font)
# plt.tick_params(axis='both', labelsize= Font)              #size of axis
# #plt.yscale('log') 
# plt.xscale('log') 
# plt.minorticks_on()             #enabling minor grid lines
# plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        #which both to plot major and minor grid lines
# plt.grid(which = 'major') 
# plt.legend(fontsize = Font)
# plt.savefig('Cleach_plot.png', format='png', bbox_inches='tight')
# plt.show()
'''
Well, I think that this is a really interesting plot! Would be good to do it for all
the elements separate, and then merge them in case of needed. Really nice plot, this
gives you a bettter sense of Cleach!
'''





#---------- Rsq ----------

# Multibar plot for rsq !
Read_JRC.ICPMS_MultiBar_plotter(Dict_rsq['< >'], Dict_rsq['std'], 
            ['Si(MR)', 'Sr(LR)','Al(MR)','Mg(MR)','Fe(MR)','S(MR)'], b = 0.1,
        Xlabel = '< Samples >', Ylabel = 'rsq [%]', 
        Title = 'Cleach (rsq version) for the S Ad Cleach exp',
        Savename = 'Cleach_rsq_SAdCleach')
'''
Here we see clearly the odd pattern for the bentonite elements, really released
for the initial samples, but then it decrease. What is this bro?
????

And a relevant question: is this new, or we have it for the cold exp as well??
??
'''


#Elem
Read_JRC.ICPMS_Plotter_mean_blk_N (
        x_list = [np.log10(Dict_el_MS['dat'].iloc[:,3:]) ], 
        std_x_list = [Dict_el_MS['std'].iloc[:,3:] / 
                      (Dict_el_MS['dat'].iloc[:,3:]*np.log(10)) ],
        y_list = [Dict_rsq['< >'].iloc[:,3:] ],
        std_y_list =[Dict_rsq['std'].iloc[:,3:] ],
    element_index = Dict_Qe['std'].index,
                x_label = 'log($C_0$ [M])', y_label ="$(C_0-C_e)/C_0 [\%]$",
                labels = ['Sar'], 
                colors= [Bent_color['Sar']],
                folder_name = 'Rsq', pre_title_plt = "rel_sorbed_am SAdLeach ", 
                pre_save_name = 'rsq', Nucl_rel= Elem_rel, plot_everything= 0)
'''
Note:
    .Tc not sorbed! In iso Qe plot more evident. In literature only 1 study
    revealing sorption can be found, Grambow2006, he had reducing conditions,
    so Tc is Tc IV. Under ox, Tc is Tc VII [Chapman1986]
'''

### Isotop
Read_JRC.ICPMS_Plotter_mean_blk_N (
        x_list = [np.log10(Dict_M_MS['dat'].iloc[:,3:]) ], 
        std_x_list = [Dict_M_MS['dat_std'].iloc[:,3:] / 
                      (Dict_M_MS['dat'].iloc[:,3:]*np.log(10)) ],
        y_list = [Dict_rsq_iso['< >'].iloc[:,3:] ],
        std_y_list =[Dict_rsq_iso['std'].iloc[:,3:] ],
    element_index = Dict_Qe_iso['std'].index,
                x_label = 'log($C_0$ [M])', y_label ="$(C_0-C_e)/C_0 [\%]$",
                labels = ['Sar'], 
                colors= [Bent_color['Sar']],
                folder_name = 'Rsq_iso', pre_title_plt = "rel_sorbed_am SAdLeach ", 
                pre_save_name = 'rsq', Nucl_rel= Iso_rel, plot_everything= 0)
'''
Well, Sr90 tricky. Negative value for Data 2, but for the other
still I have huge errorbars
'''


#Multielemental rsq (element)
'''
Okay, which plots could I do:
    i) Actinides: U, Pu, Np, Am, Cm
    ii) CL: U, Cs, La, Eu
    iii) Lanth +  analogues: Am, Cm, La, Nd, Pr, Sm
'''
#       Multielem rsq Act
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title("rsq for the S Ad leach", 
#          fontsize=22, wrap=True) #title     
            #plotting each element of the Bent_elem array except Ca(MR)
plt.errorbar(np.log10(Dict_el_MS['dat'].drop(
    MS_rem,   axis = 1).loc['U(LR)'] ), 
    xerr = (Dict_el_MS['std']/ (Dict_el_MS['dat']*np.log(10))).drop(
                  MS_rem,  axis = 1).loc['U(LR)'],
    y = Dict_rsq['< >'].drop(S_rem, axis = 1).loc['U(LR)'],
    yerr = Dict_rsq['std'].drop(S_rem, axis = 1).loc['U(LR)'], 
   fmt = 'o', label = 'U', markersize = Markersize, color = "#4682B4",
   elinewidth = Linewidth)    #U
plt.errorbar(np.log10(Dict_el_MS['dat'].drop(
    MS_rem,   axis = 1).loc['Np(LR)'] ), 
    xerr = (Dict_el_MS['std']/ (Dict_el_MS['dat']*np.log(10))).drop(
                  MS_rem,  axis = 1).loc['Np(LR)'],
    y = Dict_rsq['< >'].drop(S_rem, axis = 1).loc['Np(LR)'],
    yerr = Dict_rsq['std'].drop(S_rem, axis = 1).loc['Np(LR)'], 
   fmt = 'o', label = 'Np', markersize = Markersize, elinewidth = Linewidth,
   color = "#228B22")      #Np
plt.errorbar(np.log10(Dict_el_MS['dat'].drop(
    MS_rem,   axis = 1).loc['Pu(LR)'] ), 
    xerr = (Dict_el_MS['std']/ (Dict_el_MS['dat']*np.log(10))).drop(
                  MS_rem,  axis = 1).loc['Pu(LR)'],
    y = Dict_rsq['< >'].drop(S_rem, axis = 1).loc['Pu(LR)'],
    yerr = Dict_rsq['std'].drop(S_rem, axis = 1).loc['Pu(LR)'], 
   fmt = 'o', label = 'Pu', markersize = Markersize, elinewidth = Linewidth,
   color = "#DC143C")  #Pu
plt.errorbar(np.log10(Dict_el_MS['dat'].drop(
    MS_rem,   axis = 1).loc['Am(LR)'] ), 
    xerr = (Dict_el_MS['std']/ (Dict_el_MS['dat']*np.log(10))).drop(
                  MS_rem,  axis = 1).loc['Am(LR)'],
    y = Dict_rsq['< >'].drop(S_rem, axis = 1).loc['Am(LR)'],
    yerr = Dict_rsq['std'].drop(S_rem, axis = 1).loc['Am(LR)'], 
   fmt = 'o', label = 'Am', markersize = Markersize, elinewidth = Linewidth,
   color = '#FF8C00')      #Am
plt.errorbar(np.log10(Dict_el_MS['dat'].drop(
    MS_rem,   axis = 1).loc['Cm(LR)'] ), 
    xerr = (Dict_el_MS['std']/ (Dict_el_MS['dat']*np.log(10))).drop(
                  MS_rem,  axis = 1).loc['Cm(LR)'],
    y = Dict_rsq['< >'].drop(S_rem, axis = 1).loc['Cm(LR)'],
    yerr = Dict_rsq['std'].drop(S_rem, axis = 1).loc['Cm(LR)'], 
   fmt = 'o', label = 'Cm', markersize = Markersize, elinewidth = Linewidth,
   color = "#9932CC")      #Cm
plt.ylabel("Removal [%]", fontsize = Font)              #ylabel
plt.xlabel('log($C_0$ [M])', fontsize = Font)
plt.tick_params(axis='both', labelsize = Font)              #size of axis
#plt.yscale('log')
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth=0.5)
plt.grid(which='major')
plt.legend(fontsize = Font)
plt.savefig('Rsq_act_No' + Removed_name +'.png', format='png', bbox_inches='tight')
plt.show()

#       Multielem rsq lanth Cs
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title("rsq for the S Ad leach", 
#          fontsize=22, wrap=True) #title     
            #plotting each element of the Bent_elem array excRept Ca(MR)
plt.errorbar(np.log10(Dict_el_MS['dat'].drop(
    MS_rem,   axis = 1).loc['Cs(LR)'] ), 
    xerr = (Dict_el_MS['std']/ (Dict_el_MS['dat']*np.log(10))).drop(
                  MS_rem,  axis = 1).loc['Cs(LR)'],
    y = Dict_rsq['< >'].drop(S_rem, axis = 1).loc['Cs(LR)'],
    yerr = Dict_rsq['std'].drop(S_rem, axis = 1).loc['Cs(LR)'], 
   fmt = 'o', label = 'Cs', markersize = Markersize, elinewidth = Linewidth,
   color = "#708090")  #Cs
plt.errorbar(np.log10(Dict_el_MS['dat'].drop(
    MS_rem,   axis = 1).loc['La(LR)'] ), 
    xerr = (Dict_el_MS['std']/ (Dict_el_MS['dat']*np.log(10))).drop(
                  MS_rem,  axis = 1).loc['La(LR)'],
    y = Dict_rsq['< >'].drop(S_rem, axis = 1).loc['La(LR)'],
    yerr = Dict_rsq['std'].drop(S_rem, axis = 1).loc['La(LR)'], 
   fmt = 'o', label = 'La', markersize = Markersize, elinewidth = Linewidth,
   color = "#DAA520")       #La
plt.errorbar(np.log10(Dict_el_MS['dat'].drop(
    MS_rem,   axis = 1).loc['Nd(LR)'] ), 
    xerr = (Dict_el_MS['std']/ (Dict_el_MS['dat']*np.log(10))).drop(
                  MS_rem,  axis = 1).loc['Nd(LR)'],
    y = Dict_rsq['< >'].drop(S_rem, axis = 1).loc['Nd(LR)'],
    yerr = Dict_rsq['std'].drop(S_rem, axis = 1).loc['Nd(LR)'], 
   fmt = 'o', label = 'Nd', markersize = Markersize, elinewidth = Linewidth,
   color = "#FF69B4")      #Nd
plt.ylabel("Removal [%]", fontsize = Font)              #ylabel
plt.xlabel('log($C_0$ [M])', fontsize = Font)
plt.tick_params(axis='both', labelsize = Font)              #size of axis
#plt.yscale('log')
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth=0.5)
plt.grid(which='major')
plt.legend(fontsize = Font)
plt.savefig('Rsq_lanthCs_No' + Removed_name +'.png', format='png', bbox_inches='tight')
plt.show()

#       Multielem rsq CL
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title("rsq for the S Ad leach", 
#          fontsize=22, wrap=True) #title     
            #plotting each element of the Bent_elem array except Ca(MR)
plt.errorbar(np.log10(Dict_el_MS['dat'].drop(
    MS_rem,   axis = 1).loc['U(LR)'] ), 
    xerr = (Dict_el_MS['std']/ (Dict_el_MS['dat']*np.log(10))).drop(
                  MS_rem,  axis = 1).loc['U(LR)'],
    y = Dict_rsq['< >'].drop(S_rem, axis = 1).loc['U(LR)'],
    yerr = Dict_rsq['std'].drop(S_rem, axis = 1).loc['U(LR)'], 
   fmt = 'o', label = 'U', markersize = Markersize, color = "#4682B4",
   elinewidth = Linewidth)    #U
plt.errorbar(np.log10(Dict_el_MS['dat'].drop(
    MS_rem,   axis = 1).loc['Cs(LR)'] ), 
    xerr = (Dict_el_MS['std']/ (Dict_el_MS['dat']*np.log(10))).drop(
                  MS_rem,  axis = 1).loc['Cs(LR)'],
    y = Dict_rsq['< >'].drop(S_rem, axis = 1).loc['Cs(LR)'],
    yerr = Dict_rsq['std'].drop(S_rem, axis = 1).loc['Cs(LR)'], 
   fmt = 'o', label = 'Cs', markersize = Markersize, elinewidth = Linewidth,
   color = "#708090")  #Cs
plt.errorbar(np.log10(Dict_el_MS['dat'].drop(
    MS_rem,   axis = 1).loc['La(LR)'] ), 
    xerr = (Dict_el_MS['std']/ (Dict_el_MS['dat']*np.log(10))).drop(
                  MS_rem,  axis = 1).loc['La(LR)'],
    y = Dict_rsq['< >'].drop(S_rem, axis = 1).loc['La(LR)'],
    yerr = Dict_rsq['std'].drop(S_rem, axis = 1).loc['La(LR)'], 
   fmt = 'o', label = 'La', markersize = Markersize, elinewidth = Linewidth,
   color = "#DAA520")       #La
plt.ylabel("$(C_0 - C_e) / C_0$ [%]", fontsize = Font)              #ylabel
plt.xlabel('log($C_0$ [M])', fontsize = Font)
plt.tick_params(axis='both', labelsize = Font)              #size of axis
#plt.yscale('log')
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth=0.5)
plt.grid(which='major')
plt.legend(fontsize = Font)
plt.savefig('Rsq_CL_No' + Removed_name +'.png', format='png', bbox_inches='tight')
plt.show()


#-------------------------------------------------------------------
#----------------------------------- Ads iso, Q ----------------------------


#lin
#Read_JRC.ICPMS_Plotter_mean_blk_N (
 #   x_list = [np.log10(Dict_el_sa_avg['< >'].iloc[:,3:])],
  #  std_x_list= [np.abs(Dict_el_sa_avg['std'].iloc[:,3:]/(Dict_el_sa_avg['< >'].iloc[:,3:]*np.log(10)))],
   # y_list = [Dict_Qe['< >'].iloc[:,3:] ],
    #std_y_list = [Dict_Qe['std'].iloc[:,3:] ],
    #element_index = Dict_Qe['std'].index,
    #x_label = 'log($C_e$ [M])', y_label ="$Q_e$ [mol/kg$_{be}$]",
    #labels = ['S'],
    #folder_name = 'Isoth_Qe', pre_title_plt = "Isoth Qe SAdLeach ", 
    #pre_save_name = 'IsoQe_lin', Nucl_rel=Elem_rel, 
    #colors= [Bent_color['Sar']] , plot_everything= 0)
'''
    Not clear pattern between the released elements, Sr, Si, Al, etc
'''

#loglin
Read_JRC.ICPMS_Plotter_mean_blk_N (
    x_list = [np.log10(Dict_el_sa_avg['< >'].iloc[:,3:])],
    std_x_list= [np.abs(Dict_el_sa_avg['std'].iloc[:,3:]/(Dict_el_sa_avg['< >'].iloc[:,3:]*np.log(10)))],
    y_list = [np.log10(Dict_Qe['< >'].iloc[:,3:]) ],
    std_y_list = [np.abs(Dict_Qe['std'].iloc[:,3:]/(Dict_Qe['< >'].iloc[:,3:]*np.log(10))) ],
    element_index = Dict_Qe['std'].index,
    x_label = 'log($C_e$ [M])', y_label ="log($Q_e$ [mol/kg$_{be}$])",
    labels = ['S'],
    folder_name = 'Isoth_Qe', pre_title_plt = "Isoth Qe SAdLeach ", 
    pre_save_name = 'IsoQe_loglin', Nucl_rel=Elem_rel, 
    colors= [Bent_color['Sar']] , plot_everything= 0)

''' ########## Analysis ##############

-Madness bro xD

'''


#---- Loglin iso
Read_JRC.ICPMS_Plotter_mean_blk_N (
    x_list = [np.log10(Dict_el_sa_avg_iso['< >'].iloc[:,3:])],
    std_x_list= [np.abs(Dict_el_sa_avg_iso['std'].iloc[:,3:]/
                        (Dict_el_sa_avg_iso['< >'].iloc[:,3:]*np.log(10)))],
    y_list = [np.log10(Dict_Qe_iso['< >'].iloc[:,3:]) ],
    std_y_list = [np.abs(Dict_Qe_iso['std'].iloc[:,3:]/(Dict_Qe_iso['< >'].iloc[:,3:]*np.log(10))) ],
    element_index = Dict_Qe_iso['std'].index,
    x_label = 'log($C_e$ [M])', y_label ="log($Q_e$ [mol/kg$_{be}$])",
    labels = ['S'],
    folder_name = 'Isoth_Qe_iso', pre_title_plt = "Isoth Qe SAdLeach ", 
    pre_save_name = 'IsoQe_loglin', Nucl_rel=Iso_rel, 
    colors= [Bent_color['Sar']] , plot_everything= 0)


#Sr90plot
#Data 6 of Qe iso Sr90 is NaN!
Skip=8
Read_JRC.ICPMS_Plotter_mean_blk_N (
    x_list = [np.log10(Dict_el_sa_avg_iso['< >'].iloc[:,Skip:])],
    std_x_list= [np.abs(Dict_el_sa_avg_iso['std'].iloc[:,Skip:]/
                        (Dict_el_sa_avg_iso['< >'].iloc[:,Skip:]*np.log(10)))],
    y_list = [np.log10(Dict_Qe_iso['< >'].iloc[:,Skip:]) ],
    std_y_list = [np.abs(Dict_Qe_iso['std'].iloc[:,Skip:]/(Dict_Qe_iso['< >'].iloc[:,Skip:]*np.log(10))) ],
    element_index = Dict_Qe_iso['std'].index,
    x_label = 'log($C_e$ [M])', y_label ="log($Q_e$ [mol/kg$_{be}$])",
    labels = ['S'],
    folder_name = 'Isoth_Qe_iso', pre_title_plt = "Isoth Qe SAdLeach ", 
    pre_save_name = 'IsoQe_loglin_Sr_'+str(Skip)+'_on', Nucl_rel=['Sr90(LR)'], 
    colors= [Bent_color['Sar']] , plot_everything= 0)

'''
Okay, not so clear the Sr90, have a look, maybe because of the high rsd for some samples?
!!
!!!

Data 4,5 with huge Ce error, wtf bro! maybe just indicate not trustworthy data!!

Maybe just is not so trustworthy bro, relax!!
!!!!
!!

Both isotopes of Eu, 151, 153 showing the strange pattern, more Ba like, indicating 
that they are Ba intereferencs. Maybe Eu154 (fission) could reveal that?
'''

#--------- Multielement plot Qe, UCsLaEu (no Sr!)
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title("$Q_e$ adsorption isotherm for the S Ad leach", 
#          fontsize=22, wrap=True) #title     
            #plotting each element of the Bent_elem array except Ca(MR)
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem,   axis = 1).loc['U(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['U(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem, axis = 1).loc['U(LR)']* np.log(10) ),
    y = np.log10(Dict_Qe['< >'].drop(S_rem,  axis = 1).loc['U(LR)']),
    yerr = np.abs( (Dict_Qe['std']/ (Dict_Qe['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['U(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'U', markersize = Markersize)      #U
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem,   axis = 1).loc['Cs(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['Cs(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['Cs(LR)']* np.log(10) ),
    y = np.log10(Dict_Qe['< >'].drop(S_rem,  axis = 1).loc['Cs(LR)']),
    yerr = np.abs( (Dict_Qe['std']/ (Dict_Qe['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['Cs(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'Cs', markersize = Markersize)             #Cs
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem,   axis = 1).loc['La(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['La(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem, axis = 1).loc['La(LR)']* np.log(10) ),
    y = np.log10(Dict_Qe['< >'].drop(S_rem,  axis = 1).loc['La(LR)']),
    yerr = np.abs( (Dict_Qe['std']/ (Dict_Qe['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['La(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'La', markersize = Markersize)             #La
plt.ylabel("log($Q_e$ [mol/kg$_{be}$])", fontsize = Font)              #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font)
plt.tick_params(axis='both', labelsize = Font)              #size of axis
#plt.yscale('log')
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth=0.5)
plt.grid(which='major')
plt.legend(fontsize = Font)
plt.savefig('Iso_Qe_loglin_CL_No' + Removed_name +'.png', format='png', bbox_inches='tight')
plt.show()



#------ Multielem Qe actini
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title("$Q_e$ adsorption isotherm for the S Ad leach", 
#          fontsize=22, wrap=True) #title     
            #plotting each element of the Bent_elem array except Ca(MR)
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['U(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['U(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['U(LR)']* np.log(10) ),
    y = np.log10(Dict_Qe['< >'].drop(S_rem,  axis = 1).loc['U(LR)']),
    yerr = np.abs( (Dict_Qe['std']/ (Dict_Qe['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['U(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'U', markersize = Markersize)      #U
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['Np(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['Np(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['Np(LR)']* np.log(10) ),
    y = np.log10(Dict_Qe['< >'].drop(S_rem,  axis = 1).loc['Np(LR)']),
    yerr = np.abs( (Dict_Qe['std']/ (Dict_Qe['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['Np(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'Np', markersize = Markersize)             #Np
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem, axis = 1).loc['Pu(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['Pu(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem, axis = 1
                ).loc['Pu(LR)']* np.log(10) ),
    y = np.log10(Dict_Qe['< >'].drop(S_rem,  axis = 1).loc['Pu(LR)']),
    yerr = np.abs( (Dict_Qe['std']/ (Dict_Qe['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['Pu(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'Pu', markersize = Markersize)             #Pu
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem, axis = 1).loc['Am(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['Am(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem, axis = 1
                ).loc['Am(LR)']* np.log(10) ),
    y = np.log10(Dict_Qe['< >'].drop(S_rem,  axis = 1).loc['Am(LR)']),
    yerr = np.abs( (Dict_Qe['std']/ (Dict_Qe['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['Am(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'Am', markersize = Markersize)        #Am
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem, axis = 1).loc['Cm(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['Cm(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem, axis = 1
                ).loc['Cm(LR)']* np.log(10) ),
    y = np.log10(Dict_Qe['< >'].drop(S_rem,  axis = 1).loc['Cm(LR)']),
    yerr = np.abs( (Dict_Qe['std']/ (Dict_Qe['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['Cm(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'Cm', markersize = Markersize)         #Cm
plt.ylabel("log($Q_e$ [mol/kg$_{be}$])", fontsize = Font)              #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font)
plt.tick_params(axis='both', labelsize = Font)              #size of axis
#plt.yscale('log')
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth=0.5)
plt.grid(which='major')
plt.legend(fontsize = Font)
plt.savefig('IsoQe_loglin_Actini_No' + Removed_name +'.png', format='png', bbox_inches='tight')
plt.show()


#--------- Multielement plot Qe lantanides
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title("$Q_e$ adsorption isotherm for the S Ad leach", 
#          fontsize=22, wrap=True) #title     
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem, axis = 1).loc['Cs(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['Cs(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem, axis = 1
                ).loc['Cs(LR)']* np.log(10) ),
    y = np.log10(Dict_Qe['< >'].drop(S_rem,  axis = 1).loc['Cs(LR)']),
    yerr = np.abs( (Dict_Qe['std']/ (Dict_Qe['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['Cs(LR)']  ), elinewidth = Linewidth, 
              fmt = 'o', label = 'Cs', markersize = Markersize)     #Cs            
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['La(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['La(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['La(LR)']* np.log(10) ),
    y = np.log10(Dict_Qe['< >'].drop(S_rem,  axis = 1).loc['La(LR)']),
    yerr = np.abs( (Dict_Qe['std']/ (Dict_Qe['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['La(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'La', markersize = Markersize)      #La
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem, axis = 1).loc['Nd(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['Nd(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem, axis = 1
                ).loc['Nd(LR)']* np.log(10) ),
    y = np.log10(Dict_Qe['< >'].drop(S_rem,  axis = 1).loc['Nd(LR)']),
    yerr = np.abs( (Dict_Qe['std']/ (Dict_Qe['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['Nd(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'Nd', markersize = Markersize)     #Nd
plt.ylabel("log($Q_e$ [mol/kg$_{be}$])", fontsize = Font)              #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font)
plt.tick_params(axis='both', labelsize = Font)              #size of axis
#plt.yscale('log')
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth=0.5)
plt.grid(which='major')
plt.legend(fontsize = Font)
plt.savefig('IsoQe_loglin_lanthCs_No' + Removed_name +'.png', format='png', bbox_inches='tight')
plt.show()


#------ Qe, plot La vs La+ all actinides ---------------
#Lanthanides = La, Nd, Er, Dy, Ce, Pr, Eu, Sm, Gd, Tb, Ho, Yb, Lu
#I select the ones with better graphs

Qe_Lant = Dict_Qe['< >'].loc[['La(LR)','Ce(LR)', 'Nd(LR)', 'Pr(LR)', 'Sm(LR)'
                              ,'Am(LR)','Cm(LR)']].sum(axis = 0)
Qe_Lant_std = np.sqrt(Dict_Qe['std'].loc['La(LR)']**2 + Dict_Qe['std'].loc['Ce(LR)']**2 +
             Dict_Qe['std'].loc['Nd(LR)']**2 + Dict_Qe['std'].loc['Pr(LR)']**2  +
             Dict_Qe['std'].loc['Sm(LR)']**2 +  Dict_Qe['std'].loc['Am(LR)']**2
             +  Dict_Qe['std'].loc['Cm(LR)']**2 ) 

plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title("$Q_e$ adsorption isotherm for the S Ad leach", 
#          fontsize=22, wrap=True) #title     
            #plotting each element of the Bent_elem array except Ca(MR)
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['La(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['La(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['La(LR)']* np.log(10) ),
    y = np.log10(Dict_Qe['< >'].drop(S_rem,  axis = 1).loc['La(LR)']),
    yerr = np.abs(Dict_Qe['std'].drop(S_rem, axis = 1).loc['La(LR)'] /
        (Dict_Qe['< >'].drop(S_rem, axis = 1 ).loc['La(LR)'] * np.log(10) ) ), 
              fmt = 'o', label = 'La', markersize = Markersize)      #La
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem, axis = 1).loc['La(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['La(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem, axis = 1
                ).loc['La(LR)']* np.log(10) ),
    y = np.log10(Qe_Lant.drop(S_rem)),
    yerr = np.abs(Qe_Lant_std.drop(S_rem)  /(Qe_Lant.drop(S_rem) * np.log(10) ) ) , 
              fmt = 'o', label = 'La+Ce+Nd+Pr+Sm+Am+Cm', markersize = Markersize)        #Sum of lantha
plt.ylabel("log($Q_e$[mol/kg$_{be}$])", fontsize = Font)              #ylabel
plt.xlabel('log$([La]_e$ [M])', fontsize = Font)
plt.tick_params(axis='both', labelsize = Font)              #size of axis
#plt.yscale('log')
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth=0.5)
plt.grid(which='major')
plt.legend(fontsize = Font)
plt.savefig('IsoQe_loglin_lanth_sum_No' + Removed_name +'.png', format='png', bbox_inches='tight')
plt.show()
'''
Ojo, the x data si the Ce conc of La always!
I wanted to do this to see if I can confirm my hypothesis; La sorption is lower
in hot sorption relative to the cold one because of many more lanthanides present 
in hot conditions, which compete foir sorption!

Well, kinda, still lower, but higher the sum, closer to the cold one!
'''


#--------- Single plot Qe Sr90 
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title("$Q_e$ adsorption isotherm for the S Ad leach", 
#          fontsize=22, wrap=True) #title     
plt.errorbar(
  np.log10(Dict_el_sa_avg_iso['< >'].drop(S_rem, axis = 1
                ).loc['Sr90(LR)'] ).iloc[4:], 
    xerr = (Dict_el_sa_avg_iso['std'].drop(S_rem, 
        axis = 1).loc['Sr90(LR)']/ (Dict_el_sa_avg_iso['< >'].drop(
        S_rem, axis = 1
                ).loc['Sr90(LR)']* np.log(10) )).iloc[4:],
    y = np.log10(Dict_Qe_iso['< >'].drop(S_rem,  axis = 1
                                         ).loc['Sr90(LR)']).iloc[4:],
    yerr = np.abs(Dict_Qe_iso['std'].drop(S_rem, axis = 1
                                          ).loc['Sr90(LR)'] /
        (Dict_Qe_iso['< >'].drop(S_rem, axis = 1 
                                 ).loc['Sr90(LR)'] * np.log(10) ) ).iloc[4:], 
              fmt = 'o', label = 'Am', markersize = Markersize)        #Sr90
plt.ylabel("log($Q_e$ [mol/kg$_{be}$])", fontsize = Font)              #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font)
plt.tick_params(axis='both', labelsize = Font)              #size of axis
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth=0.5)
plt.grid(which='major')
#plt.legend(fontsize = Font)
plt.savefig('IsoQe_iso_Sr90_567891011.png', format='png', bbox_inches='tight')
plt.show()
'''
Okay, note Data 1,2,3 always removed. THe other depends:
    .Plotting data 91011: good, 3 points only [5:]
    .Plotting 891011: 4 points, the data 8 already huge errorbars! [4:]
    .Plotting 7891011: data 7 is already diffferent, no linear trend. ANd 
    huuge errorbars [3:]
    .Data 6 is NaN (not clear)
    .Including Data 5 [1:] truly fuckes up the plot. Hence we know, until 8 only!
    
Plotting all the data creates a bad plot due to the errorbars in Ce. Removing
first 4 is better, but still not perfect

4: and 5: do the same, since a NaN is there!
'''


#-----------------------------------------------------------------
#                            Ads iso, Kd 
#-----------------------------------------------------------------
Read_JRC.ICPMS_Plotter_mean_blk_N (
    x_list = [np.log10(Dict_el_sa_avg['< >'].iloc[:,3:])],
    std_x_list= [np.abs(Dict_el_sa_avg['std'].iloc[:,3:]/(Dict_el_sa_avg['< >'].iloc[:,3:]*np.log(10)))],
    y_list = [np.log10(Dict_Kd['< >'].iloc[:,3:]) ],
    std_y_list = [np.abs(Dict_Kd['std'].iloc[:,3:]/(Dict_Kd['< >'].iloc[:,3:]*np.log(10))) ],
    element_index = Dict_Kd['std'].index,
    x_label = 'log($C_e$ [M])', y_label ="log($K_d$ [L/kg$_{be}$])",
    labels = ['S'],
    folder_name = 'Isoth_Kd', pre_title_plt = "Isoth Kd SAdLeach ", 
    pre_save_name = 'IsoKd', Nucl_rel=Elem_rel, 
    colors= [Bent_color['Sar']], plot_everything= 0 )
'''
Key plots!!

*La,Am, Eu
.Patter inrease and then decrease (concave) as conc increases!

*Pu, Np, U
.Same idea, incerase and then decrease, but this time with a decrease in the 
middle, making an M like shape! wtf?

'''

#---- Loglin iso
Read_JRC.ICPMS_Plotter_mean_blk_N (
    x_list = [np.log10(Dict_el_sa_avg_iso['< >'].iloc[:,3:])],
    std_x_list= [np.abs(Dict_el_sa_avg_iso['std'].iloc[:,3:]/
                        (Dict_el_sa_avg_iso['< >'].iloc[:,3:]*np.log(10)))],
    y_list = [np.log10(Dict_Kd_iso['< >'].iloc[:,3:]) ],
    std_y_list = [np.abs(Dict_Kd_iso['std'].iloc[:,3:]/(Dict_Kd_iso['< >'].iloc[:,3:]*np.log(10))) ],
    element_index = Dict_Kd_iso['std'].index,
    x_label = 'log($C_e$ [M])', y_label ="log($K_d$ [L/kg$_{be}$])",
    labels = ['S'],
    folder_name = 'Isoth_Kd_iso', pre_title_plt = "Isoth Kd SAdLeach ", 
    pre_save_name = 'IsoKd', Nucl_rel=Iso_rel, 
    colors= [Bent_color['Sar']] , plot_everything= 0)



#--------- Single plot Kd Sr90 
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title("$Q_e$ adsorption isotherm for the S Ad leach", 
#          fontsize=22, wrap=True) #title     
plt.errorbar(
  np.log10(Dict_el_sa_avg_iso['< >'].drop(S_rem, axis = 1
                                          ).loc['Sr90(LR)'] ).iloc[5:], 
    xerr = (Dict_el_sa_avg_iso['std'].drop(S_rem, 
        axis = 1).loc['Sr90(LR)']/ (Dict_el_sa_avg_iso['< >'].drop(
            S_rem, axis = 1
                ).loc['Sr90(LR)']* np.log(10) )).iloc[5:],
    y = np.log10(Dict_Kd_iso['< >'].drop(S_rem,  axis = 1
                                         ).loc['Sr90(LR)']).iloc[5:],
    yerr = np.abs(Dict_Kd_iso['std'].drop(S_rem, axis = 1
                                          ).loc['Sr90(LR)'] /
        (Dict_Kd_iso['< >'].drop(S_rem, axis = 1 
                                 ).loc['Sr90(LR)'] * np.log(10) ) ).iloc[5:], 
              fmt = 'o', label = 'Am', markersize = Markersize)        #Sr90
plt.ylabel("log($K_d$ [L/kg$_{be}$])", fontsize = Font)              #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font)
plt.tick_params(axis='both', labelsize = Font)              #size of axis
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth=0.5)
plt.grid(which='major')
#plt.legend(fontsize = Font)
plt.savefig('IsoKd_iso_Sr90_91011.png', format='png', bbox_inches='tight')
plt.show()


#--------- Multielement plot Kd, UCsLaEu (no Sr!)
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title("$K_d$ adsorption isotherm for the S Ad leach", 
#          fontsize=22, wrap=True) #title     
            #plotting each element of the Bent_elem array except Ca(MR)
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem,   axis = 1).loc['U(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['U(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem,    axis = 1).loc['U(LR)']* np.log(10) ),
    y = np.log10(Dict_Kd['< >'].drop(S_rem,  axis = 1).loc['U(LR)']),
    yerr = np.abs( (Dict_Kd['std']/ (Dict_Kd['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['U(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'U', markersize = Markersize)      #U
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem,   axis = 1).loc['Cs(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['Cs(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem,    axis = 1).loc['Cs(LR)']* np.log(10) ),
    y = np.log10(Dict_Kd['< >'].drop(S_rem,  axis = 1).loc['Cs(LR)']),
    yerr = np.abs( (Dict_Kd['std']/ (Dict_Kd['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['Cs(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'Cs', markersize = Markersize)             #Cs
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem,   axis = 1).loc['La(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['La(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem,    axis = 1).loc['La(LR)']* np.log(10) ),
    y = np.log10(Dict_Kd['< >'].drop(S_rem,  axis = 1).loc['La(LR)']),
    yerr = np.abs( (Dict_Kd['std']/ (Dict_Kd['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['La(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'La', markersize = Markersize)             #La
plt.ylabel("log($K_d$ [L/kg$_{be}$])", fontsize = Font)              #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font)
plt.tick_params(axis='both', labelsize = Font)              #size of axis
#plt.yscale('log')
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth=0.5)
plt.grid(which='major')
plt.legend(fontsize = Font)
plt.savefig('IsoKd_loglin_CL_No' + Removed_name +'.png', format='png', bbox_inches='tight')
plt.show()

#------ Multielem Kd actini
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title("$K_d$ adsorption isotherm for the S Ad leach", 
#          fontsize=22, wrap=True) #title     
            #plotting each element of the Bent_elem array except Ca(MR)
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['U(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['U(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['U(LR)']* np.log(10) ),
    y = np.log10(Dict_Kd['< >'].drop(S_rem,  axis = 1).loc['U(LR)']),
    yerr = np.abs( (Dict_Kd['std']/ (Dict_Kd['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['U(LR)']  ), elinewidth = Linewidth,
    fmt = 'o', label = 'U', markersize = Markersize, color = '#4682B4')      #U
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['Np(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['Np(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['Np(LR)']* np.log(10) ),
    y = np.log10(Dict_Kd['< >'].drop(S_rem,  axis = 1).loc['Np(LR)']),
    yerr = np.abs( (Dict_Kd['std']/ (Dict_Kd['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['Np(LR)']  ), elinewidth = Linewidth,
    fmt = 'o', label = 'Np', markersize = Markersize, color = '#228B22') #Np
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem, axis = 1).loc['Pu(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['Pu(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem, axis = 1
                ).loc['Pu(LR)']* np.log(10) ),
    y = np.log10(Dict_Kd['< >'].drop(S_rem,  axis = 1).loc['Pu(LR)']),
    yerr = np.abs( (Dict_Kd['std']/ (Dict_Kd['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['Pu(LR)']  ), elinewidth = Linewidth,
    fmt = 'o', label = 'Pu', markersize = Markersize, color = '#DC143C')   #Pu
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem, axis = 1).loc['Am(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['Am(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem, axis = 1
                ).loc['Am(LR)']* np.log(10) ),
    y = np.log10(Dict_Kd['< >'].drop(S_rem,  axis = 1).loc['Am(LR)']),
    yerr = np.abs( (Dict_Kd['std']/ (Dict_Kd['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['Am(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'Am', markersize = Markersize, color = '#FF8C00')    
                                                    #Am
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem, axis = 1).loc['Cm(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['Cm(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem, axis = 1
                ).loc['Cm(LR)']* np.log(10) ),
    y = np.log10(Dict_Kd['< >'].drop(S_rem,  axis = 1).loc['Cm(LR)']),
    yerr = np.abs( (Dict_Kd['std']/ (Dict_Kd['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['Cm(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'Cm', markersize = Markersize, color = '#9932CC')  
                                                                    #Cm
plt.ylabel("log($K_d$ [L/kg$_{be}$])", fontsize = Font)  
plt.xlabel('log($C_e$ [M])', fontsize = Font)
plt.tick_params(axis='both', labelsize = Font)              #size of axis
#plt.yscale('log')
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth=0.5)
plt.grid(which='major')
plt.legend(fontsize = Font)
plt.savefig('IsoKd_loglin_Actini_No' + Removed_name +'.png', format='png', bbox_inches='tight')
plt.show()



#--------- Multielement plot Kd lantanides Cs
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title("$K_d$ adsorption isotherm for the S Ad leach", 
#          fontsize=22, wrap=True) #title 
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem, axis = 1).loc['Cs(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['Cs(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem, axis = 1
                ).loc['Cs(LR)']* np.log(10) ),
    y = np.log10(Dict_Kd['< >'].drop(S_rem,  axis = 1).loc['Cs(LR)']),
    yerr = np.abs( (Dict_Kd['std']/ (Dict_Kd['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['Cs(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'Cs', markersize = Markersize, color = '#708090')  #Cs    
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['La(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['La(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem,axis = 1).loc['La(LR)']* np.log(10) ),
    y = np.log10(Dict_Kd['< >'].drop(S_rem,  axis = 1).loc['La(LR)']),
    yerr = np.abs( (Dict_Kd['std']/ (Dict_Kd['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['La(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'La', markersize = Markersize, color = '#DAA520')   #La
plt.errorbar(np.log10(Dict_el_sa_avg['< >'].drop(S_rem, axis = 1).loc['Nd(LR)'] ), 
    xerr = Dict_el_sa_avg['std'].drop(S_rem, 
        axis = 1).loc['Nd(LR)']/ (Dict_el_sa_avg['< >'].drop(S_rem, axis = 1
                ).loc['Nd(LR)']* np.log(10) ),
    y = np.log10(Dict_Kd['< >'].drop(S_rem,  axis = 1).loc['Nd(LR)']),
    yerr = np.abs( (Dict_Kd['std']/ (Dict_Kd['< >'] * np.log(10) )).drop(
        S_rem, axis = 1 ).loc['Nd(LR)']  ), elinewidth = Linewidth,
              fmt = 'o', label = 'Nd', markersize = Markersize, color = '#FF69B4')  #Nd
plt.ylabel("log($K_d$ [L/kg$_{be}$])", fontsize = Font)             #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font)
plt.tick_params(axis='both', labelsize = Font)              #size of axis
#plt.yscale('log')
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth=0.5)
plt.grid(which='major')
plt.legend(fontsize = Font)
plt.savefig('IsoKd_loglin_lanthCs_No' + Removed_name +'.png', format='png', bbox_inches='tight')
plt.show()


#--------------------------------------------------------------
#%% --------------- 9) Comparison with cold experiments -----------------
#---------------------------------------------------------------------
'''
The ebst way to compare between experiments is to plot them together! 

In first instance I will only compare with the CL exp! MAybe in the future I
include all the exp!

I will remove some samples from CL, as needed

For mean plot in hot, I need to remove sample 1, since for Cs, Qe>0!! new here,
normally I was obtaining Qe<0, so no issue in including in loglin iso!
'''

SAdU = joblib.load('Dat_SAdU_20251205.pkl') 
# SAdUCs = joblib.load('Dat_SAdUCs_20251008.pkl')         #S Ad UCs basif2 data!
# SAdUCsLa = joblib.load('Dat_SAdUCsLa_20251008.pkl') 
SAdCL = joblib.load('Dat_SAdCL_20251023.pkl') 


#CL_basif1 = ['Data 10', 'Data 11']
#CL_acid =['Data 7', 'Data 8', 'Data 9']
#CL_basif2 = ['Data 12', 'Data 13', 'Data 14']

U_basif = ['Data 11']
U_acid =['Data 10']

#--------------------------------------------------------------------------
#-------------- Comparison with S Ad CL ----------------------------------
#-------------------------------------------------------------------------

#The samples removed from the Ad CL are:

CL_S_rem = ['Data 1', 'Data 2', 'Data 3', 'Data 7', 'Data 8', 'Data 9', 
         'Data 10', 'Data 11']
CL_MS_rem = ['Ad CL 0.1','Ad CL 0.2' , 'Ad CL 0.3', 'Ad CL 0.7', 'Ad CL 0.8', 
          'Ad CL 0.9', 'Ad CL 0.10', 'Ad CL 0.11'] 

####### rsq plot #######
Read_JRC.ICPMS_Plotter_mean_blk_N (
        x_list = [ np.log10(Dict_el_MS['dat'].iloc[:,3:]),
                  np.log10(SAdCL['Dict_el_MS']['dat'].drop(CL_MS_rem, axis = 1)) ], 
        std_x_list = [Dict_el_MS['std'].iloc[:,3:] / 
                      (Dict_el_MS['dat'].iloc[:,3:]*np.log(10)),
                      (SAdCL['Dict_el_MS']['std']/ 
                      (SAdCL['Dict_el_MS']['dat']*np.log(10))).drop(CL_MS_rem, axis = 1)],
        y_list = [Dict_rsq['< >'].iloc[:,3:],
                  SAdCL['<rsq[%]>']['< >'].drop(CL_S_rem, axis = 1) ],
        std_y_list =[ Dict_rsq['std'].iloc[:,3:], 
                     SAdCL['<rsq[%]>']['std'].drop(CL_S_rem, axis = 1) ],
    element_index = Dict_Qe['std'].index,
                x_label = 'log($C_0$ [M])', y_label ="$(C_0-C_e)/C_0 [\%]$",
                labels = ['Hot', 'Cold'], 
                Nucl_rel= ['Cs(LR)', 'Eu(LR)','La(LR)', 'U(LR)' ]+ Elem_rel_ben,
                folder_name = 'Hot_vs_cold/Rsq', pre_title_plt = "rel_sorbed_am ", 
                pre_save_name = 'rsq_No' + Removed_name)

#Read_JRC.ICPMS_Plotter_mean_blk_N (
 #       x_list = [ np.log10(Dict_el_MS['dat'].iloc[:,3:]),
  #                np.log10(SAdCL['Dict_el_MS']['dat'].drop(['Ad CL 0.1'], axis = 1)) ], 
   #     std_x_list = [Dict_el_MS['std'].iloc[:,3:] / 
    #                  (Dict_el_MS['dat'].iloc[:,3:]*np.log(10)),
     #                 SAdCL['Dict_el_MS']['std'].drop(['Ad CL 0.1'], axis = 1) / 
      #                (SAdCL['Dict_el_MS']['dat'].drop(['Ad CL 0.1'], axis = 1)*np.log(10))],
     #   y_list = [Dict_rsq['< >'].iloc[:,3:],
       #           SAdCL['<rsq[%]>']['< >'].drop(['Data 1'], axis = 1) ],
      #  std_y_list =[ Dict_rsq['std'].iloc[:,3:], 
        #             SAdCL['<rsq[%]>']['std'].drop(['Data 1'], axis = 1) ],
#    element_index = Dict_Qe['std'].index,
 #               x_label = 'log($C_0$ [M])', y_label ="$(C_0-C_e)/C_0 [\%]$",
  #              labels = ['Hot', 'Cold'], 
   #             Nucl_rel= ['Cs(LR)', 'Eu(LR)','La(LR)', 'U(LR)' ],
    #            folder_name = 'Hot_vs_cold/rsq', pre_title_plt = "rel_sorbed_am ", 
     #           pre_save_name = 'rsq_no1')


########### Iso Qe ############
Read_JRC.ICPMS_Plotter_mean_blk_N (
    x_list = [np.log10(Dict_el_sa_avg['< >']).drop(S_rem, axis = 1),
              np.log10(SAdCL['<Dict_M>']['< >'].drop(CL_S_rem, axis = 1) ) ],
    std_x_list= [np.abs(Dict_el_sa_avg['std']/(Dict_el_sa_avg['< >']*np.log(10))).drop(
        S_rem, axis = 1),
    np.abs(SAdCL['<Dict_M>']['std'].drop(CL_S_rem, axis = 1)/
           (SAdCL['<Dict_M>']['< >'].drop(CL_S_rem, axis = 1) *np.log(10))) ],
    y_list = [np.log10(Dict_Qe['< >']).drop(S_rem, axis = 1),
              np.log10(SAdCL['<Qe[mol/kg]>']['< >'].drop(CL_S_rem, axis = 1) ) ],
    std_y_list = [np.abs(Dict_Qe['std']/
                         (Dict_Qe['< >']*np.log(10))).drop(S_rem, axis = 1),
        np.abs(SAdCL['<Qe[mol/kg]>']['std'].drop(CL_S_rem, axis = 1)/
              (SAdCL['<Qe[mol/kg]>']['< >'].drop(CL_S_rem, axis = 1)*np.log(10)))  ],
    element_index = Dict_Qe['std'].index,
    x_label = 'log($C_e$ [M])', y_label ="log($Q_e$ [mol/kg$_{be}$])",
    labels = ['Hot', 'Cold'], Nucl_rel= ['Cs(LR)', 'Eu(LR)','La(LR)', 'U(LR)' ],
    folder_name = 'Hot_vs_cold/IsoQe', pre_title_plt = "Ads iso ", 
    pre_save_name = 'IsoQe_No' + Removed_name)


#################Iso Kd ###########
Read_JRC.ICPMS_Plotter_mean_blk_N (
    x_list = [np.log10(Dict_el_sa_avg['< >']).drop(S_rem, axis = 1),
              np.log10(SAdCL['<Dict_M>']['< >'].drop(CL_S_rem, axis = 1) ) ],
    std_x_list= [np.abs(Dict_el_sa_avg['std']/(Dict_el_sa_avg['< >']*np.log(10))).drop(
        S_rem, axis = 1),
    np.abs(SAdCL['<Dict_M>']['std'].drop(CL_S_rem, axis = 1) /
           (SAdCL['<Dict_M>']['< >'].drop(CL_S_rem, axis = 1)*np.log(10)))],
    y_list = [np.log10(Dict_Kd['< >'].drop(S_rem, axis = 1)),
              np.log10(SAdCL['<Kd [L/kg]>']['< >'].drop(CL_S_rem, axis = 1) )],
    std_y_list = [np.abs(Dict_Kd['std']/ (Dict_Kd['< >']*np.log(10))).drop(S_rem, axis = 1),
                  np.abs(SAdCL['<Kd [L/kg]>']['std'].drop(CL_S_rem, axis = 1)/
                        (SAdCL['<Kd [L/kg]>']['< >'].drop(CL_S_rem, axis = 1)*np.log(10))) ],
    element_index = Dict_Kd['std'].index,
    x_label = 'log($C_e$ [M])', y_label ="log($K_d$ [L/kg$_{be}$])",
    labels = ['Hot', 'Cold'], Nucl_rel= ['Cs(LR)', 'Eu(LR)','La(LR)', 'U(LR)' ],
    folder_name = 'Hot_vs_cold/IsoKd', pre_title_plt = "Ads iso ", 
    pre_save_name = 'IsoKd_No' + Removed_name)




#------------ Bentonite leached elements comparison -----------
'''
his could be really revealing, are there any difference in the elached bentonite
elements? These could be:
        Si, Al, Mg, Fe, Mn, Ti, S, Sr, etc

We cans tart by checking Si alone, the most abundant
'''


#--------------- Si leached ------------------

# plt.errorbar(x = Dict_el_MS['dat'].loc[el][1:], xerr =Dict_el_MS['std'].loc[el][1:],)
# #              y = -Dict_C0__Ceq['< >'].loc[el][1:],
# #              yerr = Dict_C0__Ceq['std'].loc[el][1:], #removing blank!
# #              fmt = 'o:', label = '<' + el + '>') 

# Read_JRC.ICPMS_MultiBar_plotter_N(
#     [-Dict_C0__Ceq['< >'].drop(S_rem,axis = 1), 
#                 -SAdCL['<C_0-C_eq[M]>']['< >'].drop(['Data 1','Data 2','Data 3', 'Data 12', 'Data 13',
#                                                      'Data 14'],axis = 1)],
#     [Dict_C0__Ceq['std'].drop(S_rem,axis = 1), 
#                 SAdCL['<C_0-C_eq[M]>']['std'].drop(['Data 1','Data 2','Data 3', 'Data 12', 'Data 13',
#                                                      'Data 14'],axis = 1)],
#             ['Si(MR)'], b = 0.1,
#             dataset_labels = ['Hot', 'Cold'],
#         Xlabel = '< Samples >', Ylabel = 'C_eq - C_0 [M]', 
#         Title = 'Si leach for hot vs cold exp',
#         Savename = 'Si_leach_hot_cold')


#Trial to see if adding Qe of lanthanides gives Qe of La in CL exp!
np.log10(Dict_Qe['< >'].loc[['La(LR)', 'Eu(LR)', 'Nd(LR)', 'Er(LR)', 'Dy(LR)', 
                    'Ce(LR)','Pr(LR)','Sm(LR)','Gd(LR)','Tb(LR)','Dy(LR)',
                    'Er(LR)','Yb(LR)','Lu(LR)']].sum(axis = 0) )

#Well, kinda, -4.3 obtained, -5 was only La in hot, in cold -4!
# !!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!


#--------------------------------------------------------------------------
#-------------- Comparison with S Ad U ----------------------------------
#-------------------------------------------------------------------------
'''
Since the U conc for the S ad U had lwer values, I will also compare, provided
that here I explored lower values 
'''


####### rsq plot #######
Read_JRC.ICPMS_Plotter_mean_blk_N (
        x_list = [ np.log10(Dict_el_MS['dat'].iloc[:,3:]),
                  np.log10(SAdU['Dict_el_MS']['dat'].drop(['S Ad U 0.1',
                'S Ad U 0.10','S Ad U 0.11'], axis = 1)) ], 
        std_x_list = [Dict_el_MS['std'].iloc[:,3:] / 
                      (Dict_el_MS['dat'].iloc[:,3:]*np.log(10)),
                      SAdU['Dict_el_MS']['std'].drop(['S Ad U 0.1',
                                    'S Ad U 0.10','S Ad U 0.11'], axis = 1) / 
                      (SAdU['Dict_el_MS']['dat'].drop(['S Ad U 0.1',
                         'S Ad U 0.10','S Ad U 0.11'], axis = 1)*np.log(10))],
        y_list = [Dict_rsq['< >'].iloc[:,3:],
                  SAdU['<rsq[%]>']['< >'].drop(['Data 1','Data 10','Data 11'], axis = 1) ],
        std_y_list =[ Dict_rsq['std'].iloc[:,3:], 
                     SAdU['<rsq[%]>']['std'].drop(['Data 1','Data 10','Data 11'], axis = 1) ],
    element_index = Dict_Qe['std'].index,
                x_label = 'log($C_0$ [M])', y_label ="$(C_0-C_e)/C_0 [\%]$",
                labels = ['Hot', 'S Ad U'], 
                Nucl_rel= [ 'U(LR)' ]+ Elem_rel_ben,
                folder_name = 'Hot_vs_AdU/rsq', pre_title_plt = "rel_sorbed_am ", 
                pre_save_name = 'rsq_noAcidBasif')

Read_JRC.ICPMS_Plotter_mean_blk_N (
        x_list = [ np.log10(Dict_el_MS['dat'].iloc[:,3:]),
                  np.log10(SAdU['Dict_el_MS']['dat'].drop(['S Ad U 0.1'], axis = 1)) ], 
        std_x_list = [Dict_el_MS['std'].iloc[:,3:] / 
                      (Dict_el_MS['dat'].iloc[:,3:]*np.log(10)),
                      SAdU['Dict_el_MS']['std'].drop(['S Ad U 0.1'], axis = 1) / 
                      (SAdU['Dict_el_MS']['dat'].drop(['S Ad U 0.1'], axis = 1)*np.log(10))],
        y_list = [Dict_rsq['< >'].iloc[:,3:],
                  SAdU['<rsq[%]>']['< >'].drop(['Data 1'], axis = 1) ],
        std_y_list =[ Dict_rsq['std'].iloc[:,3:], 
                     SAdU['<rsq[%]>']['std'].drop(['Data 1'], axis = 1) ],
    element_index = Dict_Qe['std'].index,
                x_label = 'log($C_0$ [M])', y_label ="$(C_0-C_e)/C_0 [\%]$",
                labels = ['Hot', 'S Ad U'], 
                Nucl_rel= [ 'U(LR)' ],
                folder_name = 'Hot_vs_AdU/rsq', pre_title_plt = "rel_sorbed_am ", 
                pre_save_name = 'rsq_no1')


########### Iso Qe ############
Read_JRC.ICPMS_Plotter_mean_blk_N (
    x_list = [np.log10(Dict_el_sa_avg['< >']).drop(S_rem, axis = 1),
              np.log10(SAdU['<Dict_M>']['< >'].drop(U_acid + U_basif, axis = 1) ) ],
    std_x_list= [np.abs(Dict_el_sa_avg['std']/(Dict_el_sa_avg['< >']*np.log(10))).drop(
        S_rem, axis = 1),
    np.abs(SAdU['<Dict_M>']['std'].drop(U_acid + U_basif, axis = 1)/
           (SAdU['<Dict_M>']['< >'].drop(U_acid + U_basif, axis = 1) *np.log(10))) ],
    y_list = [np.log10(Dict_Qe['< >']).drop(S_rem, axis = 1),
              np.log10(SAdU['<Qe[mol/kg]>']['< >'].drop(U_acid + U_basif, axis = 1) ) ],
    std_y_list = [np.abs(Dict_Qe['std']/
                         (Dict_Qe['< >']*np.log(10))).drop(S_rem, axis = 1),
        np.abs(SAdU['<Qe[mol/kg]>']['std'].drop(U_acid + U_basif, axis = 1)/
              (SAdU['<Qe[mol/kg]>']['< >'].drop(U_acid + U_basif, axis = 1)*np.log(10)))  ],
    element_index = Dict_Qe['std'].index,
    x_label = 'log($C_e$ [M])', y_label ="log($Q_e$ [mol/kg$_{be}$])",
    labels = ['Hot', 'S Ad U'], Nucl_rel= [ 'U(LR)' ],
    folder_name = 'Hot_vs_AdU/IsoQe', pre_title_plt = "Ads iso ", 
    pre_save_name = 'IsoQe_noAcidBasif' )


#################Iso Kd ###########
Read_JRC.ICPMS_Plotter_mean_blk_N (
    x_list = [np.log10(Dict_el_sa_avg['< >']).drop(S_rem, axis = 1),
              np.log10(SAdU['<Dict_M>']['< >'].drop(U_acid + U_basif, axis = 1) ) ],
    std_x_list= [np.abs(Dict_el_sa_avg['std']/(Dict_el_sa_avg['< >']*np.log(10))).drop(
        S_rem, axis = 1),
    np.abs(SAdU['<Dict_M>']['std'].drop(U_acid + U_basif, axis = 1) /
           (SAdU['<Dict_M>']['< >'].drop(U_acid + U_basif, axis = 1)*np.log(10)))],
    y_list = [np.log10(Dict_Kd['< >'].drop(S_rem, axis = 1)),
              np.log10(SAdU['<Kd[L/kg]>']['< >'].drop(U_acid + U_basif, axis = 1) )],
    std_y_list = [np.abs(Dict_Kd['std']/ (Dict_Kd['< >']*np.log(10))).drop(S_rem, axis = 1),
                  np.abs(SAdU['<Kd[L/kg]>']['std'].drop(U_acid + U_basif, axis = 1)/
                        (SAdU['<Kd[L/kg]>']['< >'].drop(U_acid + U_basif, axis = 1)*np.log(10))) ],
    element_index = Dict_Kd['std'].index,
    x_label = 'log($C_e$ [M])', y_label ="log($K_d$ [L/kg$_{be}$])",
    labels = ['Hot', 'S Ad U'], Nucl_rel= [ 'U(LR)' ],
    folder_name = 'Hot_vs_AdU/IsoKd', pre_title_plt = "Ads iso ", 
    pre_save_name = 'IsoKd_noAcidBasif' )



#------------------------------------------------------------------------------
#               Hot vs Ad U vs CL
#------------------------------------------------------------------------------
#Comparison for U only!


########### Iso Qe ############
Read_JRC.ICPMS_Plotter_mean_blk_N (
    x_list = [np.log10(Dict_el_sa_avg['< >']).drop(S_rem, axis = 1),
              np.log10(SAdCL['<Dict_M>']['< >'].drop(CL_S_rem, axis = 1) ),
              np.log10(SAdU['<Dict_M>']['< >'].drop(U_acid + U_basif, axis = 1) ) ],
    std_x_list= [np.abs(Dict_el_sa_avg['std']/(Dict_el_sa_avg['< >']*np.log(10))).drop(
        S_rem, axis = 1),
    np.abs(SAdCL['<Dict_M>']['std'].drop(CL_S_rem, axis = 1)/
           (SAdCL['<Dict_M>']['< >'].drop(CL_S_rem, axis = 1) *np.log(10))),
    np.abs(SAdU['<Dict_M>']['std'].drop(U_acid + U_basif, axis = 1) /
           (SAdU['<Dict_M>']['< >'].drop(U_acid + U_basif, axis = 1)*np.log(10)))],
    y_list = [np.log10(Dict_Qe['< >']).drop(S_rem, axis = 1),
              np.log10(SAdCL['<Qe[mol/kg]>']['< >'].drop(CL_S_rem, axis = 1) ),
              np.log10(SAdU['<Qe[mol/kg]>']['< >'].drop(U_acid + U_basif, axis = 1) )],
    std_y_list = [np.abs(Dict_Qe['std']/
                         (Dict_Qe['< >']*np.log(10))).drop(S_rem, axis = 1),
        np.abs(SAdCL['<Qe[mol/kg]>']['std'].drop(CL_S_rem, axis = 1)/
              (SAdCL['<Qe[mol/kg]>']['< >'].drop(CL_S_rem, axis = 1)*np.log(10))),
        np.abs(SAdU['<Qe[mol/kg]>']['std']/(SAdU['<Qe[mol/kg]>']['< >']*np.log(10))).drop(
            U_acid + U_basif, axis = 1)],
    element_index = Dict_Qe['std'].index,
    x_label = 'log($C_e$ [M])', y_label ="log($Q_e$ [mol/kg$_{be}$])",
    labels = ['Hot', 'S Ad CL', 'S Ad U'], Nucl_rel= [ 'U(LR)' ],
    folder_name = 'Hot_vs_cold_AdU', pre_title_plt = "Ads iso ", 
    pre_save_name = 'IsoQe_noAcidBasif1' )

########### Iso Kd ############
Read_JRC.ICPMS_Plotter_mean_blk_N (
    x_list = [np.log10(Dict_el_sa_avg['< >']).drop(S_rem, axis = 1),
              np.log10(SAdCL['<Dict_M>']['< >'].drop(CL_S_rem, axis = 1) ),
              np.log10(SAdU['<Dict_M>']['< >'].drop(U_acid + U_basif, axis = 1) ) ],
    std_x_list= [np.abs(Dict_el_sa_avg['std']/(Dict_el_sa_avg['< >']*np.log(10))).drop(
        S_rem, axis = 1),
    np.abs(SAdCL['<Dict_M>']['std'].drop(CL_S_rem, axis = 1)/
           (SAdCL['<Dict_M>']['< >'].drop(CL_S_rem, axis = 1) *np.log(10))),
    np.abs(SAdU['<Dict_M>']['std'].drop(U_acid + U_basif, axis = 1) /
           (SAdU['<Dict_M>']['< >'].drop(U_acid + U_basif, axis = 1)*np.log(10)))],
    y_list = [np.log10(Dict_Kd['< >']).drop(S_rem, axis = 1),
              np.log10(SAdCL['<Kd [L/kg]>']['< >'].drop(CL_S_rem, axis = 1) ),
              np.log10(SAdU['<Kd[L/kg]>']['< >'].drop(U_acid + U_basif, axis = 1) )],
    std_y_list = [np.abs(Dict_Kd['std']/
                         (Dict_Kd['< >']*np.log(10))).drop(S_rem, axis = 1),
        np.abs(SAdCL['<Kd [L/kg]>']['std'].drop(CL_S_rem, axis = 1)/
              (SAdCL['<Kd [L/kg]>']['< >'].drop(CL_S_rem, axis = 1)*np.log(10))),
        np.abs(SAdU['<Kd[L/kg]>']['std']/(SAdU['<Kd[L/kg]>']['< >']*np.log(10))).drop(
            U_acid + U_basif, axis = 1)],
    element_index = Dict_Qe['std'].index,
    x_label = 'log($C_e$ [M])', y_label ="log($K_d$ [L/kg$_{be}$])",
    labels = ['Hot', 'S Ad CL', 'S Ad U'], Nucl_rel= [ 'U(LR)' ],
    folder_name = 'Hot_vs_cold_AdU', pre_title_plt = "Ads iso ", 
    pre_save_name = 'IsoKd_noAcidBasif1' )


#-----------------------------------------------
#%% ------------ Fitting ---------------------------
#--------------------------------------------------
'''
Okay, complex data. Focusing on the CL eleemnts, U Cs La Eu Sr:
    .La could be this time Langmuir!
    .U, Cs could be Freundlich
    .Eu was not sorbed
    .Sr was not sorbed neither
    
Lets give it a try. Also we can include other relevants actinides present:
    Pu, Np, Am, Cm (Cm maybe not so sorbed)
    

#Brooo! for Freund, NL gives better parameters (in aprticular the constant) than
the other linear for U and Cs!!! I guess because of the lin issue, of Kd \propto 1/intercept, 
Delta \propto 1/inter^2 , mkaing Delta Kd huge!
'''


#--------------- Freundlich fit ------------------------------------
Fre_fitU = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['U(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['U(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['U(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['U(LR)'], folder_name = 'Fre_fit', 
        Color = Bent_color['Sar'], Title = 'Freundlich fit for U', 
        save_name = 'Freund_fit_U', x_label = 'log$(C_e[M])$', 
        y_label = 'log($Q_e$ [mol/kg$_{be}$])')
Fre_fitU_NL = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['U(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['U(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['U(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['U(LR)'], folder_name = 'Fre_fit', 
        Color = Bent_color['Sar'], Title = 'Freundlich fit for U', 
        save_name = 'Freund_fit_U_NL', x_label = 'log$(C_e[M])$', 
        y_label = 'log($Q_e$ [mol/kg$_{be}$])', Fit_type= 0)
        #If I remove Data 2, would be way better bro!! now r = 0.92!
        
Fre_fitCs = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Cs(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Cs(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Cs(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Cs(LR)'], folder_name = 'Fre_fit', 
        Color = Bent_color['Sar'], Title = 'Freundlich fit for Cs', 
        save_name = 'Freund_fit_Cs', x_label = 'log$(C_e[M])$', 
        y_label = 'log($Q_e$ [mol/kg$_{be}$])')
                #r = 0.99! (rounded!)
Fre_fitCs_NL = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Cs(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Cs(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Cs(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Cs(LR)'], folder_name = 'Fre_fit', 
        Color = Bent_color['Sar'], Title = 'Freundlich fit for Cs', 
        save_name = 'Freund_fit_Cs_NL', x_label = 'log$(C_e[M])$', 
        y_label = 'log($Q_e$ [mol/kg$_{be}$])', Fit_type= 0)

Fre_fitLa = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['La(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['La(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['La(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['La(LR)'], folder_name = 'Fre_fit', 
        Color = Bent_color['Sar'], Title = 'Freundlich fit for La', 
        save_name = 'Freund_fit_La', x_label = 'log$(C_e[M])$', 
        y_label = 'log($Q_e$ [mol/kg$_{be}$])')
        # r = 0.91, not bad though!
Fre_fitLa_NL = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['La(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['La(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['La(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['La(LR)'], folder_name = 'Fre_fit', 
        Color = Bent_color['Sar'], Title = 'Freundlich fit for La', 
        save_name = 'Freund_fit_La_NL', x_label = 'log$(C_e[M])$', 
        y_label = 'log($Q_e$ [mol/kg$_{be}$])', Fit_type = 0)
Fre_fitNd = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Nd(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Nd(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Nd(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Nd(LR)'], folder_name = 'Fre_fit', 
        Color = Bent_color['Sar'], Title = 'Freundlich fit for Nd', 
        save_name = 'Freund_fit_Nd', x_label = 'log$(C_e[M])$', 
        y_label = 'log($Q_e$ [mol/kg$_{be}$])')
        # r = 0.91, not bad though!
Fre_fitNd_NL = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Nd(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Nd(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Nd(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Nd(LR)'], folder_name = 'Fre_fit', 
        Color = Bent_color['Sar'], Title = 'Freundlich fit for Nd', 
        save_name = 'Freund_fit_Nd_NL', x_label = 'log$(C_e[M])$', 
        y_label = 'log($Q_e$ [mol/kg$_{be}$])', Fit_type = 0)        
Fre_fitSm_NL = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Sm(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Sm(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Sm(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Sm(LR)'], folder_name = 'Fre_fit', 
        Color = Bent_color['Sar'], Title = 'Freundlich fit for Sm', 
        save_name = 'Freund_fit_Sm_NL', x_label = 'log$(C_e[M])$', 
        y_label = 'log($Q_e$ [mol/kg$_{be}$])', Fit_type = 0)
Fre_fitSm = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Sm(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Sm(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Sm(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Sm(LR)'], folder_name = 'Fre_fit', 
        Color = Bent_color['Sar'], Title = 'Freundlich fit for Sm', 
        save_name = 'Freund_fit_Sm', x_label = 'log$(C_e[M])$', 
        y_label = 'log($Q_e$ [mol/kg$_{be}$])')
        # r = 0.91, not bad though!
Fre_fitPu = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(['Data 1','Data 2','Data 3', 'Data 11'], axis =1).loc['Pu(LR)'], 
        Dict_Qe['< >'].iloc[:,3:-1].loc['Pu(LR)'], 
        Dict_el_sa_avg['std'].drop(['Data 1','Data 2','Data 3', 'Data 11'], axis =1).loc['Pu(LR)'], 
        Dict_Qe['std'].iloc[:,3:-1].loc['Pu(LR)'], folder_name = 'Fre_fit', 
        Color = Bent_color['Sar'], Title = 'Freundlich fit for Pu', 
        save_name = 'Freund_fit_Pu_11out', x_label = 'log$(C_e[M])$', 
        y_label = 'log($Q_e$ [mol/kg$_{be}$])')
        # r =  0.88 !
Fre_fitPu_NL = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(['Data 1','Data 2','Data 3', 'Data 11'], axis =1).loc['Pu(LR)'], 
        Dict_Qe['< >'].iloc[:,3:-1].loc['Pu(LR)'], 
        Dict_el_sa_avg['std'].drop(['Data 1','Data 2','Data 3', 'Data 11'], axis =1).loc['Pu(LR)'], 
        Dict_Qe['std'].iloc[:,3:-1].loc['Pu(LR)'], folder_name = 'Fre_fit', 
        Color = Bent_color['Sar'], Title = 'Freundlich fit for Pu', 
        save_name = 'Freund_fit_Pu_NL_11out', x_label = 'log$(C_e[M])$', 
        y_label = 'log($Q_e$ [mol/kg$_{be}$])', Fit_type= 0)
Fre_fitNp = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Np(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Np(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Np(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Np(LR)'], folder_name = 'Fre_fit', 
        Color = Bent_color['Sar'], Title = 'Freundlich fit for Np', 
        save_name = 'Freund_fit_Np', x_label = 'log$(C_e[M])$', 
        y_label = 'log($Q_e$ [mol/kg$_{be}$])')
            # r = 0.98!!
Fre_fitNp_NL = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Np(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Np(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Np(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Np(LR)'], folder_name = 'Fre_fit', 
        Color = Bent_color['Sar'], Title = 'Freundlich fit for Np', 
        save_name = 'Freund_fit_Np', x_label = 'log$(C_e[M])$', 
        y_label = 'log($Q_e$ [mol/kg$_{be}$])', Fit_type= 0)
Fre_fitAm = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Am(LR)'], 
         Dict_Qe['< >'].iloc[:,3:].loc['Am(LR)'], 
         Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Am(LR)'], 
         Dict_Qe['std'].iloc[:,3:].loc['Am(LR)'], folder_name = 'Fre_fit', 
         Color = Bent_color['Sar'], Title = 'Freundlich fit for Am', 
         save_name = 'Freund_fit_Am', x_label = 'log$(C_e[M])$', 
         y_label = 'log($Q_e$ [mol/kg$_{be}$])') 
          # r = 0.95!
Fre_fitAm_NL = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Am(LR)'], 
         Dict_Qe['< >'].iloc[:,3:].loc['Am(LR)'], 
         Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Am(LR)'], 
         Dict_Qe['std'].iloc[:,3:].loc['Am(LR)'], folder_name = 'Fre_fit', 
         Color = Bent_color['Sar'], Title = 'Freundlich fit for Am', 
         save_name = 'Freund_fit_Am_NL', x_label = 'log$(C_e[M])$', 
         y_label = 'log($Q_e$ [mol/kg$_{be}$])', Fit_type=	 0) 
Fre_fitCm = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Cm(LR)'], 
          Dict_Qe['< >'].iloc[:,3:].loc['Cm(LR)'], 
          Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Cm(LR)'], 
          Dict_Qe['std'].iloc[:,3:].loc['Cm(LR)'], folder_name = 'Fre_fit', 
          Color = Bent_color['Sar'], Title = 'Freundlich fit for Cm', 
          save_name = 'Freund_fit_Cm', x_label = 'log$(C_e[M])$', 
          y_label = 'log($Q_e$ [mol/kg$_{be}$])') 
         # r = 0.98!!
Fre_fitCm_NL = Read_JRC.Fre_fit(Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Cm(LR)'], 
          Dict_Qe['< >'].iloc[:,3:].loc['Cm(LR)'], 
          Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Cm(LR)'], 
          Dict_Qe['std'].iloc[:,3:].loc['Cm(LR)'], folder_name = 'Fre_fit', 
          Color = Bent_color['Sar'], Title = 'Freundlich fit for Cm', 
          save_name = 'Freund_fit_Cm_NL', x_label = 'log$(C_e[M])$', 
          y_label = 'log($Q_e$ [mol/kg$_{be}$])', Fit_type= 0) #No converge!
#---------- Iso case, Sr90
Fre_fitSr90 = Read_JRC.Fre_fit(
    Dict_el_sa_avg_iso['< >'].drop(['Data 1','Data 2','Data 3', 'Data 4',
        'Data 5', 'Data 6', 'Data 7', 'Data 8'], axis =1).loc['Sr90(LR)'], 
          Dict_Qe_iso['< >'].drop(['Data 1','Data 2','Data 3', 'Data 4',
              'Data 5', 'Data 6', 'Data 7', 'Data 8'], axis =1).loc['Sr90(LR)'], 
          Dict_el_sa_avg_iso['std'].drop(['Data 1','Data 2','Data 3', 'Data 4',
              'Data 5', 'Data 6', 'Data 7', 'Data 8'], axis =1).loc['Sr90(LR)'], 
          Dict_Qe_iso['std'].drop(['Data 1','Data 2','Data 3', 'Data 4',
              'Data 5', 'Data 6', 'Data 7', 'Data 8'], axis =1).loc['Sr90(LR)'], 
   folder_name = 'Fre_fit', 
          Color = Bent_color['Sar'], Title = 'Freundlich fit for Sr90', 
          save_name = 'Freund_fit_Sr90_891011', x_label = 'log$(C_e[M])$', 
          y_label = 'log($Q_e$ [mol/kg$_{be}$])') 
Fre_fitSr90_NL = Read_JRC.Fre_fit(
    Dict_el_sa_avg_iso['< >'].drop(['Data 1','Data 2','Data 3', 'Data 4',
        'Data 5', 'Data 6', 'Data 7', 'Data 8'], axis =1).loc['Sr90(LR)'], 
          Dict_Qe_iso['< >'].drop(['Data 1','Data 2','Data 3', 'Data 4',
              'Data 5', 'Data 6', 'Data 7', 'Data 8'], axis =1).loc['Sr90(LR)'], 
          Dict_el_sa_avg_iso['std'].drop(['Data 1','Data 2','Data 3', 'Data 4',
              'Data 5', 'Data 6', 'Data 7', 'Data 8'], axis =1).loc['Sr90(LR)'], 
          Dict_Qe_iso['std'].drop(['Data 1','Data 2','Data 3', 'Data 4',
              'Data 5', 'Data 6', 'Data 7', 'Data 8'], axis =1).loc['Sr90(LR)'], 
   folder_name = 'Fre_fit', 
          Color = Bent_color['Sar'], Title = 'Freundlich fit for Sr90', 
          save_name = 'Freund_fit_Sr90_891011', x_label = 'log$(C_e[M])$', 
          y_label = 'log($Q_e$ [mol/kg$_{be}$])', Fit_type= 0) 

np.log10(Dict_Qe_iso['< >'].drop(S_rem,  axis = 1
                                     ).loc['Sr90(LR)']).iloc[4:]

Fre_fits= pd.DataFrame({'U': Fre_fitU, 'U_NL': Fre_fitU_NL,'Cs' : Fre_fitCs, 
    'Cs_NL' : Fre_fitCs_NL, 'La' : Fre_fitLa,'La_NL' : Fre_fitLa_NL,
    'Nd' : Fre_fitNd,'Nd_NL' : Fre_fitNd_NL,'Sm' : Fre_fitSm,'Sm_NL' : Fre_fitSm_NL,
   'Pu': Fre_fitPu,'Pu_NL': Fre_fitPu_NL, 'Np' : Fre_fitNp,'Np_NL' : Fre_fitNp_NL,
    'Am' : Fre_fitAm,'Am_NL' : Fre_fitAm_NL,'Cm' : Fre_fitCm,'Cm_NL' : Fre_fitCm_NL,
    'Sr90' : Fre_fitSr90,'Sr90_NL' : Fre_fitSr90_NL} )               
               #gathering all together!


del (Fre_fitCs, Fre_fitCs_NL, Fre_fitU, Fre_fitU_NL, Fre_fitLa, Fre_fitLa_NL, 
     Fre_fitNd, Fre_fitNd_NL,Fre_fitSm, Fre_fitSm_NL,
     Fre_fitPu, Fre_fitPu_NL,Fre_fitNp, Fre_fitNp_NL,Fre_fitAm, Fre_fitAm_NL,
     Fre_fitCm, Fre_fitCm_NL, Fre_fitSr90, Fre_fitSr90_NL)  
                                            #deleting now those variables

print('##### Freundlich fits #######')
print(Fre_fits)
print('#################')


'''
Comparing all of them:
    
    Huge %rsd for KF! for NL of U and Cs not anymore! 
    So compare always both versions!!
    
for Pu different behaviour, n < 1!! rest n >1! could be ebcause of the last value!!

Cm has n = 2!! 1st time I see sugh high n kliao!
'''

#----------         Fit = iso plot -----------------
"""
This would be a really nice thing, for reporting and papers!

I could do 2 plots:
    .Lanth + Cs
    .Actini
    
Which Fre fit to inlcude in the Fre fit? Well the best, which will be chosen 
based on:
    .low %rsd for both n and K
    .high r
    
By looking at the variable, I decide (I compared U, not so big differences xD)

U: NL
Cs: NL
La:line
Nd: lin
Sm lin
            :lanth seems from param that NL, but plot is worse, lin better!
Pu NL
Np lin
Am lin
Cm lin

Easier than diong a code xDDD
"""

#tO PLOT THE fRE FIT, LETS DEFINE THE FUNCTION:
def Fre_fit_eq(C,K, n):         #fre fit eq, non linear (original)
    return K*C**(1/n)

#           Fre fit U, Lin vs NL

#------------- Fre fit plot + ads iso Actinides
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title('Adsorption isotherm + Fre fit for the Ad CL experiment',  fontsize=22, wrap=True)           #title
#### U
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['U(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['U(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['U(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['U(LR)'], 
            fmt = 'o', label = 'U', color = '#4682B4')
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)'])) ),
    1/Fre_fits['U']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'
    ].iloc[:,3:].loc['U(LR)']),
    np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']))
    )+ np.log10( Fre_fits['U']['K[L^n/(kg*mol^{n-1})]']),
    fmt = '--', color = '#4682B4', label = "Lin Fre fit") #Fit
        #note for the len I can use the df with std or without, since size is the same!
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)'])) ),
    1/Fre_fits['U_NL']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'
    ].iloc[:,3:].loc['U(LR)']),
    np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']))
    )+ np.log10( Fre_fits['U_NL']['K[L^n/(kg*mol^{n-1})]']),
    fmt = '-', color = '#4682B4', label = "NL Fre fit") #Fit
plt.ylabel('log($Q_e$ [mol/kg$_{ben}$])', fontsize= Font)              #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font )
#plt.yscale('log')
#plt.xscale('log')
plt.legend(fontsize = Font)
plt.tick_params(axis='both', labelsize= Font) 
plt.minorticks_on()             #enabling minor grid lines
plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        
plt.grid(which = 'major')
plt.savefig('Ads_iso_Fre_fit_U_no123.png', format='png', bbox_inches='tight')
plt.show()  
'''
Okay, this was to compare the NL vs Lin Fre fit, for U. The plot look farily
similar. Consideing the paameters, they are better for the NL. But well, I
actually could choose the best, why shouldnt I do it bro? MF, you have programming
knowledge, use it mf!!
'''

#------------- Fre fit plot + ads iso Actinides
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title('Adsorption isotherm + Fre fit for the Ad CL experiment',  fontsize=22, wrap=True)           #title
#### U
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['U(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['U(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['U(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['U(LR)'], 
            fmt = 'o', label = 'U', color = '#4682B4')
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)'])) ),
    1/Fre_fits['U_NL']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'
    ].iloc[:,3:].loc['U(LR)']),
    np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']))
    )+ np.log10( Fre_fits['U_NL']['K[L^n/(kg*mol^{n-1})]']),
    fmt = '--', color = '#4682B4', label = None)                        #Fit
        #note for the len I can use the df with std or without, since size is the same!
#### Np
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Np(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Np(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Np(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Np(LR)'], 
            fmt = 'o', label = 'Np', color = '#228B22')
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Np(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Np(LR)'])) ),
    1/Fre_fits['Np']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Np(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Np(LR)']))
            )+ np.log10( Fre_fits['Np']['K[L^n/(kg*mol^{n-1})]']),
                fmt = '--', color = '#228B22', label = None) #Fit
#### Pu
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Pu(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Pu(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Pu(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Pu(LR)'], 
            fmt = 'o', label = 'Pu', color = '#DC143C')
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)'])) ),
    1/Fre_fits['Pu_NL']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)']))
            )+ np.log10( Fre_fits['Pu_NL']['K[L^n/(kg*mol^{n-1})]']),
                fmt = '--', color = '#DC143C', label = None) #Fit
#### Am
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Am(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Am(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Am(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Am(LR)'], 
            fmt = 'o', label = 'Am', color = '#FF8C00')
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)'])) ),
    1/Fre_fits['Am']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)']))
            )+ np.log10( Fre_fits['Am']['K[L^n/(kg*mol^{n-1})]']),
                fmt = '--', color = '#FF8C00', label = None) #Fit
#### Cm
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Cm(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Cm(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Cm(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Cm(LR)'], 
            fmt = 'o', label = 'Cm', color = '#9932CC')
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)'])) ),
    1/Fre_fits['Cm']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)']))
            )+ np.log10( Fre_fits['Cm']['K[L^n/(kg*mol^{n-1})]']),
                fmt = '--', color = '#9932CC', label = None) #Fit
plt.ylabel('log($Q_e$ [mol/kg$_{ben}$])', fontsize= Font)              #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font )
#plt.yscale('log')
#plt.xscale('log')
plt.legend(fontsize = Font)
plt.tick_params(axis='both', labelsize= Font) 
plt.minorticks_on()             #enabling minor grid lines
plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        
plt.grid(which = 'major')
plt.savefig('Ads_iso_Fre_fit_acti_no123.png', format='png', bbox_inches='tight')
plt.show()  



#------------- Fre fit plot + ads iso lanth + Cs
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title('Adsorption isotherm + Fre fit for the Ad CL experiment',  fontsize=22, wrap=True)           #title
#### Cs
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Cs(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Cs(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Cs(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Cs(LR)'], 
            fmt = 'o', label = 'Cs', color = '#708090')
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cs(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cs(LR)'])) ),
    1/Fre_fits['Cs_NL']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'
    ].iloc[:,3:].loc['Cs(LR)']),
    np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cs(LR)']))
    )+ np.log10( Fre_fits['Cs_NL']['K[L^n/(kg*mol^{n-1})]']),
    fmt = '--', color = '#708090', label = None) #Fit NO label!
        #note for the len I can use the df with std or without, since size is the same!
#### La
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['La(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['La(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['La(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['La(LR)'], 
            fmt = 'o', label = 'La', color = '#DAA520')
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['La(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['La(LR)'])) ),
    1/Fre_fits['La']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['La(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['La(LR)']))
            )+ np.log10( Fre_fits['La']['K[L^n/(kg*mol^{n-1})]']),
                fmt = '--', color = '#DAA520', label = None) #Fit
#### Nd
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Nd(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Nd(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Nd(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Nd(LR)'], 
            fmt = 'o', label = 'Nd', color = '#FF69B4')
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Nd(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Nd(LR)'])) ),
    1/Fre_fits['Nd']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Nd(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Nd(LR)']))
            )+ np.log10( Fre_fits['Nd']['K[L^n/(kg*mol^{n-1})]']),
                fmt = '--', color = '#FF69B4', label = None) #Fit
# #### Sr90
# plt.errorbar(np.log10(Dict_el_sa_avg_iso['< >']).iloc[:,8:].loc['Sr90(LR)'], 
#      np.log10(Dict_Qe_iso['< >']).iloc[:,8:].loc['Sr90(LR)'],
#     yerr = np.abs(Dict_Qe_iso['std'] / 
#             (Dict_Qe_iso['< >']*np.log(10))).iloc[:,8:].loc['Sr90(LR)'], 
#     xerr = np.abs(Dict_el_sa_avg_iso['std']/
#             (Dict_el_sa_avg_iso['< >']*np.log(10))).iloc[:,8:].loc['Sr90(LR)'], 
#             fmt = 'o', label = 'Sr90', color = 'b')
# plt.errorbar(
#     np.log10(np.linspace(np.min(Dict_el_sa_avg_iso['< >'].iloc[:,8:].loc['Sr90(LR)']),
#             np.max(Dict_el_sa_avg_iso['< >'].iloc[:,8:].loc['Sr90(LR)'])) ),
#     1/Fre_fits['Sr90']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg_iso['< >'
#     ].iloc[:,8:].loc['Sr90(LR)']),
#     np.max(Dict_el_sa_avg_iso['< >'].iloc[:,8:].loc['Sr90(LR)']))
#     )+ np.log10( Fre_fits['Sr90']['K[L^n/(kg*mol^{n-1})]']),
#     fmt = '--', label = None, color = 'b') #Fit NO label!
# #### Sm
# plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Sm(LR)'], 
#      np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Sm(LR)'],
#     yerr = np.abs(Dict_Qe['std'] / 
#             (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Sm(LR)'], 
#     xerr = np.abs(Dict_el_sa_avg['std']/
#             (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Sm(LR)'], 
#             fmt = 'o', label = 'Sm', color = '#32CD32')
# plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Sm(LR)']),
#             np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Sm(LR)'])) ),
#     1/Fre_fits['Sm']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Sm(LR)']),
#             np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Sm(LR)']))
#             )+ np.log10( Fre_fits['Sm']['K[L^n/(kg*mol^{n-1})]']),
#                fmt = '--', color = , label = None) #Fit
plt.ylabel('log($Q_e$ [mol/kg$_{ben}$])', fontsize= Font)              #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font )
#plt.yscale('log')
#plt.xscale('log')
plt.legend(fontsize = Font)
plt.tick_params(axis='both', labelsize= Font) 
plt.minorticks_on()             #enabling minor grid lines
plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        
plt.grid(which = 'major')
plt.savefig('Ads_iso_Fre_fit_lanthCs_no123.png', format='png', bbox_inches='tight')
plt.show()  


#------------- Fre fit + ads iso for Sr 90!
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title('Adsorption isotherm + Fre fit for the Ad CL experiment',  fontsize=22, wrap=True)           #title
#### Sr90
plt.errorbar(np.log10(Dict_el_sa_avg_iso['< >']).iloc[:,8:].loc['Sr90(LR)'], 
     np.log10(Dict_Qe_iso['< >']).iloc[:,8:].loc['Sr90(LR)'],
    yerr = np.abs(Dict_Qe_iso['std'] / 
            (Dict_Qe_iso['< >']*np.log(10))).iloc[:,8:].loc['Sr90(LR)'], 
    xerr = np.abs(Dict_el_sa_avg_iso['std']/
            (Dict_el_sa_avg_iso['< >']*np.log(10))).iloc[:,8:].loc['Sr90(LR)'], 
            fmt = 'o', label = 'Cs', color = 'b')
plt.errorbar(
    np.log10(np.linspace(np.min(Dict_el_sa_avg_iso['< >'].iloc[:,8:].loc['Sr90(LR)']),
            np.max(Dict_el_sa_avg_iso['< >'].iloc[:,8:].loc['Sr90(LR)'])) ),
    1/Fre_fits['Sr90']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg_iso['< >'
    ].iloc[:,8:].loc['Sr90(LR)']),
    np.max(Dict_el_sa_avg_iso['< >'].iloc[:,8:].loc['Sr90(LR)']))
    )+ np.log10( Fre_fits['Sr90']['K[L^n/(kg*mol^{n-1})]']),
    fmt = '--', label = None, color = 'b') #Fit NO label!
        #note for the len I can use the df with std or without, since size is the same!
plt.ylabel('log($Q_e$ [mol/kg$_{ben}$])', fontsize= Font)              #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font)
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(fontsize = Font)
plt.tick_params(axis='both', labelsize= Font) 
plt.minorticks_on()             #enabling minor grid lines
plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        
plt.grid(which = 'major')
plt.savefig('Ads_iso_Fre_fit_Sr90_91011.png', format='png', bbox_inches='tight')
plt.show()  
'''
Okay, this plot I might include in the report. I put 4 datasets. I might
consider to remove the 1st one (data 8), due to the massive errorbar

Furtherworktodohere!!!!!
'''

#----------------------------------
#%%             Lang fit
#--------------------------------\


'''
Remember that for the Langmuir fit, there are 2 linearizations. We will try
both of them actually

We wil try only for La, since from the Qe iso plots, we saw that it could be
Langmuir actually!

Maybe for U also works, since Kd decrease, a sign of limit in sorption!
'''


Lang_fitLa_1 = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['La(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['La(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['La(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['La(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir lin1 fit for La', 
        save_name = 'lang_fit_La_1', x_label = '$C_e [M]$', 
             y_label = '$C_e/Q_e$ [kg$_{be}$/L]', Fit_type= 1)
        # r = 0.96 ! better than Freund bro, nice!
Lang_fitLa_2 = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['La(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['La(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['La(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['La(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir lin2 fit for La', 
        save_name = 'lang_fit_La_2', x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$Q_e/C_e$ [kg$_{be}$/L]', Fit_type= 2)
        # not so good, r= 0.5, so we focus on 1 xD
        
Lang_fitLa_NL = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['La(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['La(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['La(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['La(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir fit for La', 
        save_name = 'lang_fit_La_NL', x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$Q_e/C_e$ [kg$_{be}$/L]', Fit_type= 0)
            #Really good fit, r = 0.98!
#Nd
Lang_fitNd_1 = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Nd(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Nd(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Nd(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Nd(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir lin1 fit for Nd', 
        save_name = 'lang_fit_Nd_1_11out',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$C_e/Q_e$ [kg$_{be}$/L]', Fit_type= 1)
            #bad fit, r = 0.31!!
Lang_fitNd_2 = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Nd(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Nd(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Nd(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Nd(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir lin2 fit for Nd', 
        save_name = 'lang_fit_Nd_2_11out',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$Q_e/C_e$ [kg$_{be}$/L]', Fit_type= 2)
                #horrible, r = 0.2!
Lang_fitNd_NL = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Nd(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Nd(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Nd(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Nd(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir fit for Nd', 
        save_name = 'lang_fit_Nd_NL_11out',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$Q_e/C_e$ [kg$_{be}$/L]', Fit_type= 0)
#Pu
Lang_fitPu_1 = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(['Data 1','Data 2','Data 3', 'Data 11'], axis =1).loc['Pu(LR)'], 
        Dict_Qe['< >'].iloc[:,3:-1].loc['Pu(LR)'], 
        Dict_el_sa_avg['std'].drop(['Data 1','Data 2','Data 3', 'Data 11'], axis =1).loc['Pu(LR)'], 
        Dict_Qe['std'].iloc[:,3:-1].loc['Pu(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir lin1 fit for Pu', 
        save_name = 'lang_fit_Pu_1_11out',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$C_e/Q_e$ [kg$_{be}$/L]', Fit_type= 1)
            #bad fit, r = 0.31!!
Lang_fitPu_2 = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(['Data 1','Data 2','Data 3', 'Data 11'], axis =1).loc['Pu(LR)'], 
        Dict_Qe['< >'].iloc[:,3:-1].loc['Pu(LR)'], 
        Dict_el_sa_avg['std'].drop(['Data 1','Data 2','Data 3', 'Data 11'], axis =1).loc['Pu(LR)'], 
        Dict_Qe['std'].iloc[:,3:-1].loc['Pu(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir lin2 fit for Pu', 
        save_name = 'lang_fit_Pu_2_11out',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$Q_e/C_e$ [kg$_{be}$/L]', Fit_type= 2)
                #horrible, r = 0.2!
Lang_fitPu_NL = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(['Data 1','Data 2','Data 3', 'Data 11'], axis =1).loc['Pu(LR)'], 
        Dict_Qe['< >'].iloc[:,3:-1].loc['Pu(LR)'], 
        Dict_el_sa_avg['std'].drop(['Data 1','Data 2','Data 3', 'Data 11'], axis =1).loc['Pu(LR)'], 
        Dict_Qe['std'].iloc[:,3:-1].loc['Pu(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir fit for Pu', 
        save_name = 'lang_fit_Pu_NL_11out',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$Q_e/C_e$ [kg$_{be}$/L]', Fit_type= 0)
        #not so ggod, good r, 0.92, but horrible values!
#U
Lang_fitU_1 = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['U(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['U(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['U(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['U(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir lin1 fit for U', 
        save_name = 'lang_fit_U_1', x_label = '$C_e [M]$', 
             y_label = '$C_e/Q_e$ [kg$_{be}$/L]', Fit_type= 1)
            #r = 0.89!!
Lang_fitU_2 = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['U(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['U(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['U(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['U(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir lin2 fit for U', 
        save_name = 'lang_fit_U_2',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$Q_e/C_e$ [kg$_{be}$/L]', Fit_type= 2)                
            #r = 0.56
            
Lang_fitU_NL = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['U(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['U(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['U(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['U(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir fit for U', 
        save_name = 'lang_fit_U_NL',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$Q_e/C_e$ [kg$_{be}$/L]', Fit_type= 0) 
        # r = 0.942!!!

#Cm
Lang_fitCm_NL = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Cm(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Cm(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Cm(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Cm(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir fit for Cm', 
        save_name = 'lang_fit_Cm_NL',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$Q_e/C_e$ [kg$_{be}$/L]', Fit_type= 0)
                #r = 0.905, without last sample possibly better!!!
Lang_fitCm_1 = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Cm(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Cm(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Cm(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Cm(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir lin1 fit for Cm', 
        save_name = 'lang_fit_Cm_1',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$C_e/Q_e$ [kg$_{be}$/L]', Fit_type= 1)
        #R = 0.857! could be good
Lang_fitCm_2 = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Cm(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Cm(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Cm(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Cm(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir lin2 fit for Cm', 
        save_name = 'lang_fit_Cm_2',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$Q_e/C_e$ [kg$_{be}$/L]', Fit_type= 2)
#Am
Lang_fitAm_NL = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Am(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Am(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Am(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Am(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir fit for Am', 
        save_name = 'lang_fit_Am_NL',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$Q_e/C_e$ [kg$_{be}$/L]', Fit_type= 0)
                #r = 0.905, without last sample possibly better!!!
Lang_fitAm_1 = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Am(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Am(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Am(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Am(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir lin1 fit for Am', 
        save_name = 'lang_fit_Am_1',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$C_e/Q_e$ [kg$_{be}$/L]', Fit_type= 1)
        #R = 0.857! could be good
Lang_fitAm_2 = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Am(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Am(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Am(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Am(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir lin2 fit for Am', 
        save_name = 'lang_fit_Am_2',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$Q_e/C_e$ [kg$_{be}$/L]', Fit_type= 2)
Lang_fitNp_NL = Read_JRC.Lang_fit(
        Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Np(LR)'], 
        Dict_Qe['< >'].iloc[:,3:].loc['Np(LR)'], 
        Dict_el_sa_avg['std'].drop(S_rem, axis =1).loc['Np(LR)'], 
        Dict_Qe['std'].iloc[:,3:].loc['Np(LR)'], folder_name = 'Lang_fit', 
        Color = Bent_color['Sar'], Title = 'Langmuir fit for Np', 
        save_name = 'lang_fit_Np_NL',x_label = '$Q_e$ [mol/kg$_{be}$]', 
             y_label = '$Q_e/C_e$ [kg$_{be}$/L]', Fit_type= 0)

Lang_fits= pd.DataFrame({'La_1': Lang_fitLa_1, 'La_2' : Lang_fitLa_2, 
        'La_NL' : Lang_fitLa_NL,
        'Nd_1': Lang_fitNd_1, 'Nd_2' : Lang_fitNd_2, 
                'Nd_NL' : Lang_fitNd_NL,
        'U_1': Lang_fitU_1, 'U_2' : Lang_fitU_2, 'U_NL' : Lang_fitU_NL,  
        'Pu_1': Lang_fitPu_1, 'Pu_2' : Lang_fitPu_2, 'Pu_NL' : Lang_fitPu_NL,
        'Np_NL' : Lang_fitNp_NL,
        'Am_1' : Lang_fitAm_1,'Am_2' : Lang_fitAm_2,'Am_NL' : Lang_fitAm_NL,
        'Cm_1' : Lang_fitCm_1,'Cm_2' : Lang_fitCm_2,'Cm_NL' : Lang_fitCm_NL} ) 
                         #gathering all together!
del (Lang_fitLa_1, Lang_fitLa_2, Lang_fitLa_NL, 
     Lang_fitNd_1, Lang_fitNd_2, Lang_fitNd_NL, Lang_fitU_NL, Lang_fitU_1,
     Lang_fitU_2,Lang_fitAm_NL, Lang_fitAm_1, Lang_fitAm_2, 
     Lang_fitCm_1,Lang_fitCm_2, Lang_fitCm_NL,
     Lang_fitNp_NL,
     Lang_fitPu_NL, Lang_fitPu_1, Lang_fitPu_2) 


''' ------------ Fit discussion -------
Take care with non linear fits!! because r2 is misleading in quadratic fits, 
 and use also %rsd as metric! Considering that, NL fit for U not so nice!
##Fre
U, Np, Pu not so so clear, since they do a kind of fluctuating pattern, not so
clear


U is Freundl, as has beeen always r 0.98 for NL, 0.92 for lin fit, best pararm NL!
Cs also,
La good, r 0.9 lin, 0.91 NL, NL best param!
Pu r = 0.88 lin, NL worse paramteres. But, if data 11 removed, no converge!
Np, r = 0.98 lin, really good! 
Am r = 0.95 linm the best (r and %rsd)
Cm, lin r = 0.98

### Lang

La, La 1 r = 0.96, good!
U, U1 o NL best, r = 0.89 lin, 
Pu bad fit. But, removing the data 11, NL fit absolutely great!
Np bad fit
Am ciould be, Am 1 0.86, but parameters meh also does not look good!
Cm could be, but does not look good, so discarding it!

La here Langmuir better, but for the cold experiments, langmuir was also
applicable for the CL ==> effect of having other lanthanides?



WOuld be nice to show this in plot version, but how?
    i) FOr each element, Fre/Lang
    ii) For 3/4 elements, Fre, and appart, Lang
                                                             
What could be the best?
'''


#           Lang fit + ads iso
"""
Okay, this could me more tricky.

Well, We kinda agree that La, Pu have really lang fit, so I could try to plot it


Q_e = Q_max * K_L* C_e/(1+K_L* C_e) 


"""

#TO apply the function, lets define the Lang eq:
def Lang_fit_eq(Ce, Qmax, K):      #equaition fo langmuir fit
            #Q_e = Q_max * K_L/(1+K_L* C_e) * C_e
    return Qmax * K * Ce /(1+K*Ce) 

#------------- Lang fit plot + ads iso Pu, U, La, Nd
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title('Adsorption isotherm + Fre fit for the Ad CL experiment',  fontsize=22, wrap=True)           #title
#### La
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['La(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['La(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['La(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['La(LR)'], 
            fmt = 'o', label = 'La', color = '#DAA520')
plt.errorbar(
    np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['La(LR)']),
                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['La(LR)'])) ),
    np.log10(
   Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['La(LR)']),
               np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['La(LR)'])),
    Lang_fits['La_1']['Q_max[mol/kg_be]'], Lang_fits['La_1']['K_L[1/M]']) ),
                fmt = '--', color = '#DAA520', label = None)   # Lang Fit
#### Nd
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Nd(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Nd(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Nd(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Nd(LR)'], 
            fmt = 'o', label = 'Nd', color = "#FF69B4")
plt.errorbar(
    np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Nd(LR)']),
                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Nd(LR)'])) ),
    np.log10(
   Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Nd(LR)']),
               np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Nd(LR)'])),
    Lang_fits['Nd_1']['Q_max[mol/kg_be]'], Lang_fits['Nd_1']['K_L[1/M]']) ),
                fmt = '--', color = "#FF69B4", label = None)           #Lang Fit
#### U
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['U(LR)'], 
      np.log10(Dict_Qe['< >']).iloc[:,3:].loc['U(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['U(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['U(LR)'], 
            fmt = 'o', label = 'U', color = '#4682B4')
plt.errorbar(
    np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']),
                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)'])) ),
    np.log10(
    Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']),
                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)'])),
    Lang_fits['U_1']['Q_max[mol/kg_be]'], Lang_fits['U_1']['K_L[1/M]']) ),
                fmt = '--', color = '#4682B4', label = None)           #Lang Fit
#### Pu
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Pu(LR)'], 
      np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Pu(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Pu(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Pu(LR)'], 
            fmt = 'o', label = 'Pu', color = '#DC143C')
plt.errorbar(
    np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)']),
                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)'])) ),
    np.log10(
    Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)']),
                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)'])),
    Lang_fits['Pu_1']['Q_max[mol/kg_be]'], Lang_fits['Pu_1']['K_L[1/M]']) ),
                fmt = '--', color = '#DC143C', label = None)           #Lang Fit
# #### Am
# plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Am(LR)'], 
#      np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Am(LR)'],
#     yerr = np.abs(Dict_Qe['std'] / 
#             (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Am(LR)'], 
#     xerr = np.abs(Dict_el_sa_avg['std']/
#             (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Am(LR)'], 
#             fmt = 'o', label = 'Am', color = '#FF8C00')
# plt.errorbar(
#     np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)']),
#                 np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)'])) ),
#     np.log10(
#    Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)']),
#                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)'])),
#     Lang_fits['Am_1']['Q_max[mol/kg_be]'], Lang_fits['Am_1']['K_L[1/M]']) ),
#                 fmt = '--', color = '#FF8C00', label = None)           #Lang Fit
# #### Cm
# plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Cm(LR)'], 
#      np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Cm(LR)'],
#     yerr = np.abs(Dict_Qe['std'] / 
#             (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Cm(LR)'], 
#     xerr = np.abs(Dict_el_sa_avg['std']/
#             (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Cm(LR)'], 
#             fmt = 'o', label = 'Cm', color = '#4682B4')
# plt.errorbar(
#     np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)']),
#                 np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)'])) ),
#     np.log10(
#    Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)']),
#                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)'])),
#     Lang_fits['Cm_1']['Q_max[mol/kg_be]'], Lang_fits['Cm_1']['K_L[1/M]']) ),
#                 fmt = '--', color = '#4682B4', label = None)           #Lang Fit

plt.ylabel('log($Q_e$ [mol/kg$_{ben}$])', fontsize= Font)              #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font )
#plt.yscale('log')
#plt.xscale('log')
plt.legend(fontsize = Font)
plt.tick_params(axis='both', labelsize= Font) 
plt.minorticks_on()             #enabling minor grid lines
plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        
plt.grid(which = 'major')
plt.savefig('Ads_iso_Lang_fit_LaNdPuU_no123.png', format='png', bbox_inches='tight')
plt.show()  

#------------- Lang fit plot + ads iso La + Nd
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title('Adsorption isotherm + Fre fit for the Ad CL experiment',  fontsize=22, wrap=True)           #title
#### La
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['La(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['La(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['La(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['La(LR)'], 
            fmt = 'o', label = 'La', color = '#DAA520')
plt.errorbar(
    np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Sm(LR)']),
                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Sm(LR)'])) ),
    np.log10(
   Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['La(LR)']),
               np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['La(LR)'])),
    Lang_fits['La_1']['Q_max[mol/kg_be]'], Lang_fits['La_1']['K_L[1/M]']) ),
                fmt = '--', color = '#DAA520', label = None)   # Lang Fit
#### Nd
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Nd(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Nd(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Nd(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Nd(LR)'], 
            fmt = 'o', label = 'Nd', color = "#FF69B4")
plt.errorbar(
    np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Nd(LR)']),
                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Nd(LR)'])) ),
    np.log10(
   Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Nd(LR)']),
               np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Nd(LR)'])),
    Lang_fits['Nd_1']['Q_max[mol/kg_be]'], Lang_fits['Nd_1']['K_L[1/M]']) ),
                fmt = '--', color = "#FF69B4", label = None)           #Lang Fit
# #### Pu
# plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Pu(LR)'], 
#      np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Pu(LR)'],
#     yerr = np.abs(Dict_Qe['std'] / 
#             (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Pu(LR)'], 
#     xerr = np.abs(Dict_el_sa_avg['std']/
#             (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Pu(LR)'], 
#             fmt = 'o', label = 'Pu', color = '#DC143C')
# plt.errorbar(
#     np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)']),
#                 np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)'])) ),
#     np.log10(
#    Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)']),
#                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)'])),
#     Lang_fits['Pu_1']['Q_max[mol/kg_be]'], Lang_fits['Pu_1']['K_L[1/M]']) ),
#                 fmt = '--', color = '#DC143C', label = None)           #Lang Fit
# #### U
# plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['U(LR)'], 
#      np.log10(Dict_Qe['< >']).iloc[:,3:].loc['U(LR)'],
#     yerr = np.abs(Dict_Qe['std'] / 
#             (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['U(LR)'], 
#     xerr = np.abs(Dict_el_sa_avg['std']/
#             (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['U(LR)'], 
#             fmt = 'o', label = 'U', color = '#4682B4')
# plt.errorbar(
#     np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']),
#                 np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)'])) ),
#     np.log10(
#    Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']),
#                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)'])),
#     Lang_fits['U_1']['Q_max[mol/kg_be]'], Lang_fits['U_1']['K_L[1/M]']) ),
#                 fmt = '--', color = '#4682B4', label = None)           #Lang Fit
# #### Am
# plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Am(LR)'], 
#      np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Am(LR)'],
#     yerr = np.abs(Dict_Qe['std'] / 
#             (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Am(LR)'], 
#     xerr = np.abs(Dict_el_sa_avg['std']/
#             (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Am(LR)'], 
#             fmt = 'o', label = 'Am', color = '#FF8C00')
# plt.errorbar(
#     np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)']),
#                 np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)'])) ),
#     np.log10(
#    Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)']),
#                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)'])),
#     Lang_fits['Am_1']['Q_max[mol/kg_be]'], Lang_fits['Am_1']['K_L[1/M]']) ),
#                 fmt = '--', color = '#FF8C00', label = None)           #Lang Fit
# #### Cm
# plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Cm(LR)'], 
#      np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Cm(LR)'],
#     yerr = np.abs(Dict_Qe['std'] / 
#             (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Cm(LR)'], 
#     xerr = np.abs(Dict_el_sa_avg['std']/
#             (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Cm(LR)'], 
#             fmt = 'o', label = 'Cm', color = '#4682B4')
# plt.errorbar(
#     np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)']),
#                 np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)'])) ),
#     np.log10(
#    Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)']),
#                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)'])),
#     Lang_fits['Cm_1']['Q_max[mol/kg_be]'], Lang_fits['Cm_1']['K_L[1/M]']) ),
#                 fmt = '--', color = '#4682B4', label = None)           #Lang Fit

plt.ylabel('log($Q_e$ [mol/kg$_{ben}$])', fontsize= Font)              #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font )
#plt.yscale('log')
#plt.xscale('log')
plt.legend(fontsize = Font)
plt.tick_params(axis='both', labelsize= Font) 
plt.minorticks_on()             #enabling minor grid lines
plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        
plt.grid(which = 'major')
plt.savefig('Ads_iso_Lang_fit_LaNd_no123.png', format='png', bbox_inches='tight')
plt.show()  

#------------- Lang fit plot + ads iso La + Nd
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title('Adsorption isotherm + Fre fit for the Ad CL experiment',  fontsize=22, wrap=True)           #title

#### Pu
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Pu(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Pu(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Pu(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Pu(LR)'], 
            fmt = 'o', label = 'Pu', color = '#DC143C')
plt.errorbar(
    np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)']),
                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)'])) ),
    np.log10(
   Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)']),
               np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)'])),
    Lang_fits['Pu_1']['Q_max[mol/kg_be]'], Lang_fits['Pu_1']['K_L[1/M]']) ),
                fmt = '--', color = '#DC143C', label = None)           #Lang Fit
#### U
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['U(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['U(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['U(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['U(LR)'], 
            fmt = 'o', label = 'U', color = '#4682B4')
plt.errorbar(
    np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']),
                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)'])) ),
    np.log10(
   Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']),
               np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)'])),
    Lang_fits['U_1']['Q_max[mol/kg_be]'], Lang_fits['U_1']['K_L[1/M]']) ),
                fmt = '--', color = '#4682B4', label = None)           #Lang Fit
#### Am
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Am(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Am(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Am(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Am(LR)'], 
            fmt = 'o', label = 'Am', color = '#FF8C00')
plt.errorbar(
    np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)']),
                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)'])) ),
    np.log10(
   Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)']),
               np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)'])),
    Lang_fits['Am_1']['Q_max[mol/kg_be]'], Lang_fits['Am_1']['K_L[1/M]']) ),
                fmt = '--', color = '#FF8C00', label = None)           #Lang Fit
#### Cm
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Cm(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Cm(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Cm(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Cm(LR)'], 
            fmt = 'o', label = 'Cm', color = '#9932CC')
plt.errorbar(
    np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)']),
                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)'])) ),
    np.log10(
   Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)']),
               np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)'])),
    Lang_fits['Cm_1']['Q_max[mol/kg_be]'], Lang_fits['Cm_1']['K_L[1/M]']) ),
                fmt = '--', color = '#9932CC', label = None)           #Lang Fit

plt.ylabel('log($Q_e$ [mol/kg$_{ben}$])', fontsize= Font)              #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font )
#plt.yscale('log')
#plt.xscale('log')
plt.legend(fontsize = Font)
plt.tick_params(axis='both', labelsize= Font) 
plt.minorticks_on()             #enabling minor grid lines
plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        
plt.grid(which = 'major')
plt.savefig('Ads_iso_Lang_fit_actin_no123.png', format='png', bbox_inches='tight')
plt.show()  
'''
This is a really good plot. WE can see clearly how Pu and La are really will 
fitted via Langmuir model. Ofc, last Pu sample was removed

For U, more or less. Still Fre is not a perfect fit, but from parameters Fre
is slightly better.
'''


#               Isotherms  plots + Lang and Fre fit together

#------------- Fre fit plot + ads iso Actinides
plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
#plt.title('Adsorption isotherm + Fre fit for the Ad CL experiment',  fontsize=22, wrap=True)           #title
#### U
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['U(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['U(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['U(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['U(LR)'], 
            fmt = 'o', label = 'U', color = '#4682B4')
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)'])) ),
    1/Fre_fits['U_NL']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'
    ].iloc[:,3:].loc['U(LR)']),
    np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']))
    )+ np.log10( Fre_fits['U_NL']['K[L^n/(kg*mol^{n-1})]']),
    fmt = '--', color = '#4682B4', label = None)                        #Fit
plt.errorbar(
    np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']),
                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)'])) ),
    np.log10(
   Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)']),
               np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['U(LR)'])),
    Lang_fits['U_1']['Q_max[mol/kg_be]'], Lang_fits['U_1']['K_L[1/M]']) ),
                fmt = '-.', color = '#4682B4', label = None)           #Lang Fit
        #note for the len I can use the df with std or without, since size is the same!
#### Pu
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Pu(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Pu(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Pu(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Pu(LR)'], 
            fmt = 'o', label = 'Pu', color = '#DC143C')
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)'])) ),
    1/Fre_fits['Pu_NL']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)']))
            )+ np.log10( Fre_fits['Pu_NL']['K[L^n/(kg*mol^{n-1})]']),
                fmt = '--', color = '#DC143C', label = None)        # Fre Fit
plt.errorbar(
    np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)']),
                np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)'])) ),
    np.log10(
   Lang_fit_eq(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)']),
               np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Pu(LR)'])),
    Lang_fits['Pu_1']['Q_max[mol/kg_be]'], Lang_fits['Pu_1']['K_L[1/M]']) ),
                fmt = '-.', color = '#DC143C', label = None)           #Lang Fit
#### Np
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Np(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Np(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Np(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Np(LR)'], 
            fmt = 'o', label = 'Np', color = '#228B22')
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Np(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Np(LR)'])) ),
    1/Fre_fits['Np']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Np(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Np(LR)']))
            )+ np.log10( Fre_fits['Np']['K[L^n/(kg*mol^{n-1})]']),
                fmt = '--', color = '#228B22', label = None) #Fit
#### Am
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Am(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Am(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Am(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Am(LR)'], 
            fmt = 'o', label = 'Am', color = '#FF8C00')
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)'])) ),
    1/Fre_fits['Am']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Am(LR)']))
            )+ np.log10( Fre_fits['Am']['K[L^n/(kg*mol^{n-1})]']),
                fmt = '--', color = '#FF8C00', label = None) #Fit
#### Cm
plt.errorbar(np.log10(Dict_el_sa_avg['< >']).iloc[:,3:].loc['Cm(LR)'], 
     np.log10(Dict_Qe['< >']).iloc[:,3:].loc['Cm(LR)'],
    yerr = np.abs(Dict_Qe['std'] / 
            (Dict_Qe['< >']*np.log(10))).iloc[:,3:].loc['Cm(LR)'], 
    xerr = np.abs(Dict_el_sa_avg['std']/
            (Dict_el_sa_avg['< >']*np.log(10))).iloc[:,3:].loc['Cm(LR)'], 
            fmt = 'o', label = 'Cm', color = '#9932CC')
plt.errorbar(np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)'])) ),
    1/Fre_fits['Cm']['n']*np.log10(np.linspace(np.min(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)']),
            np.max(Dict_el_sa_avg['< >'].iloc[:,3:].loc['Cm(LR)']))
            )+ np.log10( Fre_fits['Cm']['K[L^n/(kg*mol^{n-1})]']),
                fmt = '--', color = '#9932CC', label = None) #Fit
plt.ylabel('log($Q_e$ [mol/kg$_{ben}$])', fontsize= Font)              #ylabel
plt.xlabel('log($C_e$ [M])', fontsize = Font )
#plt.yscale('log')
#plt.xscale('log')
plt.legend(fontsize = Font)
plt.tick_params(axis='both', labelsize= Font) 
plt.minorticks_on()             #enabling minor grid lines
plt.grid(which = 'minor', linestyle=':', linewidth=0.5)        
plt.grid(which = 'major')
plt.savefig('Ads_iso_FreLang_fit_acti_no123.png', format='png', bbox_inches='tight')
plt.show()  
"""
Well, nice idea to keep it more compact and create less figures, but it is not
so clear, possibly because I am plotting lots of elements together.

To be further improved, if needed....
"""

#-----------------------------------------
#%%             D-R fitting model
#------------------------------------
"""
In order to see if I can get some info on the sorption mechanism, I will try
these model, for all the elements of interest!
"""

DR_U = Read_JRC.D_R_fit(
    Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['U(LR)'], 
    Dict_Qe['< >'].iloc[:,3:].loc['U(LR)'], 
            delta_Ce=Dict_el_sa_avg['std'].drop(S_rem, 
                                    axis =1).loc['U(LR)'], 
            delta_Qe = Dict_Qe['std'].iloc[:,3:].loc['U(LR)'], 
            folder_name = 'DR_Fits', Title = 'D-R fit U', 
            save_name = 'DR_fit_U')

DR_Np = Read_JRC.D_R_fit(
    Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Np(LR)'], 
    Dict_Qe['< >'].iloc[:,3:].loc['Np(LR)'], 
            delta_Ce=Dict_el_sa_avg['std'].drop(S_rem, 
                                    axis =1).loc['Np(LR)'], 
            delta_Qe = Dict_Qe['std'].iloc[:,3:].loc['Np(LR)'], 
            folder_name = 'DR_Fits', Title = 'D-R fit Np', 
            save_name = 'DR_fit_Np')

DR_Pu = Read_JRC.D_R_fit(
    Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Pu(LR)'], 
    Dict_Qe['< >'].iloc[:,3:].loc['Pu(LR)'], 
            delta_Ce=Dict_el_sa_avg['std'].drop(S_rem, 
                                    axis =1).loc['Pu(LR)'], 
            delta_Qe = Dict_Qe['std'].iloc[:,3:].loc['Pu(LR)'], 
            folder_name = 'DR_Fits', Title = 'D-R fit Pu', 
            save_name = 'DR_fit_Pu')

DR_Am = Read_JRC.D_R_fit(
    Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Am(LR)'], 
    Dict_Qe['< >'].iloc[:,3:].loc['Am(LR)'], 
            delta_Ce=Dict_el_sa_avg['std'].drop(S_rem, 
                                    axis =1).loc['Am(LR)'], 
            delta_Qe = Dict_Qe['std'].iloc[:,3:].loc['Am(LR)'], 
            folder_name = 'DR_Fits', Title = 'D-R fit Am', 
            save_name = 'DR_fit_Am')
DR_Cm = Read_JRC.D_R_fit(
    Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Cm(LR)'], 
    Dict_Qe['< >'].iloc[:,3:].loc['Cm(LR)'], 
            delta_Ce=Dict_el_sa_avg['std'].drop(S_rem, 
                                    axis =1).loc['Cm(LR)'], 
            delta_Qe = Dict_Qe['std'].iloc[:,3:].loc['Cm(LR)'], 
            folder_name = 'DR_Fits', Title = 'D-R fit Cm', 
            save_name = 'DR_fit_Cm')
DR_Cs = Read_JRC.D_R_fit(
    Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['Cs(LR)'], 
    Dict_Qe['< >'].iloc[:,3:].loc['Cs(LR)'], 
            delta_Ce=Dict_el_sa_avg['std'].drop(S_rem, 
                                    axis =1).loc['Cs(LR)'], 
            delta_Qe = Dict_Qe['std'].iloc[:,3:].loc['Cs(LR)'], 
            folder_name = 'DR_Fits', Title = 'D-R fit Cs', 
            save_name = 'DR_fit_Cs')

DR_La = Read_JRC.D_R_fit(
    Dict_el_sa_avg['< >'].drop(S_rem, axis =1).loc['La(LR)'], 
    Dict_Qe['< >'].iloc[:,3:].loc['La(LR)'], 
            delta_Ce=Dict_el_sa_avg['std'].drop(S_rem, 
                                    axis =1).loc['La(LR)'], 
            delta_Qe = Dict_Qe['std'].iloc[:,3:].loc['La(LR)'], 
            folder_name = 'DR_Fits', Title = 'D-R fit La', 
            save_name = 'DR_fit_La')

DR_fits= pd.DataFrame({'U': DR_U, 'Np': DR_Np,'Pu': DR_Pu,'Am': DR_Am,
                       'Cm': DR_Cm,'Cs': DR_Cs,'La': DR_La,} ) 
                         #gathering all together!
del (DR_U,DR_Np,DR_Pu,DR_Am,DR_Cm,DR_Cs,DR_La) 
"""
Almost all with Ion exchange/soft chemisorption

Note that Cs and Np have not so different F, this could be revealing
"""