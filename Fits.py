#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 09:09:07 2021

@author: Daniel López Aires// danlopair@gmail.com
"""


####### 0) General packages useful#############33
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize
from statsmodels.formula.api import ols 

Font = 18               #Fontsize, for the plots (labels, ticks, legends, etc)   
Markersize = 7

#%% ###### 1) Linear fit function ############################
####################################################################

def LinearRegression(x, y, delta_x =0, delta_y =0, npo = 100, x_label = 'x', y_label = 'y', 
                     Title = 'Linear fit',  Color = 'b', 
                     x_legend = 'x', y_legend = 'y', save_name = ''):
    '''
    Function that makes a linear regression of the 2 list (numpy preferred) X, Y
    and returns, if the fit equation is y = a x + b:
        a, \Delta{a}, b, \Delta{b}, r (correlation coefficient, aka r^2)
        

    *Inputs:
        .x, y = 1D numpy arrays containing the data to fit/ pd Series, with maching indexes!
        .delta_x/y = 0: error of y and x, in case you want to plot the fit with it.
        .npo = 'number of points of the linspace for the fit plotting. Default value = 100
        .x_label, y_label= x and y label, for the plot. Default value: 'x' and 'y'
        .x_legend, y_legend: y and y label, to appear in the legend. Default: 'x' and 'y'
        .Title = title of the plot. Default: 'Linear fit'
        .save_name = filename of the plot, if it wants to be save. Default value = '' ==> no saving.
            this variable is followed by .png for savinf
        .Color = 'b': color for the plot
        
    *Outputs:
        .df Series with fit parameters: slope and its error, intercept and its error, correlation
        	coefficient
	'''
    

    ################ 1.1) Fit ##################
    
    N = len(x)                                  #vector length
    data = pd.DataFrame({'X': x, 'Y': y})  
                                        #data in pandas style to do the fit
    fit= ols("Y ~ X", data).fit()                                #fit
                    #Important to write "Y ~ X". If you write "X ~ Y", it will do
                    #the opposite analysis
                    
    a = fit.params[1]                         #slope
    b = fit.params[0]                     #intercept with the X axis
    r = fit.rsquared                             #correlation coefficient R
    delta_a = fit.bse[1]                   #std(a) = delta(a)
    delta_b = fit.bse[0]                    
    
    '''
    Note those error values are very similar to the excel values (19/9/23), the ones
    from my formulas were different, we keep with this!
    '''
    #print(fit.summary())
    
    
    ################## 2) Error calculation ##############
    #This cformulas come from 1st course physics, so I better use the ones from the model xD
    #delta_aa = 3 * np.sqrt(a**2/(N-2) * (1/r**2 - 1))
    #delta_bb = 3 * np.sqrt( sum(x**2) * delta_aa**2 / N)
  
    
    ################ 3) Storing #########################
    values = {'a' : a, '\Delta(a)' : delta_a,
              'b' : b, '\Delta(b)' : delta_b, 'r' : r }
    Ser_values = pd.Series(values, name = Title)      #gathering output in a df Series
            #naming the column like the post_title variable, since this variable is an isotope: U238
    
    
    ############# 4) Plot of the fit##########
    x_vector = np.linspace(min(x),max(x),npo)         #for the fit plotting
    
    fig = plt.figure(figsize=(11,8))  #width, heigh 6.4*4.8 inches by default
    ax = fig.add_subplot(111)
    ax.errorbar(x, y, delta_y, delta_x, 'o', color = Color, markersize = Markersize, label = 'Data')
    ax.plot(x_vector, linear(x_vector, a, b),'--', color = Color,
            label= 'Fit: ' + y_legend + f' = {a:.1e} ' + x_legend + f'+{b:.1e}' + ',\n r= ' + f'{r:.5f}')      #fit
            #.2f to show 2 decimals on the coefficients!
            #2e for scientific notation with 2 significative digits
    ax.set_title(Title, fontsize=22)          #title
    ax.set_xlabel(x_label, fontsize = Font)                                    #xlabel
    ax.set_ylabel(y_label, fontsize= Font)                                    #ylabel
    ax.tick_params(axis='both', labelsize= Font)            #size of tick labels  
    ax.grid(True)                                              #show grid
    ax.legend(fontsize = Font)             #legend
                    #Plot of the fit equation. (0,0) is lower-left corner, and (1,1) the upper right
    plt.savefig(save_name +'.png', format='png', bbox_inches='tight')                
                    ###This require some thoughts!!!!! to automatize the show of the equation!!!!!!!!!!!
    
    
    ################# 5) Return of values############
        #the values will be returned in a dictionary indicating what is each
        #value

    return Ser_values




#####Fit function    
def linear(x, a, b):       #Definition of the function to use to fit the data
    '''
    Linear function:
        f(x) = a*x + b
    '''
	
    return a * x + b 




'''Long way to calculate the sums (matlab influenced)
    suma_x2 = 0                                     #inicializaicon
    suma_y2 = 0                                     #inicializaicon
    suma_xy = 0                                     #inicializaicon              
    
        """
        Way to calculate the sum x^2 and sum y^2, but not sum x*y (Matlab style)
    for element in x:               #calc of sum x^2
        suma_x2 = suma_x2 + element**2
    
    for element in y:                #calc of sum y^2
        suma_y2 = suma_y2 + element**2
        """

    for value in range(0,N):
        suma_x2 = suma_x2 + x[value]**2
        suma_y2 = suma_y2 + y[value]**2
        suma_xy = suma_xy + x[value]*y[value]
'''
    

#%%  #### 2) Quadratic fit function ############################
####################################################################


def QuadraticRegression(x, y, npo = 100):
    '''
    Function that makes a linear regression of the 2 list (numpy preferred) X, Y
    and returns, if the fit equation is y = ax**2 + bx + c:
        a,b,c, \Delta{a,b,c}, n, r (correlation coefficient, aka r^2)
        
    *Inputs:
        .x, y = 1D numpy arrays containing the data to fit / pd series (same index names)
        .npo = 'number of points of the linspace for the fit plotting. Default value = 100
        
    *Outputs:
        .pd Series with the fit parameters and their errors, correlation
        	coefficient
        '''

    
    ######## 2.1) Fit ############
    N = len(x)                                  #vector length


	#initial = [max(y_data), x_data[0], (x_data[1] - x_data[0]) * 5]
                #initial guesses for the fit. If None, this does not work, so this
                #is very important when having an offset! Thank you 
                #Lucas Hermann Negri (PeakUtils)
                
    cuadratic_fit = scipy.optimize.curve_fit(cuadratic, x, y)#, initial)
                

    opt_values = cuadratic_fit[0]   #optimal values of the function to fit the data
    cov_of_opt_val = cuadratic_fit[1]            #covariances of the optimal values
    #the diagonal are the variance of the parameter to estimate.
    
    a = opt_values[0]  
    b = opt_values[1]
    c = opt_values[2]

    perr = np.sqrt(np.diag(cov_of_opt_val))        #standard deviation error (el 
                                                #error de toa la via vamos)
                                                
    delta_a = perr[0]                                       #\Delta{a}
    delta_b = perr[1]                                       #\Delta{b}
    delta_c = perr[2]                                       #\Delta{c}


    #find r-squared of polynomial model with degree = 3
    r = r_square(x, y, 2)


    #### 2.2) Storing###
    values = {'a' : a, 'b' : b, 'c' : c,
              'r' : r, '\Delta{a}' : delta_a,  '\Delta{b}' : delta_b, 
              '\Delta{c}' : delta_c}
    
    Ser_values = pd.Series(values)      #gathering output in a df Series
    
    
    ####3) Plot of the fit####
    x_vector = np.linspace(min(x),max(x),npo)         #for the fit plotting
    plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
    plt.plot(x, y, 'bo', markersize = Markersize, label = 'Data' )            #original data
    plt.plot(x_vector, cuadratic(x_vector, a, b, c), 'r--', linewidth=2, 
             label=f'Fit: y = {a:.1e}x^2 + {b:.1e}x + {c:.1e}' + ', r= ' + f'{r:.5f}')      #fit
                        #Like that I put the fit eq into the legend plot ;)
    plt.title('Quadratic fit', fontsize=22)          #title
    plt.xlabel("X ", fontsize= Font )                                    #xlabel
    plt.ylabel('Y', fontsize= Font)                                    #ylabel
    plt.tick_params(axis='both', labelsize= Font)            #size of tick labels  
    plt.grid(True)                                              #show grid
    plt.legend(fontsize= Font)             #legend
    #plt.text(2.7,300, 'y(x) = {0:1.3f}x^2 + {1:1.3f}x + {2:1.3f} ; r = {3:1.3f}'
       #  .format(a, b, c,r_cuadratic_fit), fontsize=14) #10 default size
       #plt.xlim(0,15)
       #plt.savefig('sigma2_vs_peaknumber_cuadratic_fit_py.png', format='png')
   
        
   
    ######## 2.3) Return of values ###########
        #the values will be returned in a dictionary indicating what is each
        #value

    return Ser_values


#####Fit function    

def cuadratic(x, a, b, c):       
    '''
    Definition of the function to use to fit the data:
        f(x) = a x**2 + b*x + c
    '''
    return a * x**2 + b*x + c 


def r_square(x, y, degree):
    '''
    Define function to calculate r-squared 
    (https://www.statology.org/quadratic-regression-python/)
    '''
    coeffs = np.polyfit(x, y, degree)
    p = np.poly1d(coeffs)
    
    #calculate r-squared
    yhat = p(x)
    ybar = np.sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)
    sstot = np.sum((y - ybar)**2)
    results = ssreg / sstot

    return results


#%%  ###### 3) Gaussian fit function ############################
####################################################################
    
def Gaussian_fit(x, y, index_df = 0, N = 100):
    """GAUSSIAN FIT

    This script contains a function that fits a set of data to a gaussian function,
    giving the statistical parameters and plotting the fit. The gaussian function is:
    
    Gauss(x) = Heigh * np.exp(- (x-Mean)**2 / (2 * Std_dev**2)),
    
        Heigh = amplitude of the gaussian function
        Mean = mean value of the gaussian
        Std_dev = standard deviation of the gaussian

    *Inputs:
        .x, y = 1D numpy arrays containing the data to fit
        .index_df = index name of the dataframe containing the outputs
        .N = 'number of points of the linspace for the fit plotting. Default value = 100
        
    *Outputs:
        .fit parameters: heigh, amplitude and standar deviation and their error
        .FHWM and its error


    @Watch out! 
	i)Python indexes begin at 0, while usually the indexes start at 1,
		that is, if the line in your .txt is Z, for python you should type Z-1
	ii) Sometimes computer error occurs and Std_dev is <0, although the fit is good.
		This is just a computer error, ignore it. This issue can not be solved
		by fitting with the variance = std_dev^2 for the moment, in the future I
		may think how to fix this :)

    """
    ########## 3.1) Fit calc ################3
    
    
    initial = [max(y), x[0], (x[1] - x[0]) * 5 ]
                #initial guesses for the fit. If None, this does not work, so this
                #is very important when having an offset! Thank you 
                #Lucas Hermann Negri (PeakUtils) 
                
    fit = scipy.optimize.curve_fit(gaussian, x, y, initial)       #fit


    #Obtaining the values from the fit:
    opt_values = fit[0]             #optimal values of the function to fit the data
    cov_of_opt_val = fit[1]                     #covariances of the optimal values
                    #the diagonal are the variance of the parameter to estimate.   
    
    heigh = opt_values[0]                       		#heigh of the fit
    mean = opt_values[1]                        	#mean value of the fit 
    sigma = opt_values[2]                       #variance of the fit STD DEV!!
    #sigma = np.sqrt(opt_values[2]) 			#std deviation of the fit if using variance as parameter
    
    perr = np.sqrt(np.diag(cov_of_opt_val))        #standard deviation error ('el 
                                 #error de toa la via vamos')
    
    Delta_mean = perr[1]                #error of the mean
    Delta_heigh = perr[0]               #error of the heigh
    #Delta_sigma = 1 / (2 * sigma) * perr[2]              #error of the standar deviation if using variance as parameter   
    Delta_sigma = perr[2] 					#error of the standar deviation
  
  #source: 
  #https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html


    #Now lets compute other useful parameters:                  

    FWHM = 2 * np.sqrt(2 * np.log(2)) * sigma                   #FWHM of the peak
    Delta_FWHM = 2 * np.sqrt(2 * np.log(2)) * Delta_sigma     #error of the FWHM
    #print('FWHM: ' + str(FWHM) + ' +/- ' + str(Delta_FWHM) + ' MeV')

    Res = 100 * FWHM / mean                               #Resolution [%]
    Delta_Res = Res * np.sqrt( (Delta_FWHM / FWHM)**2 + (Delta_mean / mean)**2 )
                                                    #Error of the resolution [%]


    ########## 3.2) Plot of the fit################3
    x_vector = np.linspace(min(x),max(x),N)         #for the fit plotting
    
    plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
    plt.plot(x, y, 'bo', label = 'data', markersize = Markersize)     #original data
    plt.plot(x_vector, gaussian(x_vector, heigh, mean, sigma), 'r--', 
      label= 'Fit: y = ' + f' = {heigh:.2f} * exp('  + f' -(x-{mean:.2f})^2 / 2*{sigma:.2f}^2)' + '\nRes= ' + f'{Res:.2f}') #fit
             #.2f to show 2 decimals on the coefficients!
             #2e for scientific notation with 2 significative digits)           #fit
    plt.legend(fontsize = 12)
    plt.title('Gaussian fit of the data', fontsize=22)                      #title
    plt.xlabel("X", fontsize= Font)                                    #xlabel
    plt.ylabel("Y", fontsize= Font)                                    #ylabel
    plt.tick_params(axis='both', labelsize= Font)                  #size of tick labels  
    plt.grid(True)                                              #show grid
    #plt.xlim(5.35,5.55)                                         #limits of x axis
    
    
    
   #3.3 ) Return of values########################
   #the values will be returned in a DataFrame indicating what is each
   #value
    aux = {'heigh' : heigh, '\Delta(heigh)' : Delta_heigh, 
              'mean' : mean, '\Delta(mean)' : Delta_mean,  
              'sigma' : sigma, '\Delta(sigma)' : Delta_sigma, 
              'FWHM' : FWHM, '\Delta(FWHM)' : Delta_FWHM,
              'Res[%]' : Res, '\Delta(Res[%])' : Delta_Res
              }     #variable containing everyting
    
    #Let´s print that so it appears in the command line:
    print('\n#######################\n')
    print('Gaussian fit parameters ########\n')
    print(aux)
    print('\n ###############')
    
    values = pd.DataFrame(aux, index = [index_df])  #dataframe creation
    return values



######Fit function    
def gaussian(x, Heigh, Mean, Std_dev):
	return Heigh * np.exp(- (x-Mean)**2 / (2 * Std_dev**2)) 


