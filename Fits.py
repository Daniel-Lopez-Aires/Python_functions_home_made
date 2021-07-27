#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 09:09:07 2021

@author: Daniel López Aires// danlopair@gmail.com
"""


#######0) General packages useful#############33
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize
from statsmodels.formula.api import ols 

#%%  ###############################################################
############### 1) Linear fit function ############################
####################################################################

def LinearRegression(x, y, npo = 100):
    '''
    Function that makes a linear regression of the 2 list (numpy preferred) X, Y
    and returns, if the fit equation is y = m x + n:
        m, \Delta{m}, n, \Delta{n}, r (correlation coefficient, aka r^2)
        

    *Inputs:
        .x, y = 1D numpy arrays containing the data to fit
        .npo = 'number of points of the linspace for the fit plotting. Default value = 100
        
    *Outputs:
        .fit parameters: slope and its error, intercept and its error, correlation
        	coefficient
	'''
    

    ################ 1.1) Fit ##################
    
    N = len(x)                                  #vector length
    data = pd.DataFrame({'X': x, 'Y': y})  
                                        #data in pandas style to do the fit
    fit= ols("Y ~ X", data).fit()                                #fit
                    #Important to write "Y ~ X". If you write "X ~ Y", it will do
                    #the opposite analysis
                    
    slope = fit.params[1]                         #slope
    intercept = fit.params[0]                     #intercept with the X axis
    r = fit.rsquared                             #correlation coefficient R
    
    ################## 2) Error calculation ##############

    
    delta_slope = 3 * np.sqrt(slope**2/(N-2) * (1/r**2 - 1))
    delta_intercept = 3 * np.sqrt( sum(x**2) * delta_slope**2 / N)
  
    ####3) Storing###
    values = {'Slope' : slope, 'Intercept' : intercept, 'r' : r, 
              '\Delta{slope}' : delta_slope, 
              '\Delta{intercept}' : delta_intercept}
    
    
    
    #############3) Plot of the fit##########
    x_vector = np.linspace(min(x),max(x),npo)         #for the fit plotting
    
    plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
    plt.plot(x,y, 'r*')
    plt.plot(x_vector, linear(x_vector, values['Slope'], values['Intercept']) )      #fit
    plt.title('Linear regression', fontsize=22)          #title
    plt.xlabel("x ", fontsize=14)                                    #xlabel
    plt.ylabel('y', fontsize=14)                                    #ylabel
    plt.tick_params(axis='both', labelsize=14)            #size of tick labels  
    plt.grid(True)                                              #show grid
    plt.legend(['data','linear fit'], fontsize=14)             #legend
     
    
    
    
    #################4) Return of values############
        #the values will be returned in a dictionary indicating what is each
        #value

    return values




#####Fit function    
def linear(x, m, n):       #Definition of the function to use to fit the data
	return m * x + n 




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
    

#%%  ###############################################################
############### 2) Quadratic fit function ############################
####################################################################


def QuadraticRegression(x, y, npo = 100):
    '''
    Function that makes a linear regression of the 2 list (numpy preferred) X, Y
    and returns, if the fit equation is y = m x + n:
        m, \Delta{m}, n, \Delta{n}, r (correlation coefficient, aka r^2)
        
    *Inputs:
        .x, y = 1D numpy arrays containing the data to fit
        .npo = 'number of points of the linspace for the fit plotting. Default value = 100
        
    *Outputs:
        .fit parameters and their errors, correlation
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
    
    
    
    
    ####3) Plot of the fit####
    x_vector = np.linspace(min(x),max(x),npo)         #for the fit plotting
    plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
    plt.plot(x, y, 'r*', linewidth=3 )                         #original data
    plt.plot(x_vector, cuadratic(x_vector, a, b, c), linewidth=3)      #fit

    plt.title('Quadratic fit', fontsize=22)          #title
    plt.xlabel("X ", fontsize=14)                                    #xlabel
    plt.ylabel('Y', fontsize=14)                                    #ylabel
    plt.tick_params(axis='both', labelsize=14)            #size of tick labels  
    plt.grid(True)                                              #show grid
    plt.legend(['data', 'quadratic fit'], fontsize=16)             #legend
    #plt.text(2.7,300, 'y(x) = {0:1.3f}x^2 + {1:1.3f}x + {2:1.3f} ; r = {3:1.3f}'
       #  .format(a, b, c,r_cuadratic_fit), fontsize=14) #10 default size
       #plt.xlim(0,15)
       #plt.savefig('sigma2_vs_peaknumber_cuadratic_fit_py.png', format='png')
   
        
   
    ######## 2.3) Return of values ###########
        #the values will be returned in a dictionary indicating what is each
        #value

    return values


#####Fit function    

def cuadratic(x, a, b, c):       
    '''
    Definition of the function to use to fit the data
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


#%%  ###############################################################
############### 3) Gaussian fit function ############################
####################################################################
    
def Gaussian_fit(x,y,N=100):
    """GAUSSIAN FIT

    This script contains a function that fits a set of data to a gaussian function,
    giving the statistical parameters and plotting the fit. The gaussian function is:
    
    Gauss(x) = Heigh * np.exp(- (x-Mean)**2 / (2 * Std_dev**2)),
    
        Heigh = amplitude of the gaussian function
        Mean = mean value of the gaussian
        Std_dev = standard deviation of the gaussian

    *Inputs:
        .x, y = 1D numpy arrays containing the data to fit
        .npo = 'number of points of the linspace for the fit plotting. Default value = 100
        
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

    R = 100 * FWHM / mean                               #Resolution [%]
    Delta_R = R * np.sqrt( (Delta_FWHM / FWHM)**2 + (Delta_mean / mean)**2 )
                                                    #Error of the resolution [%]


    ########## 3.2) Plot of the fit################3
    x_vector = np.linspace(min(x),max(x),N)         #for the fit plotting
    
    plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
    plt.plot(x, y, 'b.')                         #original data
    plt.plot(x_vector, gaussian(x_vector, heigh, mean, sigma), 'r--')           #fit
    plt.title('Gaussian fit of the data', fontsize=20)                      #title
    #plt.xlabel("E (MeV)", fontsize=10)                                    #xlabel
    #plt.ylabel("Cuentas", fontsize=10)                                    #ylabel
    plt.legend(['data', 'gaussian fit'], fontsize=10) 
    plt.tick_params(axis='both', labelsize=10)                  #size of tick labels  
    plt.grid(True)                                              #show grid
    #plt.xlim(5.35,5.55)                                         #limits of x axis
    
    
    
   #3.3 ) Return of values########################
   #the values will be returned in a dictionary indicating what is each
   #value
    values = {'heigh' : heigh, '\Delta(heigh)' : Delta_heigh, 
              'mean' : mean, '\Delta(mean)' : Delta_mean,  
              'sigma' : sigma, '\Delta(sigma)' : Delta_sigma, 
              'FWHM' : FWHM, '\Delta(FWHM)' : Delta_FWHM,
              'R[%]' : R, '\Delta(R[%])' : Delta_R
              }
    return values



######Fit function    
def gaussian(x, Heigh, Mean, Std_dev):
	return Heigh * np.exp(- (x-Mean)**2 / (2 * Std_dev**2)) 


