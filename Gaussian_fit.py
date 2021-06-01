#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 14:19:15 2021

@author: dla

GAUSSIAN FIT

This script contains a function that fit a set of data to a gaussian function,
giving the statistical parameters and plotting the fit. The gaussian function is:
    
    Gauss(x) = Heigh * np.exp(- (x-Mean)**2 / (2 * Std_dev**2)),
    
        Heigh = amplitude of the gaussian function
        Mean = mean value of the gaussian
        Std_dev = standard deviation of the gaussian

Inputs: 1D array (numbers), x and y. 

Whatch out, because python indexes begin at 0, while normaly the indexes start at 1,
#s that is the line in your .txt is Z, for python you should type Z-1

"""

def Gaussian_fit(x,y):

    
    #######0) General packages useful###########

    import matplotlib.pyplot as plt  #for simplicity, to not write matplotlib.pyplot
        #everytime we want to plot something
    import scipy.optimize              #to do the fit. doing only import scipy sometimes
                                #gives an error, so have to do this
    import numpy as np          #np contain linspaces as np.linspace(a,b,N)
    ####
    


    #######1) Fit ###########
     
     #Fit function:
     
    def gaussian(x, Heigh, Mean, Std_dev):
        return Heigh * np.exp(- (x-Mean)**2 / (2 * Std_dev**2)) 
     
    
    #Data:
    x_data = np.array(x)    
    y_data = np.array(y)

    #Fit:
    initial = [max(y_data), x_data[0], (x_data[1] - x_data[0]) * 5 ]
                #initial guesses for the fit. If None, this does not work, so this
                #is very important when having an offset! Thank you 
                #Lucas Hermann Negri (PeakUtils) 	#std ved!!!!
                
    fit = scipy.optimize.curve_fit(gaussian, x_data, y_data, initial)       #fit


    #Obtaining the values from the fit:
    opt_values = fit[0]             #optimal values of the function to fit the data
    cov_of_opt_val = fit[1]                     #covariances of the optimal values
                    #the diagonal are the variance of the parameter to estimate.   
    
    heigh = opt_values[0]                       #heigh of the fit
    mean = opt_values[1]                        #mean value of the fit 
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



    ########## 2)Plot of the fit################3
    plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
    plt.plot(x_data, y_data, label = 'data')                         #original data
    plt.plot(x_data, gaussian(x_data, heigh, mean, sigma), 'ro')           #fit
    plt.title('Gaussian fit of the data', fontsize=20)                      #title
    #plt.xlabel("E (MeV)", fontsize=10)                                    #xlabel
    #plt.ylabel("Cuentas", fontsize=10)                                    #ylabel
    plt.legend(['data', 'gaussian fit'], fontsize=10) 
    plt.tick_params(axis='both', labelsize=10)                  #size of tick labels  
    plt.grid(True)                                              #show grid
    #plt.xlim(5.35,5.55)                                         #limits of x axis
    
    
    
   #3) Return of values########################
   #the values will be returned in a dictionary indicating what is each
   #value
    values = {'heigh' : heigh, '\Delta(heigh)' : Delta_heigh, 
              'mean' : mean, '\Delta(mean)' : Delta_mean,  
              'sigma' : sigma, '\Delta(sigma)' : Delta_sigma, 
              'FWHM' : FWHM, '\Delta(FWHM)' : Delta_FWHM}
    return values
    
   #4) Comments
   #sometimes, the std_dev returned is negative, although the fit is perfect. This is just
   #a computer error, and sadly do not know how to avoid it. If trying to use the variance as
   #the fitting parameter, the fit does not work, so for the moment I wont do that, simply will
   #leave this as an open question.
