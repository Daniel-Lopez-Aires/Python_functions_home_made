#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 09:54:19 2021

@author: dla
"""

#####-2) Modules needed#######
import scipy.optimize              #to do the fit. doing only import scipy sometimes
                                #gives an error, so have to do this
import numpy as np
import matplotlib.pyplot as plt  #for simplicity, to not write matplotlib.pyplot
        #everytime we want to plot something    
        
    

def QuadraticRegression(x, y):
    '''
    Function that makes a linear regression of the 2 list (numpy preferred) X, Y
    and returns, if the fit equation is y = m x + n:
        m, \Delta{m}, n, \Delta{n}, r (correlation coefficient, aka r^2)
        '''
    

    
    
    #####-1) Debug, values of x and y to test it####
    #x = [1, 2]
    #y = [3, 4]
    

    #####0) Preliminary work####
    x = np.array(x)                         #conversion to np array, in case the list
                                    #is not an numpy array
    y = np.array(y)    
    N = len(x)                                  #vector length
    
    #1) Fit


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




    
    ####2) Error calculation###

    ####3) Storing###
    values = {'a' : a, 'b' : b, 'c' : c,
              'r' : r, '\Delta{a}' : delta_a,  '\Delta{b}' : delta_b, 
              '\Delta{c}' : delta_c}
    
    ####4) Plot of the fit###

    plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
    plt.plot(x, y, 'r*', linewidth=3 )                         #original data
    plt.plot(x, cuadratic(x, a, b, c), linewidth=3)      #fit

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
      
    ####5) Return of values####
        #the values will be returned in a dictionary indicating what is each
        #value

    return values


#####Fit function    

def cuadratic(x, a, b, c):       
    '''Definition of the function to use to fit the data
    '''
    return a * x**2 + b*x + c 


def r_square(x, y, degree):
    '''#define function to calculate r-squared 
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

    
