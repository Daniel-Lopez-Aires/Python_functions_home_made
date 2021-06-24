#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 09:54:19 2021

@author: dla
"""

#####-2) Modules needed#######
import pandas as pd   
from statsmodels.formula.api import ols 
import numpy as np
import matplotlib.pyplot as plt  #for simplicity, to not write matplotlib.pyplot
        #everytime we want to plot something    
        
    

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
    

    
    
    #####-1) Debug, values of x and y to test it####
    #x = [1, 2]
    #y = [3, 4]
    

    #####0) Preliminary work####

    N = len(x)                                  #vector length
    
    #1) Fit
    data = pd.DataFrame({'X': x, 'Y': y})  
                                        #data in panda style to do the fit
    fit= ols("Y ~ X", data).fit()                                #fit
                    #Important to write "Y ~ X". If you write "X ~ Y", it will do
                    #the opposite analysis
                    
    slope = fit.params[1]                         #slope
    intercept = fit.params[0]                     #intercept with the X axis
    r = fit.rsquared                             #correlation coefficient R
    
    ####2) Error calculation###

    
    delta_slope = 3 * np.sqrt(slope**2/(N-2) * (1/r**2 - 1))
    delta_intercept = 3 * np.sqrt( sum(x**2) * delta_slope**2 / N)
  
    ####3) Storing###
    values = {'Slope' : slope, 'Intercept' : intercept, 'r' : r, 
              '\Delta{slope}' : delta_slope, 
              '\Delta{intercept}' : delta_intercept}
    
    ####4) Plot of the fit###
    x_vector = np.linspace(min(x),max(x),npo)         #for the fit plotting
    
    plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
    plt.plot(x,y, 'r*')
    plt.plot(x, linear(x, values['Slope'], values['Intercept']) )      #fit
    plt.title('Linear regression', fontsize=22)          #title
    plt.xlabel("x ", fontsize=14)                                    #xlabel
    plt.ylabel('y', fontsize=14)                                    #ylabel
    plt.tick_params(axis='both', labelsize=14)            #size of tick labels  
    plt.grid(True)                                              #show grid
    plt.legend(['data','linear fit'], fontsize=14)             #legend
      
    ####5) Return of values####
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
    
