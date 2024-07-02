''' 
This is a module for post-processing of numerical and eCerimental t histories, 
primarily to calculate lag t
'''

import numpy as np
import scipy.constants as sc
from scipy.optimize import fsolve, curve_fit
from scipy import interpolate


def lag_time(t, C, frac=0.5):
    
    '''
    t: time, x-axis, 1D array
    C: Concentration, y-axis, 1D array with same length as t
    frac: fraction of t (from the end, not the beginning) to use for line fit
    ''' 
    
    Nfit = int(len(t) * frac)
    t_fit = t[-Nfit:]
    y_fit = C[-Nfit:]
    
    def fitfun(x, m, b): return m * x + b
    
    popt, pcov = curve_fit(fitfun, t_fit, y_fit) 
    m,b = popt 
    lag = -b/m
    
    perr = np.sqrt(np.diag(pcov))
    lag_error = np.sqrt( (perr[0] * (b/m**2))**2  +  (perr[1] * (-1/m))**2 )
    
    return dict(lag_time=lag, lag_error=lag_error, popt=popt, fiterr = perr)

