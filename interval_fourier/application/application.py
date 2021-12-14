
"""
    Created June 2021
    
    Author:    Marco Behrendt
               Leibniz Universität Hannover, Germany
               University of Liverpool, United Kingdom

    https://github.com/marcobehrendt/Projecting-interval-uncertainty-through-the-discrete-Fourier-transform
    
"""

import numpy
from numpy import (arange, cos, exp, linspace, mean, pi,  sin, zeros) # for convenience
from matplotlib import pyplot, cm
from scipy.spatial import ConvexHull

from fourier.number.number import Interval, IntervalVector 

# The code in this file should comply to PEP-8: https://realpython.com/python-pep8/
    
def subplots(figsize=(16,8),size=None): # wrapper of the matplotlib.pyplot figure gen
    if size is None:
        fig, ax  = pyplot.subplots(figsize=figsize)
    else:
        fig, ax = pyplot.subplots(figsize=figsize,size=size)
    return fig,ax

def plot_signal(signal,figsize=(18,6),xlabel=r'#$x$',ylabel=r'$x$',color=None,lw=1,title=None,ax=None,label=None):
    x = list(range(len(signal)))
    y = signal
    if ax is None:
        fig = pyplot.figure(figsize=figsize)
        ax = fig.subplots()
        ax.grid()
    ax.plot(x,y,marker='.',color=color,lw=lw,label=label)  # https://matplotlib.org/3.1.0/gallery/color/named_colors.html
    ax.set_xlabel(xlabel,fontsize=20)
    ax.set_ylabel(ylabel,fontsize=20)
    ax.tick_params(direction='out', length=6, width=2, labelsize=14)
    if title is not None:
        ax.set_title(title,fontsize=20)
    return None
  
def jonswap_spectrum(w,alpha,w_p,gamma,sigma1,sigma2):
    g = 9.81
    N = len(w)
    spectrum = numpy.zeros(N)
    r = numpy.zeros(N)
    for x in range(len(w)):
        if w[x] == 0:
            spectrum[x] = 0
        else:
            if w[x] <= w_p:          
                r[x] = exp(-(w[x]-w_p)**2/(2*sigma1**2*w_p**2))
            else:           
                r[x] = exp(-(w[x]-w_p)**2/(2*sigma2**2*w_p**2))   
            spectrum[x] = alpha * g**2 / w[x]**5 * exp( -5/4 * (w_p/w[x])**4 ) * gamma**r[x]
    return spectrum
    
def stochastic_process(spectrum, w, t):
    Nt = len(t)
    Nw = len(w)
    dw = w[1] - w[0]
    signal = numpy.zeros(Nt)
    for w_n in range(Nw):
        if w[w_n] == 0:
            A = 0
        else:
            A = (2*spectrum[w_n]*dw)**0.5
        phi = 2*pi*numpy.random.random_sample()
        signal += 2**0.5 * A * cos(w[w_n] * t + phi)
    return signal      
        
def wind_turbine(R,r,h_pile,rho_steel,c,k):
    A_steel = (R**2 - r**2)*pi
    V_steel = A_steel * h_pile
    m_steel = rho_steel * V_steel
    w0 = (k/m_steel)**0.5
    xi = c/(w0*2*m_steel)
    return w0,xi
    
def frequency_response_interval(w,spectrum,w0,xi):
    ai_low=[ai.lo() for ai in spectrum]
    ai_high=[ai.hi() for ai in spectrum]
    H_low = ai_low * abs(1/(w0**2 - w**2 + 2 * xi * w0*w*1j))**2
    H_high = ai_high * abs(1/(w0**2 - w**2 + 2 * xi * w0*w*1j))**2
    return H_low,H_high
    
def frequency_response(w,spectrum,w0,xi):
    H = spectrum * abs(1/(w0**2 - w**2 + 2 * xi * w0*w*1j))**2
    return H 

def periodogram(spectrum,t, dt):
    for x in range(len(spectrum)):    
        spectrum[x] = spectrum[x]**2*dt**2/t[len(t)-1]/(2*pi)
    return spectrum 
   
def plot_line(x,y,figsize=(18,6),xlabel=r'#$x$',ylabel='$x$',color=None,lw=1,title=None,ax=None,label=None):
    if ax is None:
        fig = pyplot.figure(figsize=figsize)
        ax = fig.subplots()
        ax.grid()
    ax.plot(x,y,marker='.',color=color,lw=lw,label=label)  # https://matplotlib.org/3.1.0/gallery/color/named_colors.html
    ax.set_xlabel(xlabel,fontsize=20)
    ax.set_ylabel(ylabel,fontsize=20)
    ax.tick_params(direction='out', length=6, width=2, labelsize=14)
    if title is not None:
        ax.set_title(title,fontsize=20)  
    return ax
    
def plot_bounds(x,bounds,color=None,alpha=None,ax=None):
    if ax is None:
        fig = pyplot.figure(figsize=(18,6))
        ax = fig.subplots()
        ax.grid()
    bounds_low=[ai.lo() for ai in bounds]
    bounds_high=[ai.hi() for ai in bounds]
    ax.fill_between(x,bounds_low,bounds_high,alpha=alpha,label='Interval',edgecolor='blue',lw=2,color=color)  
