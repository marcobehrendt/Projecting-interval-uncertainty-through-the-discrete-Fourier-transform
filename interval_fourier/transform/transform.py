
"""
    Created September 2019
    Author: Marco De Angelis
    University of Liverpool
    github.com/marcodeangelis
    GNU LGPL v3.0

    Code for the propagation of intervals through the discrete Fourier transform. 

    This code lives at: https://github.com/marcodeangelis/Fourier-transform

    If you use this code in your work/research and want to cite it, you can cite  
    the research paper at: https://arxiv.org/abs/2012.09778
"""

import numpy
from numpy import (arange, cos, exp, linspace, mean, pi,  sin, zeros) # for convenience
from matplotlib import pyplot, cm
from scipy.spatial import ConvexHull

from fourier.number.number import Interval, IntervalVector 

# The code in this file should comply to PEP-8: https://realpython.com/python-pep8/

def intervalize(signal, plusminus=0.): # intervalizes the signal. Default returns an interval vector with 0 widths.
    if is_iterable(plusminus): # plusminus is iterable
        intervalsignal = [(s-u, s+u) for s,u in zip(signal,plusminus)]
    else:
        u = plusminus
        intervalsignal = [(s-u, s+u) for s in signal]
    return IntervalVector(intervalsignal)

def basic_stationary_spectrum(NN,wu,d=0.6,nf=12,scale=1):    # this code was ported from Matlab
    N = NN/2 #summation times
    dw = wu/N
    w = arange(0,wu-dw,dw) #w=(0:dw:wu-dw)
    fr = w
    S = (nf**4+(4*((d*nf*fr)**2)))/(((nf**2)-(fr**2))**2+(2*d*fr*nf)**2)/scale
    S[0] = 0 
    return S, dw # returns two numpy arrays
def gen_stat_sig(spectrum, time, freq, noise=0.): # this code was ported from Matlab
    Fn = len(spectrum)
    N = len(time)
    delta_f = freq[1] - freq[0]
    signal = numpy.zeros(N)
    for i in range(Fn):
        u = numpy.random.random_sample() # random number in (0,1]
        signal += ((2*spectrum[i]*delta_f)**0.5)*(cos(freq[i]*time+2*pi*u))    # convert 1-sided power spectrum Y from power to amplitude  (i.e. sqrt(2*power)=amplitude)
    if noise != 0.:
        noise_arr = noise * numpy.random.random_sample(N)
        noise_arr = noise_arr - numpy.mean(noise)
        signal = signal + noise_arr
    return signal  # returns a numpy array
def generate_signal(N = 6): 
    NN = 2**N # NN = 2**6 #64
    t = linspace(0, 2*pi-2*pi/NN, NN)
    w = linspace(0, NN//2-1, NN//2)
    Sx, _ =  basic_stationary_spectrum(NN, NN//2-1, 0.6, 12, 1) #% this generates a power spectrum
    f = gen_stat_sig(Sx,t,w) #% this draws a signal based on the spectrum (this is our assumed recording)
    return f,w
def generate_intervalized_signal(N=6, plusminus=0.):
    NN = 2**N # NN = 2**6 #64
    t = linspace(0, 2*pi-2*pi/NN, NN)
    w = linspace(0, NN//2-1, NN//2)
    Sx, _ =  basic_stationary_spectrum(NN, NN//2-1, 0.6, 12, 1) #% this generates a power spectrum
    f = gen_stat_sig(Sx,t,w) #% this draws a signal based on the spectrum (this is our assumed recording)
    ff = intervalize(f,plusminus=plusminus)
    return ff,w

def Fourier_transform(signal): # precise version of the periodogram, this is the DFT of the signal
    F=[]
    N = len(signal)
    for w_i in range(N//2):
        f = 0
        for n in range(N):
            theta_n = (-2 * pi * n * w_i) / N
            f += signal[n] * exp(1j*theta_n) # 
        F.append(f)
    return F # always outputs a Python list

def Fourier_amplitude(signal): 
    F=[]
    N = len(signal) # Needs to be a power of 2
    for w_i in range(N//2):
        f = 0
        for n in range(N):
            theta_n = (-2 * pi * n * w_i) / N
            f += signal[n] * exp(1j*theta_n) 
        F.append(abs(f))
    return F  # always outputs a Python list

def compute_amplitude_bounds_givenfrequency(intsignal, frequency): 
    def in_complex_hull(chull,point=(0,0)): # default check the origin of the complex plain
            chull_RI=[(c.real,c.imag) for c in chull]
            ch_add = chull_RI+[point] 
            chull_add=ConvexHull(ch_add)
            cchh_v = [ch_add[k][0]+1j*ch_add[k][1] for k in chull_add.vertices]
            return set(cchh_v) == set(chull)
    N = len(intsignal)
    pair = [exp(-2*pi*1j*0*(frequency)/N) * intsignal[0].lo(), exp(-2*pi*1j*0*(frequency)/N) * intsignal[0].hi()]
    convhull = pair
    ci=exp(-2*pi*1j*0*(frequency)/N) * intsignal[0]
    for k in range(1,N):
        pair =  [exp(-2*pi*1j*k*(frequency)/N) * intsignal[k].lo(),exp(-2*pi*1j*k*(frequency)/N) * intsignal[k].hi()]
        ci += exp(-2*pi*1j*k*(frequency)/N) * intsignal[k]
        convhull = [[pair[0] + ch for ch in convhull],[pair[1] + ch for ch in convhull]]
        convhull = convhull[0] + convhull[1]
        pairs_RI = [[p.real, p.imag] for p in convhull]
        hull = ConvexHull(pairs_RI)
        convhull = [pairs_RI[k][0] + 1j*pairs_RI[k][1] for k in hull.vertices] # turns CH into a list of complex values
    ch_max = numpy.argmax([abs(h) for h in convhull])
    ch_min = numpy.argmin([abs(h) for h in convhull])
    ci_ll = abs(ci.real().lo() +1j*ci.imag().lo())
    ci_lh = abs(ci.real().lo() +1j*ci.imag().hi())
    ci_hl = abs(ci.real().hi() +1j*ci.imag().lo())
    ci_hh = abs(ci.real().hi() +1j*ci.imag().hi())
    bI_hi = max(ci_ll,ci_lh,ci_hl,ci_hh)
    if all(ci.stradzero()):
        bI_lo = 0
    elif ci.stradzero()[0]: # real component straddle zero
        bI_lo = min(abs(ci.imag().lo()),abs(ci.imag().hi()))
    elif ci.stradzero()[1]: # imag component straddle zero
        bI_lo = min(abs(ci.real().lo()),abs(ci.real().hi()))
    else:
        bI_lo = min(ci_ll,ci_lh,ci_hl,ci_hh)
    boundI = Interval(bI_lo,bI_hi)
    if in_complex_hull(convhull):
        boundC = Interval(0,abs(convhull[ch_max]))
    else:
        boundC = Interval(abs(convhull[ch_min]),abs(convhull[ch_max]))
    return boundI, boundC, (convhull[ch_min], convhull[ch_max])

def compute_amplitude_bounds(intervalsignal): 
    N = len(intervalsignal)
    BOUNDS_I = []
    BOUNDS_C = []
    for frequency in range(1,N//2):
        bI,bC,_ = compute_amplitude_bounds_givenfrequency(intervalsignal,frequency)
        BOUNDS_I.append(bI)
        BOUNDS_C.append(bC)
    return BOUNDS_I, BOUNDS_C

def compute_amplitude_bounds_bruteforce(intervalsignal, frequency):
    pass

def compute_amplitude_bounds_selective(intervalsignal, frequency):
    pass

def is_iterable(x):
    tt = ['float','int']
    m=[]
    for t in tt:
         m+=[f"'{t}' object is not iterable"]
    iterable = True
    try:
        iter(x)
    except TypeError as e:
        if str(e) in m:
            iterable = False
    return iterable

def subplotting(N):
    sqroot = N**.5
    if abs(int(sqroot)-sqroot)<1e-9: # N is a square
        n,m = int(sqroot), int(sqroot)
    else:
        n = int(sqroot) # cast the square root into the its corresponding integer
        m = n+1
        while m*n < N:
            m+=1
    return n,m

def final_box_convexhull(intsignal, freq): # port of Liam's Matlab code
    N = len(intsignal)
    pair = [exp(-2*pi*1j*0*(freq)/N) * intsignal[0].lo(), exp(-2*pi*1j*0*(freq)/N) * intsignal[0].hi()]
    convhull = pair
    ci=exp(-2*pi*1j*0*(freq)/N) * intsignal[0]
    for k in range(1,N):
        pair =  [exp(-2*pi*1j*k*(freq)/N) * intsignal[k].lo(),exp(-2*pi*1j*k*(freq)/N) * intsignal[k].hi()]
        ci += exp(-2*pi*1j*k*(freq)/N) * intsignal[k]
        convhull = [[pair[0] + ch for ch in convhull],[pair[1] + ch for ch in convhull]]
        convhull = convhull[0] + convhull[1]
        pairs_RI = [[p.real, p.imag] for p in convhull]
        hull = ConvexHull(pairs_RI)
        convhull = [pairs_RI[k][0] + 1j*pairs_RI[k][1] for k in hull.vertices]
    return ci, convhull

def verify_selective_with_plot(intervalsignal: IntervalVector,freq,figsize=(21,12),save=None,aspect=None,sharexy=None):
    def in_complex_hull(chull,point=(0,0)): # default check the origin of the complex plain
        chull_RI=[(c.real,c.imag) for c in chull]
        ch_add = chull_RI+[point] 
        chull_add=ConvexHull(ch_add)
        cchh_v = [ch_add[k][0]+1j*ch_add[k][1] for k in chull_add.vertices]
        return len(cchh_v) == len(chull)
    ww_iterable = is_iterable(freq)
    if ww_iterable:
        n,m = subplotting(len(freq))
        if sharexy is not None:
            fig, ax = pyplot.subplots(n, m,figsize=figsize,sharey=sharexy[1],sharex=sharexy[0])
        else:
            fig, ax = pyplot.subplots(n, m,figsize=figsize)
        ax = ax.flatten()
        for i,w in enumerate(freq):
            ci,ch = final_box_convexhull(intervalsignal,w)
            ch_max = numpy.argmax([abs(h) for h in ch])
            ch_min = numpy.argmin([abs(h) for h in ch])
            if in_complex_hull(ch):
                bC_lo = 0+1j*0
                bC_hi = ch[ch_max]
            else:
                bC_lo = ch[ch_min]
                bC_hi = ch[ch_max]
            xx = [c.real for c in ch]
            yy = [c.imag for c in ch]
            xx = xx + [xx[0]]
            yy = yy + [yy[0]]
            ax[i].plot(xx,yy, color='red',lw=2,label='Covex hull')
            ax[i].scatter(xx,yy,marker='.')
            ax[i].plot([min(0,min(xx)),max(max(xx),0)],[0,0],color='blue',lw=1) # x axis
            ax[i].plot([0,0],[min(0,min(yy)),max(max(yy),0)],color='blue',lw=1,label='Coordinate axes') # y axis
            ax[i].scatter(bC_lo.real,bC_lo.imag,marker='v',s=100,label='Min ampl.')
            ax[i].scatter(bC_hi.real,bC_hi.imag,marker='^',s=100,label='Max ampl.')
            ax[i].plot([0,bC_hi.real],[0,bC_hi.imag],lw=1)
            ax[i].plot([0,bC_lo.real],[0,bC_lo.imag],lw=1)
            r1,r2,i1,i2 = ci.lo().real,ci.hi().real,ci.lo().imag,ci.hi().imag
            ax[i].scatter([r1,r1,r2,r2],[i1,i2,i1,i2],marker='.')
            ax[i].fill_between(x=[r1,r2],y1=[i1,i1],y2=[i2,i2],alpha=0.4, color='gray',label='Bounding box')
            p1 = float('%g'%bC_hi.real)
            p2 = float('%g'%bC_hi.imag)
            ax[i].text(bC_hi.real,bC_hi.imag,f'max amplitude \n ({p1}, {p2}i)',fontdict=None,fontsize=16)
            p3 = float('%.4g'%bC_lo.real)
            p4 = float('%.4g'%bC_lo.imag)
            ax[i].text(bC_lo.real,bC_lo.imag,f'min amplitude \n ({p3}, {p4}i)',fontdict=None,fontsize=16)
            ax[i].set_xlabel(r'real($z_{k='+str(w)+'}$)',fontsize=20)
            ax[i].set_ylabel(r'imag($z_{k='+str(w)+'}$)',fontsize=20)
            if aspect is not None:
                ax[i].set_aspect(aspect) # 'equal', 'auto'
            ax[i].tick_params(direction='out', length=6, width=2, labelsize=14)
            ax[i].grid()
        ax[2].legend(fontsize=14,loc='upper left') # https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.legend.html
    else:
        w=freq
        fig = pyplot.figure(figsize=figsize)
        ax = fig.subplots()
        ci,ch = final_box_convexhull(intervalsignal,w)
        ch_max = numpy.argmax([abs(h) for h in ch])
        ch_min = numpy.argmin([abs(h) for h in ch])
        if in_complex_hull(ch):
            bC_lo = 0+1j*0
            bC_hi = ch[ch_max]
        else:
            bC_lo = ch[ch_min]
            bC_hi = ch[ch_max]
        xx = [c.real for c in ch]
        yy = [c.imag for c in ch]
        xx = xx + [xx[0]]
        yy = yy + [yy[0]]
        ax.plot(xx,yy, color='red',lw=2,label='Covex hull')
        ax.scatter(xx,yy,marker='.')
        ax.plot([min(0,min(xx)),max(max(xx),0)],[0,0],color='blue',lw=1) # x axis
        ax.plot([0,0],[min(0,min(yy)),max(max(yy),0)],color='blue',lw=1,label='Coordinate axes') # y axis
        ax.scatter(bC_lo.real,bC_lo.imag,marker='v',s=100,label='Min ampl.')
        ax.scatter(bC_hi.real,bC_hi.imag,marker='^',s=100,label='Max ampl.')
        ax.plot([0,bC_hi.real],[0,bC_hi.imag],lw=1)
        ax.plot([0,bC_lo.real],[0,bC_lo.imag],lw=1)
        r1,r2,i1,i2 = ci.lo().real,ci.hi().real,ci.lo().imag,ci.hi().imag
        ax.scatter([r1,r1,r2,r2],[i1,i2,i1,i2],marker='.')
        ax.fill_between(x=[r1,r2],y1=[i1,i1],y2=[i2,i2],alpha=0.4, color='gray',label='Bounding box')
        p1 = float('%g'%bC_hi.real)
        p2 = float('%g'%bC_hi.imag)
        ax.text(bC_hi.real,bC_hi.imag,f'max amplitude \n ({p1}, {p2}i)',fontdict=None,fontsize=16)
        p3 = float('%.4g'%bC_lo.real)
        p4 = float('%.4g'%bC_lo.imag)
        ax.text(bC_lo.real,bC_lo.imag,f'min amplitude \n ({p3}, {p4}i)',fontdict=None,fontsize=16)
        ax.set_xlabel(r'real($z_{k='+str(w)+'}$)',fontsize=20)
        ax.set_ylabel(r'imag($z_{k='+str(w)+'}$)',fontsize=20)
        if aspect is not None:
            ax.set_aspect(aspect) # 'equal', 'auto'
        ax.tick_params(direction='out', length=6, width=2, labelsize=14)
        ax.grid()
        ax.legend(fontsize=14,loc='upper right') # https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.legend.html
    fig.tight_layout()
    if save is not None:
        pyplot.savefig(save)
    pyplot.show()



# ---------------------------- #
# ----- Plots start here ----- #

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

def plot_y(y,figsize=(18,6),xlabel=r'#$x$',ylabel=r'$x$',color=None,lw=1,title=None,ax=None,label=None):
    x = list(range(len(y)))
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

def plot_xy(x,y,figsize=(18,6),xlabel=r'$x$',ylabel=r'$y$',color=None,lw=1,title=None,ax=None):
    pass

# ------ Plots end here ------ #
# ---------------------------- #
