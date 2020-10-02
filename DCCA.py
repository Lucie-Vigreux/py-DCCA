import glob
import os
import json
import pandas as pd
import numpy as np
import time as ti
import matplotlib.pyplot as plt
import math



 #Calculation of the DCCA coefficient
 
 #n+1 between len(x)/4 and len(x)/3 and greater than 10
 #cuting in sub-sequences overlaying each other
 #calculation of the polynomials associated to each sub-sequence
 
 #Inputs :
 #x : time serie
 #y : time serie compared to x
 #n (optionnal input, n=lenght(x)*7//24-1 by default) : length(x) - n is the number of sub-sequences
 
 #Output : 
 #Cdcca : Detrended Cross-Correlation Analysis Coefficient
 
def Dcca(x,y,n=0):
    L=len(x)
    if np.std(x)!=0 or np.mean(x)!=0:
         x=standardize_signal_vector(x) #standardization of signals 
    if np.std(y)!=0 or np.mean(y)!=0:
         y=standardize_signal_vector(y) 
    if n==0:
        n=L*7//24-1
    X=np.array(subsequences(serie_int(x),n)) #cuting in sub-sequences
    Y=np.array(subsequences(serie_int(y),n))
    Rx=np.array([Polynome(x) for x in X])
    Ry=np.array([Polynome(y) for y in Y])
    X_detrend=X-Rx #detrending
    Y_detrend=Y-Ry
    F2dcca=Fxy2(X_detrend,Y_detrend,L,n) #calculation of the covariance
    fx=Fxy2(X_detrend,X_detrend,L,n) #calculation of the variance of the x serie
    fy=Fxy2(Y_detrend,Y_detrend,L,n)
    Fdfa_x=fx**0.5 #standard deviation
    Fdfa_y=fy**0.5
    Cdcca=F2dcca/(Fdfa_x*Fdfa_y)
    return Cdcca



#Calculation of the integrated serie

#Input :
#s: time serie

#Output :
#RES : integrated serie of s

def serie_int(s):
    s=np.array(s)
    RES=[s[0]]
    for i in range(1,len(s)):
        RES.append(RES[i-1]+s[i])
    return RES



#Cutting of overlaying sub-sequences procedure

#Input : 
#s : time serie
#n + 1 : number of points in each overlaying sub-sequence

#Output :
#RES : list of length(s) - n overlaying sub-sequences from s

def subsequences(s,n):
    RES=[]
    L=len(s)
    k=(L//(n+1))
    for i in range(k):
        for j in range(n+1):
            sequence=s[i*(n+1)+j:i*(n+1)+j+n+1]
            if len(sequence)==n+1:
                RES.append(sequence)
    return RES


#Polynomization procedure

#Input :
#s : time serie

#Output : 
#Poly : Polynomial associated with the signal 

def Polynome(s):
    Poly=np.poly1d(np.polyfit(np.linspace(0,len(s)-1,len(s)),s,1),variable='z')(np.linspace(0,len(s)-1,len(s)))
    return Poly



def Fxy2(x,y,L,n):
    cov=[]
    k=[]
    for v in range(L-n):
       k.append(covariance(x[v],y[v]))
    K=np.mean(k)
    return K


def covariance(x,y):
    return sum(x*y)/(len(x)-1)



#Standardization procedure

#Input:
#signal : time serie

#Output:
#signal_stand : standardize signal

def standardize_signal_vector(signal):
    try:
        signal_stand = []
        signal_mean = np.mean(signal)
        signal_std = np.std(signal)
        for i, x in enumerate(signal):
            signal_stand.append(float((x - signal_mean) / signal_std))
        return signal_stand
    except Exception as e:
        print(e)
        return float('nan')




# Coarse-grained procedure

#Inputs :
# signal : original signal
# scaleFactor : scale factor (tau)

#Output :
#y : coarse-grained time series

def coarseGrained(signal, scaleFactor):
    N = len(signal)
    y_size = int(N / scaleFactor)
    y = []
    for j in range(0, y_size):
        y_i=0
        for i in range(j * scaleFactor, (j + 1) * scaleFactor):
            y_i += signal[i]
        y_i /= scaleFactor
        y.append(y_i)
    return y



# Time shift procedure

#Inputs :
# signal : original signal
# scaleFactor : scale factor (tau)

# Output :
#y : time shifted signals (# scaleFactor)

def timeShift(signal, scaleFactor):
    N = len(signal)
    y = []
    for k in range(0, scaleFactor):
        y_lim = k + math.floor((N-1 - k) / scaleFactor) * scaleFactor
        y_temp = []
        for i in range(k, y_lim+1, scaleFactor):
            y_temp.append(signal[i])
        y.append(y_temp)
    return y



