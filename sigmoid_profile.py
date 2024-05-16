# -*- coding: utf-8 -*-
"""
Created on Wed May 15 16:15:59 2024

@author: avste
"""
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pickle


nrows=567 #number of data points in each vector
cellsize=700 #m

def make_profile(x, cellsize,slope, width_fac, height_fac, position_fac, vert_shift):
    L = np.amax(x)
    pos = position_fac * L
    width = width_fac * L
    y = (x - pos) / width 
    height = height_fac * L
    z = 0.5 * height * (1.0 - np.tanh(y))
    z += slope * (L - x)
    z += vert_shift*cellsize
    #z=z/cellsize
    return z

x = np.arange(0,nrows)*cellsize # m
z = make_profile(x, cellsize, 0, 0.15, .003, 0,0) #m
z2 = make_profile(x, cellsize, 0, 0.4, .005, 0,0) #m

fig,ax=plt.subplots(1,figsize=(5,5))
ax.set_title('Sigmoid Profile')
plt.plot(np.linspace(0,700,100),1000*np.ones(100),'--b',label='Ice Stability Line')
plt.plot(x*.001,z,'-m',label='Sigmoid 1') #Sigmoid S=0,WF=.15,HF=.003
plt.plot(x*.001,z2,'-k',label='Sigmoid 2') #Sigmoid S=0,WF=.4,HF=.005
plt.xlim([0,400])
plt.xlabel('Distance Along Watershed [km]')
plt.ylabel('Elevation [m]')
#plt.ylim([mini,maxi]) 
plt.legend()
plt.show()