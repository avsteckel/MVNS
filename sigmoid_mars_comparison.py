# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 09:57:35 2022

@author: avste
"""
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pickle

# Run importARCdata first to load the river valley networks from ArcGIS data. 
# Averaged elevation (avgC) and distance (y) is saved into pickle file
#C:\Users\avste\Local_Documents\Geology\MRVN_Code\python_scripts\data
#C:\Users\Geology\Desktop\python_scripts\data
local_path=r'C:\Users\avste\local_doc\Mars\data'
with open(local_path+r'\mars_river_valleys.pkl', 'rb') as f:
    test = pickle.load(f)
coords=test[0]
y=test[1]
avgC=test[2]
n=np.size(coords) #number of river valley networks
color = cm.rainbow(np.linspace(0, 1, n))
nrows=2000 #number of data points in each vector
cellsize=200 #m
maxi=3500
mini=-2000

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


plt.figure(1)
#matplotlib.rcParams.update({'font.size': 22})
x = np.arange(0,nrows)*cellsize # m
z = make_profile(x, cellsize, 0, 0.6, .003, 0.5,-7) #m
z2 = make_profile(x, cellsize, 0, 0.4, .005, 0.5,-1.5) #m
z3 = make_profile(x, cellsize, 0, 0.5, .005, 0.5,4.5) #m

fig,ax=plt.subplots(1,figsize=(10,10))
ax.set_title('Sigmoid Profile')
for i in range(n):
    plt.plot(y[i],avgC[i],c=color[i],label=coords[i])
plt.plot(np.linspace(0,700,100),1000*np.ones(100),'--b',label='Ice Stability Line')
#plt.plot(x*.001,z,'-r',zorder=2,label='Calculated')
plt.plot(x*.001,z[::-1],'-m',label='Sigmoid S=0,WF=.15,HF=.003')
plt.plot(x*.001,z2[::-1],'-k',label='Sigmoid S=0,WF=.4,HF=.005')
plt.plot(x*.001,z3,'-g',label='Sigmoid S=0,WF=.4,HF=.005')
#plt.plot(x[0:50]*.001,1480*np.ones(nrows)[0:50])
#plt.plot(x[150:nrows]*.001,580*np.ones(nrows)[150:nrows])
#plt.plot(x[0:50]*.001,905*np.ones(nrows)[0:50])
#plt.plot(x[150:nrows]*.001,5*np.ones(nrows)[150:nrows])
plt.xlabel('Distance Along Watershed [km]')
plt.ylabel('Elevation [m]')
plt.ylim([mini,maxi]) 
plt.legend()
plt.show()


plot2=plt.figure(2)
fig,ax=plt.subplots(1,figsize=(10,10))
ax.set_title('Sigmoid Profile')
plt.xlim(700, 0)
for i in range(n):
    if i==3:
        plt.plot(y[i],avgC[i][::-1],c=color[i],label=coords[i])
        print(i)
    else:
        plt.plot(y[i],avgC[i],c=color[i],label=coords[i])
        print(i)
plt.plot(np.linspace(0,700,100),1000*np.ones(100),'--b',label='Ice Stability Line')
#plt.plot(x*.001,z,'-r',zorder=2,label='Calculated')
plt.plot(x*.001,z[::-1],'-m',label='Sigmoid WF=.15,HF=.003')
plt.plot(x*.001,z2[::-1],'-k',label='Sigmoid WF=.4,HF=.005')
plt.plot(x*.001,z3[::-1],'-k',label='Sigmoid WF=.4,HF=.005')
#plt.plot(x[0:50]*.001,1480*np.ones(nrows)[0:50])
#plt.plot(x[150:nrows]*.001,580*np.ones(nrows)[150:nrows])
#plt.plot(x[0:50]*.001,905*np.ones(nrows)[0:50])
#plt.plot(x[150:nrows]*.001,5*np.ones(nrows)[150:nrows])
plt.xlabel('Distance Along Watershed [km]')
plt.ylabel('Elevation [m]')
plt.ylim([mini,maxi]) 
#plt.xlim(max(x*.001), min(x*.001))
plt.legend()
plt.show()

plot3=plt.figure(3)
z4 = make_profile(x, cellsize, 0, 0.5, .005, 0.5,2) #m
fig,ax=plt.subplots(1,figsize=(10,10))
ax.set_title('Sigmoid Profile')
plt.xlim(700, 0)
for i in range(n):
    if i==3:
        plt.plot(y[i],avgC[i][::-1],c=color[i],label=coords[i],marker='o')
        print(i)
    else:
        plt.plot(y[i],avgC[i],c=color[i],label=coords[i],marker='o')
        print(i)
plt.plot(np.linspace(0,700,100),1000*np.ones(100),'--b',label='Ice Stability Line')
#plt.plot(x*.001,z,'-r',zorder=2,label='Calculated')
plt.plot(x*.001,z[::-1],'-m',label='Sigmoid WF=.15,HF=.003')
plt.plot(x*.001,z4[::-1],'-k',label='Sigmoid WF=.4,HF=.005')
#plt.plot(x[0:50]*.001,1480*np.ones(nrows)[0:50])
#plt.plot(x[150:nrows]*.001,580*np.ones(nrows)[150:nrows])
#plt.plot(x[0:50]*.001,905*np.ones(nrows)[0:50])
#plt.plot(x[150:nrows]*.001,5*np.ones(nrows)[150:nrows])
plt.xlabel('Distance Along Watershed [km]')
plt.ylabel('Elevation [m]')
plt.ylim([mini,maxi]) 
#plt.xlim(max(x*.001), min(x*.001))
plt.legend()
plt.show()



### Create a best fit for sigmoid 1
x1=np.concatenate([y[0],y[4],y[5],y[6]],axis=0)
y1=np.concatenate([avgC[0],avgC[4],avgC[5],avgC[6]],axis=0)
p1=np.polyfit(x1, y1, 2)

x2=np.concatenate([y[1],y[2],y[3]],axis=0)
y2=np.concatenate([avgC[1],avgC[2],avgC[3][::-1]],axis=0)
p2=np.polyfit(x2, y2, 2)

plot3=plt.figure(3)
fig,ax=plt.subplots(1,figsize=(10,10))
ax.set_title('Best Fit Slope')
plt.xlim(700, 0)
plt.plot(x1,y1,c='k',ls='',marker='.')
plt.plot(x2,y2,c='g',ls='',marker='.')
plt.plot(np.linspace(0,700,100),1000*np.ones(100),'--b',label='Ice Stability Line')
plt.plot(np.unique(x1), np.poly1d(np.polyfit(x1, y1, 2))(np.unique(x1)),c='r')
plt.plot(np.unique(x2), np.poly1d(np.polyfit(x2, y2, 2))(np.unique(x2)),c='c')
#plt.plot(x*.001,z,'-r',zorder=2,label='Calculated')
#plt.plot(x*.001,z[::-1],'-m',label='Sigmoid WF=.15,HF=.003',lw=5)
#plt.plot(x*.001,z4[::-1],'-k',label='Sigmoid WF=.4,HF=.005')
plt.xlabel('Distance Along Watershed [km]')
plt.ylabel('Elevation [m]')
plt.ylim([mini,maxi]) 
plt.show()