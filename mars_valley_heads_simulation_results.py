# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 13:14:52 2022

@author: avste
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.colors
from matplotlib import cm
import pickle

# Location folders are saved on local computer
# C:\Users\Geology\OneDrive - UCB-O365\Desktop\ARC_Documents\Results
localPath=r'C:\Users\Geology\OneDrive - UCB-O365\Desktop\ARC_Documents'
localPath2=r'C:\Users\Geology\OneDrive - UCB-O365\Desktop\arc_rasters'

filepath1=[localPath+r'\Results\15n30e.txt',
           localPath+r'\Results\43_12.txt',
           localPath+r'\Results\0_23.txt',
           localPath+r'\Results\6s45e.txt',
           localPath+r'\Results\3s5e.txt',
           localPath+r'\Results\10s14w.txt',
           localPath+r'\Results\7s3e.txt']

filepath2=[localPath2+r'\mask_15n30e.txt',
           localPath2+r'\12n43e_mask.txt',
           localPath2+r'\0n23e_mask.txt',
           localPath2+r'\6s45e_mask.txt',
           localPath2+r'\mask_3s5e.txt',
           localPath2+r'\mask_10s14w.txt',
           localPath2+r'\7s3emask.txt']

filepath3=[localPath+r'\Results\tables\15n30e.txt',
           localPath+r'\Results\tables\12n43e.txt',
           localPath+r'\Results\tables\0n23e.txt',
           localPath+r'\Results\tables\6s45e.txt',
           localPath+r'\Results\tables\3s5e.txt',
           localPath+r'\Results\tables\10s14w.txt',
           localPath+r'\Results\tables\7e3s.txt',]

# dem_mask is all the valley networks. This file is too large to use
filepath4=localPath+r'\Results\tables\dem.txt' 

# table with 1km bins for  the region of interest only
filepath5=localPath+r'\Results\tables\kmbins_roi.txt' 

# Simulation Results
filepath6=r'C:\Users\Geology\OneDrive - UCB-O365\Desktop\Simulation_Results\ww_nets.txt'
filepath7=r'C:\Users\Geology\OneDrive - UCB-O365\Desktop\Simulation_Results\ww_ascii.asc'


coords=['15N_30E','12N_43E','0N_23E','6S_45E','3S_5E','10S_14W','7S_3E']

# select beginning and end point in ArcGIS of valley net in degree units
longARC=np.array([[29.4512,32.0174],[42.1289,44.3517],[22.4106,22.9043],[44.4593,44.0912],[2.61685,7.295844],[-13.6979,-13.6348],[-1.5379,8.417 ]])
latARC=np.array([[18.4392,10.7211],[12.9265,9.7546],[3.6185,-2.5061],[-3.5612, -8.8692],[-1.671917,-4.154649],[-6.5579,-14.3066],[-3.7659,-9.7707]])
topoARC=np.array([[941,-2241],[1109,-609],[2018,-237],[3408,730],[291,-2383],[-481,-2690],[615,-2441]])

n=np.size(coords)
#n=1
#Initialize variables
avgC=[]
distance=[]
vbininfo=[] # stop / end of each bin (elevations)
vcount=[] # number  of heads in each bin
vharea=[] #area of each bin
toggle=0
for i in range(0,7):
    
    x=[]
    y=[]
    # Local max and min elevation
    maxi=topoARC[i,0]
    mini=topoARC[i,1]
    
    long=longARC[i]
    lat=latARC[i]
    
    with open(filepath1[i], 'r') as file1:
        file_header=file1.readlines()[:6]
    
    # Split the header
    file_header=[item.strip().split()[-1] for item in file_header]
    ncols=int(file_header[0])
    nrows=int(file_header[1])
    xllcorner=float(file_header[2]) # keep units of deg
    yllcorner=float(file_header[3])
    cellsize=float(file_header[4])
    nodata=float(file_header[5])
    
    #(grid,z)=read_esri_ascii(filepath1[i], name='topographic__elevation',halo=0) # Halo adds perimiter of nodata values around DEM
    
    #print('Header information for file',coords[i])
    #print('Number of Columns:',ncols)
    #print('Number of Rows:', nrows)
    #print('Location of top left corner in Long, Lat:',xllcorner,yllcorner)
    
    
    # Read in the array
    z=np.loadtxt(filepath1[i],dtype=float,skiprows=6)
    z[z==nodata] = np.nan # Set nodata = nan for plotting
    extent_z=[xllcorner, xllcorner+ncols*cellsize,yllcorner,yllcorner+nrows*cellsize]
    
    #import csv
    data=pd.read_csv(filepath3[i], sep=',')   
    firstorder=data.loc[data['VNOrder']==1]
    x_cords=list(firstorder['START_X'])
    y_cords=list(firstorder['START_Y'])
    z_cords=list(firstorder['START_Z'])
    #    temp=csv.reader(headLoc, delimeter='\t')
    #    d=list(temp)
    
    n_bins=5
    area_conversion=200*200 # fix area conversion
    bininfo=np.histogram(z_cords,bins=n_bins)
    vbininfo.append(bininfo)
    count=np.histogram(z[~np.isnan(z)],bins=bininfo[1]) 
    vcount.append(count)
    harea=bininfo[0]/(count[0]*area_conversion) 
    vharea.append(harea)
    
    if toggle:
        #from landlab.components import StreamPowerEroder, FlowAccumulator, LinearDiffuser, DepressionFinderAndRouter, FastscapeEroder, DrainageDensity
        from landlab.io.esri_ascii import read_esri_ascii
        (mask,val)=read_esri_ascii(filepath2[i], name='valley__mask',halo=0) # Halo adds perimiter of nodata values around DEM
        #m=read_esri_ascii(filepath3, name='Mask',halo=0)
        ids=np.where(mask.at_node['valley__mask']==1)
        
        
        strtitle=coords[i]
        plt.figure()
        plt.gcf().set_size_inches(6, 6)
        plt.imshow(z,extent=extent_z,cmap='pink',vmin=mini,vmax=maxi)
        plt.plot(x_cords,y_cords,'*b')
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.plot(mask.x_of_node[ids],mask.y_of_node[ids],',')
        #plt.imshow(zMask,extent=extent_zMask,zorder=2,vmin=1,cmap='gray')
        cbar = plt.colorbar(label='Elevation [meters]')
        cbar.ax.tick_params(labelsize=20)
        for t in cbar.ax.get_yticklabels():
            t.set_fontsize(20)
        #plt.colorbar(label='Elevation [meters]',fontsize=20)
        #imshow_grid(mask,'valley__mask',vmin=0,vmax=1,cmap='gray_r')
        #imshow_grid(grid,'topographic__elevation',vmin=mini,vmax=maxi)
        plt.xlabel('Degrees (Longitude)',fontsize=20)
        plt.ylabel('Degrees (Latitude)',fontsize=20)
        #plt.ylim([np.min(yy),np.max(yy)]) 
        #plt.xlim([np.min(x),np.max(x)]) 
        plt.title(strtitle,fontsize=20)
        
    
        #strtitle2='Head Distributions '+coords[i]
        ## create elevation bins
        ##minbin=np.nanmin(z)
        ##maxbin=np.nanmax(z)
        ##bin_space=(maxbin-minbin)/n_bins
        ##binvec=np.arange(minbin, maxbin, bin_space)
        #plt.figure(2)
        #bininfo=plt.hist(z_cords,bins=n_bins)
        #plt.xlabel('Elevation (meters)')
        #plt.ylabel('Number of heads')
        #plt.title(strtitle2)
        #plt.show()
        
        #plt.figure()
        #plt.hist(count[1][:-1],count[1],weights=count[0]*area_conversion)
        #plt.xlabel('Elevation (meters)')
        #plt.ylabel('Area in bin (meters^2)')
        #plt.title('Area bins of landscape (number of pixels*pixel area)')
        #plt.show()
        
        #number of heads in each bin divided by area in each bin
        strtitle2='Normalized Head Distributions '+coords[i]
        #plt.figure()
        #plt.hist(count[1][:-1],bins=bininfo[1],weights=harea)
        #plt.ylim((0,np.max(harea)))
        #plt.xlabel('Elevation (meters)')
        #plt.ylabel('Number of Heads/Area [m^-2]')
        #plt.title(strtitle2)
        #plt.show()
        
plt.show()
data=pd.read_csv(filepath4, sep=',')   
firstorder=data.loc[data['VNOrder']==1]
x_global=list(firstorder['START_X'])
y_global=list(firstorder['START_Y'])
xz_global=list(firstorder['START_Z'])

# Import simulation results for Warm Wet
local_path=r'C:\Users\Geology\OneDrive - UCB-O365\Desktop\Simulation_Results'
ww_table=local_path+r'\ww_stream_table.txt'
ww_order=local_path+r'\ww_nets.txt'
ww_raster=local_path+r'\raster_ww_stream.txt'
ww_ascii=r'C:\Users\Geology\OneDrive - UCB-O365\Desktop\ARC_Documents\simulations\ww_s1_ascii.asc'
ww_ascii2=r'C:\Users\Geology\OneDrive - UCB-O365\Desktop\ARC_Documents\simulations\ww_s2_ascii.asc'
WWS1=r'C:\Users\Geology\OneDrive - UCB-O365\Desktop\Simulation_Results\Point_WW_S1_Table.csv'
WWS2=r'C:\Users\Geology\OneDrive - UCB-O365\Desktop\Simulation_Results\Point_WW_S2_Table.csv'

with open(ww_raster, 'r') as file1:
    file_header=file1.readlines()[:6]
file_header=[item.strip().split()[-1] for item in file_header]
ncols=int(file_header[0])
nrows=int(file_header[1])
xllcorner=int(file_header[2])
yllcorner=int(file_header[3])
cellsize=int(file_header[4])
nodata=int(file_header[5])
extent_ww=[xllcorner,xllcorner+ncols*cellsize,yllcorner,yllcorner+nrows*cellsize]

data=pd.read_csv(ww_table, sep=',')
firstorder=data.loc[data['grid_code']==1]
x_1=list(firstorder['x_1'])
y_1=list(firstorder['y_1'])
#z_cords=list(firstorder['START_Z'])

#fig=plt.figure(figsize=(20,10))
#z=np.loadtxt(ww_raster,skiprows=6)
#z[z==nodata]=np.nan
#plt.imshow(z,extent=extent_ww)
#plt.plot(x_1,y_1,'*b',markersize=2)
#plt.show()

with open(ww_ascii, 'r') as file1:
    file_header=file1.readlines()[:9]
elev=np.loadtxt(ww_ascii,skiprows=9)
y_cords=np.array(y_1)
x_cords=np.array(x_1)
x_cells=y_cords//cellsize
y_cells=x_cords//cellsize

x_cells.astype(int)
y_cells.astype(int)
l=np.shape(x_cells)
wwheads=[]
for i in range(0,l[0]):
     wwheads.append(elev[int(x_cells[i]),int(y_cells[i])])

with open(ww_ascii2, 'r') as file2:
    file_header2=file2.readlines()[:9]
elev2=np.loadtxt(ww_ascii2,skiprows=9)
wwheads2=[]
for i in range(0,l[0]):
     wwheads2.append(elev2[int(x_cells[i]),int(y_cells[i])])


nbins=10


# Excel data exported from ArcGIS
head_WWS1=pd.read_csv(WWS1)
head_WWS2=pd.read_csv(WWS2)

HWWS1=np.asarray(head_WWS1)
HWWS2=np.asarray(head_WWS2)

wwbininfo1=plt.hist(HWWS1,bins=nbins)
plt.title('WW Sigmoid 1 Head Elevations')
plt.xlabel('Elevation (m)')
plt.ylabel('Counts')
plt.show()

wwbininfo2=plt.hist(HWWS2,bins=nbins)
plt.title('WW Sigmoid 2 Head Elevations')
plt.xlabel('Elevation (m)')
plt.ylabel('Counts')
plt.show()

wwarea_conversion=cellsize*cellsize
wwcount1=np.histogram(elev,bins=wwbininfo1[1])
wwharea1=wwbininfo1[0]/(wwcount1[0]*wwarea_conversion)

wwcount2=np.histogram(elev2,bins=wwbininfo2[1])
wwharea2=wwbininfo2[0]/(wwcount2[0]*wwarea_conversion)


plt.hist(wwcount1[1][:-1],bins=wwbininfo1[1],weights=wwharea1)
plt.title('Warm-Wet Sim Normalized heads / area [m2]')
plt.xlabel('Elevation (m)')
plt.ylabel('Counts/area [m2]')
plt.show()

plt.hist(wwcount2[1][:-1],bins=wwbininfo2[1],weights=wwharea2)
plt.title('Warm-Wet Sim Normalized heads / area [m2]')
plt.xlabel('Elevation (m)')
plt.ylabel('Counts/area [m2]')
plt.show()

# Manually create Icy Cold Histogram
ICcounts=nrows # One "head" created manually at each row
ICheads=np.ones(nrows)*1000 # set elevation of input water to be at 1000m
ic_ascii=local_path+r'\ic_ascii.asc'
elev_ic=np.loadtxt(ic_ascii,skiprows=9)
ICbininfo=plt.hist(ICheads,bins=wwbininfo1[1])
plt.title('Icy-Cold Head Elevations')
plt.xlabel('Elevation (m)')
plt.ylabel('Counts')
plt.show()

iccount=np.histogram(elev_ic,bins=ICbininfo[1])
icharea=ICbininfo[0]/(iccount[0]*wwarea_conversion) # area conversion is same as ww
plt.hist(iccount[1][:-1],bins=ICbininfo[1],weights=icharea)
plt.title('Icy-Cold Sim Normalized heads / area [m2]')
plt.xlabel('Elevation (m)')
plt.ylabel('Counts/area [m2]')
plt.show()

#n_bins_global=40
#plt.figure(20)
#plt.hist(xz_global,bins=n_bins_global)
#plt.xlabel('Elevation (meters)')
#plt.ylabel('Number of heads')
#plt.title('Global Head Distributions')
#plt.show()

# 1km bins from ArcGIS
roi_data=pd.read_csv(filepath5, sep=';',decimal=',')   
roi_area=list(roi_data['Shape_Leng'])

# Subset region of interest (ROI)
roi=np.where((np.asarray(y_global)>-20) & (np.asarray(y_global)<20) & (np.asarray(x_global)>-30))
roi=roi[0]
xzcords_roi=[]
for i in range(np.size(roi)):
    xzcords_roi.append(xz_global[roi[i]])
    
with open('roi_file.pkl','rb') as f:
    hist_roi=pickle.load(f)



plt.figure(20)
roicount=plt.hist(xzcords_roi,bins=[-6000,-5000,-4000,-3000,-2000,-1000,0,1000,2000,3000,4000,5000])
plt.xlabel('Elevation (meters)')
plt.ylabel('Number of heads')
plt.title('Region of Interest Head Distributions')
plt.show()


plt.figure(21)
norm_roi=plt.hist(roicount[1][:-1],bins=roicount[1],weights=roicount[0]/roi_area)
plt.xlabel('Elevation (meters)')
plt.ylabel('Number of heads / area [m^-2]')
plt.xlim((-6000,5000))
plt.title('Normalized ROI Distributions')
plt.show()

#colors=[cm.jet(0), cm.jet, (0, .1, .2), (0, .1, .3), (0, .2, .4), (0, .3, .5), (.1, .3, .6)]

center_bins=[]
for i in range(0,7):
    temp=[]
    for j in range(0,n_bins):
        temp.append(np.average([vbininfo[i][1][j],vbininfo[i][1][j+1]]))
    center_bins.append(temp)
        
axis_break1 = 2.5*10**(-8)
axis_break2 = 3.25*10**(-8)

fig, ax1 = plt.subplots()
#plt.hist(roicount[1][:-1],bins=roicount[1],weights=roicount[0]/roi_area,color="w",label="ROI",edgecolor='black')
plt.yscale('log')
#plt.hist(hist_roi[1][:-1],bins=hist_roi[1],weights=hist_roi[0],color='w',label="ROI",edgecolor="black")
plt.ylabel('Number of heads / area [m^-2] in ROI')
#ax1.cla() # clear the axis 
ax2 = ax1.twinx()
for i in range(0,7):
    plt.hist(vcount[i][1][:-1],bins=vbininfo[i][1],weights=vharea[i],color=cm.jet(i*.15),alpha=0.2)
    plt.hist(vcount[i][1][:-1],bins=vbininfo[i][1],weights=vharea[i],label=coords[i],color=cm.jet(i*.15),histtype=u'step')
    plt.plot(center_bins[i],vharea[i],marker='*',color=cm.jet(i*.15),linestyle = 'None')
#plt.hist(wwcount[1][:-1],bins=wwbininfo[1],weights=wwharea,label='Warm-Wet')
#plt.hist(iccount[1][:-1],bins=ICbininfo[1],weights=icharea,label='Icy-Cold')  
plt.xlabel('Elevation (meters)')
plt.ylabel('Number of heads / area [m^-2] in Individual Valleys')
plt.xlim((-4000,5000))

plt.title('Normalized Head Distributions')
fig.legend(loc="upper left",bbox_to_anchor =(.12,.89))
plt.show()


fig, ax1 = plt.subplots()
#plt.hist(roicount[1][:-1],bins=roicount[1],weights=roicount[0]/roi_area,color="w",label="ROI",edgecolor='black')
plt.yscale('log')
#plt.hist(hist_roi[1][:-1],bins=hist_roi[1],weights=hist_roi[0],color='w',label="ROI",edgecolor="black")
plt.ylabel('Number of heads / area [m^-2]')
#ax1.cla() # clear the axis 
#ax2 = ax1.twinx()
for i in range(0,7):
    plt.hist(vcount[i][1][:-1],bins=vbininfo[i][1],weights=vharea[i],color=cm.jet(i*.15),alpha=0.2)
    plt.hist(vcount[i][1][:-1],bins=vbininfo[i][1],weights=vharea[i],label=coords[i],color=cm.jet(i*.15),histtype=u'step')
    plt.plot(center_bins[i],vharea[i],marker='*',color=cm.jet(i*.15),linestyle = 'None')
plt.hist(wwcount1[1][:-1],bins=wwbininfo1[1],weights=wwharea1,label='Warm-Wet Sigmoid 1')
plt.hist(wwcount2[1][:-1],bins=wwbininfo2[1],weights=wwharea2,label='Warm-Wet Sigmoid 2')
#plt.hist(iccount[1][:-1],bins=ICbininfo[1],weights=icharea,label='Icy-Cold')  
plt.xlabel('Elevation (meters)')
#plt.ylabel('Number of heads / area [m^-2] in Individual Valleys')
plt.xlim((-4000,5000))

plt.title('Normalized Head Distributions')
fig.legend(loc="upper left",bbox_to_anchor =(.12,.89))
plt.show()

fig, ax1 = plt.subplots()
#plt.hist(roicount[1][:-1],bins=roicount[1],weights=roicount[0]/roi_area,color="w",label="ROI",edgecolor='black')
plt.yscale('log')
#plt.hist(hist_roi[1][:-1],bins=hist_roi[1],weights=hist_roi[0],color='w',label="ROI",edgecolor="black")
plt.ylabel('Number of heads / area [m^-2]')
#ax1.cla() # clear the axis 
#ax2 = ax1.twinx()
i=2
plt.hist(vcount[i][1][:-1],bins=vbininfo[i][1],weights=vharea[i],color=cm.jet(i*.15),alpha=0.2)
plt.hist(vcount[i][1][:-1],bins=vbininfo[i][1],weights=vharea[i],label=coords[i],color=cm.jet(i*.15),histtype=u'step')
plt.plot(center_bins[i],vharea[i],marker='*',color=cm.jet(i*.15),linestyle = 'None')
plt.hist(wwcount1[1][:-1],bins=wwbininfo1[1],weights=wwharea1,label='Warm-Wet Sigmoid 1')
plt.hist(wwcount2[1][:-1],bins=wwbininfo2[1],weights=wwharea2,label='Warm-Wet Sigmoid 2')
#plt.hist(iccount[1][:-1],bins=ICbininfo[1],weights=icharea,label='Icy-Cold')  
plt.xlabel('Elevation (meters)')
#plt.ylabel('Number of heads / area [m^-2] in Individual Valleys')
plt.xlim((-4000,5000))

plt.title('Normalized Head Distributions')
fig.legend(loc="upper left",bbox_to_anchor =(.12,.89))
plt.show()

fig, ax1 = plt.subplots()
#plt.hist(roicount[1][:-1],bins=roicount[1],weights=roicount[0]/roi_area,color="w",label="ROI",edgecolor='black')
plt.yscale('log')
#plt.hist(hist_roi[1][:-1],bins=hist_roi[1],weights=hist_roi[0],color='w',label="ROI",edgecolor="black")
plt.ylabel('Number of heads / area [m^-2]')
#ax1.cla() # clear the axis 
#ax2 = ax1.twinx()
i=2
plt.hist(vcount[i][1][:-1],bins=vbininfo[i][1],weights=vharea[i],color=cm.jet(i*.15),alpha=0.2)
plt.hist(vcount[i][1][:-1],bins=vbininfo[i][1],weights=vharea[i],label=coords[i],color=cm.jet(i*.15),histtype=u'step')
plt.plot(center_bins[i],vharea[i],marker='*',color=cm.jet(i*.15),linestyle = 'None')
plt.hist(wwcount1[1][:-1],bins=wwbininfo1[1],weights=wwharea1,label='Warm-Wet Sigmoid 1')
plt.hist(wwcount2[1][:-1],bins=wwbininfo2[1],weights=wwharea2,label='Warm-Wet Sigmoid 2')
plt.plot(1000*np.ones(10),np.linspace(0,6*1e-9,10),'--',label = 'Ice Stability Line')
#plt.hist(iccount[1][:-1],bins=ICbininfo[1],weights=icharea,label='Icy-Cold')  
plt.xlabel('Elevation (meters)')
#plt.ylabel('Number of heads / area [m^-2] in Individual Valleys')
plt.xlim((-200,1800))
plt.ylim((1e-10,3e-8))
plt.title('Simulation and Mars Comparison')
fig.legend(loc="upper left",bbox_to_anchor =(.12,.89))
plt.show()

fig, ax1 = plt.subplots()
#plt.hist(roicount[1][:-1],bins=roicount[1],weights=roicount[0]/roi_area,color="w",label="ROI",edgecolor='black')
#plt.yscale('log')
#plt.hist(hist_roi[1][:-1],bins=hist_roi[1],weights=hist_roi[0],color='w',label="ROI",edgecolor="black")
plt.ylabel('Number of heads / area [m^-2]')
#ax1.cla() # clear the axis 
#ax2 = ax1.twinx()
i=2
plt.hist(vcount[i][1][:-1],bins=vbininfo[i][1],weights=vharea[i],color=cm.jet(i*.15),alpha=0.2)
plt.hist(vcount[i][1][:-1],bins=vbininfo[i][1],weights=vharea[i],label=coords[i],color=cm.jet(i*.15),histtype=u'step')
plt.plot(center_bins[i],vharea[i],marker='*',color=cm.jet(i*.15),linestyle = 'None')
plt.hist(wwcount1[1][:-1],bins=wwbininfo1[1],weights=wwharea1,label='Warm-Wet Sigmoid 1')
plt.hist(wwcount2[1][:-1],bins=wwbininfo2[1],weights=wwharea2,label='Warm-Wet Sigmoid 2')
#plt.hist(iccount[1][:-1],bins=ICbininfo[1],weights=icharea,label='Icy-Cold')  
plt.plot(1000*np.ones(10),np.linspace(0,6*1e-9,10),'--',label = 'Ice Stability Line')
plt.xlabel('Elevation (meters)')
#plt.ylabel('Number of heads / area [m^-2] in Individual Valleys')
plt.xlim((-500,2000))
#plt.ylim((10e-10,10e-8))

plt.title('Simulation and Mars Comparison')
fig.legend(loc="upper left",bbox_to_anchor =(.12,.89))
plt.show()