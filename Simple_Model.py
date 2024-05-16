#!/usr/bin/env python
# coding: utf-8

# # Make an idealized initial topography
import numpy as np #Amanda
import matplotlib.pyplot as plt
from landlab import RasterModelGrid, imshow_grid
from landlab.components import FlowAccumulator, LinearDiffuser, FastscapeEroder, DrainageDensity
import pickle

runName='ic_26'

WW=0

# Landscape Parameters
nrows = 972#840
ncols = 567#490#140
dx = 700#200#700.0 # meters
slope_default = 0
wf_default = 0.15
hf_default = .003
pf_default = 0.5
noise_ampl = 20.0
runtime = 600000 # years
dt = 3000 # years
Ksp=.00025
k=.001

# Set Rainfall Rate based on Warm Wet or Icy Cold scenario
if WW:
    rr=0.0143 #m/yr warm wet
else:
    rr=2 #m/yr icy cold

variables=[nrows,ncols,dx,slope_default,wf_default,hf_default,pf_default,noise_ampl,runtime,dt,Ksp,k,rr]

total_discharge=0

def make_profile(x, slope, width_fac, height_fac, position_fac):
    L = np.amax(x)
    pos = position_fac * L
    width = width_fac * L
    y = (x - pos) / width
    height = height_fac * L
    z = 0.5 * height * (1.0 - np.tanh(y))
    z += slope * (L - x)
    return z

plt.figure(1)
x = np.arange(0.0, ncols)
z1 = make_profile(x, 0, 0.15, 0.0065, 0.54)
z2= make_profile(x,slope_default,wf_default,hf_default,pf_default)
plt.plot(x, z1,label='AGU DEM')
plt.plot(x,z2,label=runName,zorder=2)
plt.legend()
plt.show()

def make_init_topo(grid,
                   slope,
                   width_fac,
                   height_fac,
                   position_fac,
                   noise_ampl
                  ):
    nc = grid.number_of_node_columns
    z = grid.add_zeros('topographic__elevation', at='node')
    pz = make_profile(grid.x_of_node[:nc],
                      slope, width_fac, height_fac, position_fac,
                     )
    np.random.seed(123)
    for i in range(grid.number_of_node_rows):
        start = i * nc
        z[start:(start+nc)] = pz + noise_ampl * np.random.rand(nc)


# Grid Parameters
grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
make_init_topo(grid,
               slope_default,
               wf_default,
               hf_default,
               pf_default,
               noise_ampl
              )
z = grid.at_node['topographic__elevation']
grid.add_zeros("surface_water__discharge",at="node")
initial_topo=grid.add_zeros("initial_topography",at="node")
initial_topo[:]=grid.at_node['topographic__elevation']
grid.add_zeros("water__unit_flux_in",at="node")
grid.add_zeros("cumulative_erosion__depth",at="node")
cumulative_erosion = grid.at_node['cumulative_erosion__depth']

grid.set_closed_boundaries_at_grid_edges(False, True, True, True)

K_values = grid.ones(at='node')*Ksp
wse = grid.add_zeros('water_surface__elevation', at='node', clobber=True)

diff=LinearDiffuser(grid,linear_diffusivity=k)  # m2/yr
fa = FlowAccumulator(grid,
                     flow_director='FlowDirectorD8',
                     depression_finder='LakeMapperBarnes',
                     method='D8',
                     fill_surface=wse,
                     redirect_flow_steepest_descent=True,
                     reaccumulate_flow=True
                    )
erode = FastscapeEroder(grid,
                        K_sp=K_values,
                        m_sp=0.25,
                        n_sp = 1,
                        discharge_field = 'surface_water__discharge')

max_node=grid.number_of_elements('node')


## Change this to do warm wet vs/\. icy cold
if WW:
    rain_area=list(range(0,grid.number_of_nodes)) # Warm-Wet rain uniformly over all nodes
else:
    start_node=148400//dx
    print("Rain input at x=",grid.x_of_node[start_node],"m")
    rain_area=list(np.arange(start_node,max_node,ncols))  # Icy Cold rain over x=20000m np.where(grid_ww.x_of_node==49000)


## Set the area that recieves the rainfall rate
for i in range(0,max_node):
    if i in rain_area:
        grid.at_node['water__unit_flux_in'][i] = rr
    else:
        grid.at_node['water__unit_flux_in'][i] = 10**-20

def myQuadPlots(grid,name):
    fig=plt.figure(figsize=(16,22))

    plt.subplot(4,1,1)
    imshow_grid(grid,'water__unit_flux_in',colorbar_label='m/yr') #,shrink=.5
    plt.title('Water Input',fontsize=16)
    #plt.xlabel('Distance in meters',fontsize=16)
    plt.ylabel('Distance in meters',fontsize=16)


    plt.subplot(4,1,2)
    imshow_grid(grid,z,colorbar_label='Elevation in meters')
    plt.title('Eroded Topography',fontsize=16)
    #plt.xlabel('Distance in meters',fontsize=16)
    plt.ylabel('Distance in meters',fontsize=16)

    plt.subplot(4,1,3)
    imshow_grid(grid,'surface_water__discharge',colorbar_label='$m^3/yr$')
    plt.title('Surface Water Discharge',fontsize=16)
    #plt.xlabel('Distance in meters',fontsize=16)
    plt.ylabel('Distance in meters',fontsize=16)

    plt.subplot(4,1,4)
    imshow_grid(grid,cumulative_erosion,colorbar_label='meters')
    plt.title('Cumulative Erosion',fontsize=16)
    plt.xlabel('Distance in meters',fontsize=16)
    plt.ylabel('Distance in meters',fontsize=16)

    fig.savefig(name)
    plt.show()

for i in range(0,runtime+dt,dt):
    diff.run_one_step(dt)
    fa.run_one_step()
    erode.run_one_step(dt)
    cumulative_erosion[:] = (grid.at_node['initial_topography'] - grid.at_node['topographic__elevation'])
    total_discharge += np.sum(grid.field_values('node','water__unit_flux_in'))

Q=[]
fileName=runName+'.png'
Q.append(grid.at_node['surface_water__discharge'])
myQuadPlots(grid,fileName)

print('Total Discharge: ',total_discharge)

cross_section=np.where(grid.x_of_node==49000)
#print(grid.at_node.keys())
#print(grid['node']['topographic__elevation'][5])

h=grid.at_node['topographic__elevation'][cross_section]
yy=grid.y_of_node[cross_section]

fig=plt.figure()
plt.plot(yy,h,'-b')
plt.xlabel('Distance in Meters')
plt.ylabel('Surface Height in Meters')
plt.title('Channel Depths'+runName)
plt.show()

discharge_threshold = 10000000 #.3 m3/s
Q_mask=np.zeros(max_node,dtype=np.uint8)
for i in range(0,max_node):
    if Q[0][i] > discharge_threshold:
        Q_mask[i] = 1
    else:
        Q_mask[i] = 0

imshow_grid(grid,Q_mask)
plt.title('Channel Mask',fontsize=16)
plt.show()

path=r'C:\Users\avste\Local_Documents\Geology\Code\python_scripts\data'
with open(runName+'.pkl', 'wb') as f:
    pickle.dump([grid,Q_mask,variables,total_discharge], f)

    dd = DrainageDensity(grid, channel__mask=Q_mask)
    mean_drainage_density = dd.calculate_drainage_density()
    print(mean_drainage_density)
    
    plt.figure(figsize=(16,8))
    imshow_grid(grid,'surface_to_channel__minimum_distance')
    plt.title('Drainage Density',fontsize=16)
    plt.show()
