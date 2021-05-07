# VERSION 1.0.            MOHAMMED EL-SAYED                                      
# COPYRIGHT (C) 2020       POTSDAM 
# MY PERSONAL CODE FOR THE FRIENDS-OF-FRIENDS (FOF) CLUSTERING ALGORITHM TO IDENTIFY THE STRUCTRE OF COSMIC V-WEB cosmic velocity (V) web from Simulation Illustris-TNG




import numpy as np 
import sys 
import scipy.io 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.colors as mcol
from scipy import ndimage

file='CIC_135_128_2.0_ShearTensor_Eigen_cell_web.sav' 
s=scipy.io.readsav(file) 
web_type=s['web_type']

#sys.stdout = open("cosmic.txt", "w") #cosmic_web_catalogue_knots

def get_cosmic_web_catalogue_3d(web_type):
 
      neighboring_cells=[] #neighboring_cells = np.zeros((128,128,128))
      for i in range(len(web_type)): 
         for j in range (len(web_type[i])): 
             for k, cell in enumerate(web_type[i][j]):
                 #if cell==2:           
                      
                        if i == 0 or i == len(web_type) - 1 or j == 0 or j == len(web_type[i]) - 1 or k == 0 or k==len(web_type[i][j]) - 1 : 

                           # Edges Or Corners
                           # Link Each Center Cell With Any Other Cell With Which It Shares An Edge Or A Corner
                           neu_neighboring_cells =[]
                           id_grp=[]

                           if i != 0:
                              neu_neighboring_cells.append(web_type[i - 1][k][j])
                              id_grp.append(j * len(web_type[i-1][j]) + k)  

                           if i != len(web_type) - 1:
                              neu_neighboring_cells.append(web_type[i + 1][k][j])  
                              id_grp.append(j * len(web_type[i+1][j]) + k)
                             
                           if k != 0:
                              neu_neighboring_cells.append(web_type[i][k - 1][j])  
                              id_grp.append(j * len(web_type[i][j]) + (k-1))

                           if k != len(web_type[i][j]) - 1:
                                
                              neu_neighboring_cells.append(web_type[i][k + 1][j])
                              id_grp.append(j * len(web_type[i][j]) + (k+1))
                           if j != 0:
                              neu_neighboring_cells.append(web_type[i][k][j - 1])  
                              id_grp.append((j-1) * len(web_type[i][j]) + k+1)                        
                           if j != len(web_type[i]) - 1:
                              neu_neighboring_cells.append(web_type[i][k][j + 1])  
                              id_grp.append( (j+1) * len(web_type[i][j+1]) + k)
                           if i != 0 and k != 0 and j != 0:
                              neu_neighboring_cells.append(web_type[i-1][k-1][j-1])
                              id_grp.append((j-1) * len(web_type[i-1][j-1]) + (k-1))   
                               
                           if i != 0 and k != 0 and j != len(web_type[i]) - 1:
                              neu_neighboring_cells.append(web_type[i-1][k-1][j+1])     
                              id_grp.append((j+1) * len(web_type[i-1][j+1]) + (k-1))
                           if i != 0 and k != len(web_type[i][j])-1 and j != 0:   
                              neu_neighboring_cells.append(web_type[i-1][k+1][j-1])   
                              id_grp.append((j-1) * len(web_type[i-1][j-1]) + (k+1) )
                           if i != 0 and k != len(web_type[i][j])-1 and j != len(web_type[i])-1:
                              neu_neighboring_cells.append(web_type[i-1][k+1][j+1])   
                              id_grp.append((j+1) * len(web_type[i-1][j+1]) + (k+1) )
                           if i != len(web_type)-1 and k != 0 and j != 0:
                              neu_neighboring_cells.append(web_type[i+1][k-1][j-1])   
                              id_grp.append((j-1) * len(web_type[i+1][j-1]) + (k-1))
                           if i != len(web_type)-1 and k != 0 and j != len(web_type[i])-1:
                              neu_neighboring_cells.append(web_type[i+1][k-1][j+1])  
                              id_grp.append((j+1) * len(web_type[i+1][j+1]) + (k-1) )
                    
                           if i != len(web_type)-1 and k != len(web_type[i][j])-1 and j != 0:
                              neu_neighboring_cells.append(web_type[i+1][k+1][j-1]) 
                              id_grp.append( (j-1) * len(web_type[i+1][j-1]) + (k+1) ) 
                    
                           if i != len(web_type)-1 and k != len(web_type[i][j])-1 and j != len(web_type[i])-1 :  
                              neu_neighboring_cells.append(web_type[i+1][k+1][j+1])   
                              id_grp.append((j+1) * len(web_type[i+1][j+1]) + (k+1) )
                           if i != 0 and k != 0:
                              neu_neighboring_cells.append(web_type[i-1][k-1][j])  
                              id_grp.append((j) * len(web_type[i-1][j]) + (k-1) )
                              
                           if i != 0 and  k != len(web_type[i][j])-1:  
                              neu_neighboring_cells.append(web_type[i - 1][k+1][j])
                              id_grp.append((j) * len(web_type[i-1][j]) + (k+1) )
                           if i != len(web_type)-1 and k != 0:
                              neu_neighboring_cells.append(web_type[i + 1][k-1][j])
                              id_grp.append((j) * len(web_type[i+1][j]) + (k-1) )
                           if i != len(web_type)-1 and k != len(web_type[i][j])-1:
                              neu_neighboring_cells.append(web_type[i + 1][k+1][j])
                              id_grp.append((j) * len(web_type[i+1][j]) + (k+1) )
                           if i != 0 and j != 0: 
                              neu_neighboring_cells.append(web_type[i - 1][k][j-1])
                              id_grp.append((j-1) * len(web_type[i-1][j-1]) + k )
                           if i != 0 and j !=len(web_type[i])-1: 
                              neu_neighboring_cells.append(web_type[i - 1][k][j+1])
                              id_grp.append((j+1) * len(web_type[i-1][j+1]) + k )
                           if i != len(web_type)-1 and j != 0:  
                              neu_neighboring_cells.append(web_type[i + 1][k][j-1])
                              id_grp.append((j-1) * len(web_type[i+1][j-1]) + k )
                           if i != len(web_type)-1 and j !=len(web_type[i])-1:   
                              neu_neighboring_cells.append(web_type[i + 1][k][j+1])
                              id_grp.append((j+1) * len(web_type[i+1][j+1]) + k)
                              
                           if k != len(web_type[i][j])-1 and j !=len(web_type[i])-1:  
                              neu_neighboring_cells.append(web_type[i][k+1][j+1])
                              id_grp.append((j+1) * len(web_type[i][j+1]) + (k+1) )
                           if k != len(web_type[i][j])-1 and j != 0:  
                              neu_neighboring_cells.append(web_type[i][k+1][j-1])
                              id_grp.append((j-1) * len(web_type[i][j-1]) + (k+1))
                           if k != 0 and j != 0:  
                              neu_neighboring_cells.append(web_type[i][k-1][j-1])
                              id_grp.append((j-1) * len(web_type[i][j-1]) + (k-1) )
                           if k != 0 and j !=len(web_type[i])-1:  
                              neu_neighboring_cells.append(web_type[i][k-1][j+1])    
                              id_grp.append((j+1) * len(web_type[i][j+1]) + (k-1) )
                           #neu_neighboring_cells=np.array(neu_neighboring_cells)         
                           
                        else: 
                   # Add 26 Neighboring_Cells: Identify If Any One Of The 26 (Ie 3x3x3-1) Neighbouring Cells
                   # Each Cell Has 26 Neighboring Cells(a cube Whose size is  3x3x3 cells). 128^3≈2M Is The Total Number Of Cells In The Entire Simulation Box 
                   # Linking List Of Neighboring Cells     
                   
                          neu_neighboring_cells = np.array([ 
                              web_type[i - 1][k][j],       # j * len(web_type[i-1][j]) + k
                              web_type[i + 1][k][j],  
                              web_type[i][k - 1][j],  
                              web_type[i][k + 1][j],    
                              web_type[i][k][j - 1],   
                              web_type[i][k][j + 1],  
                              web_type[i-1][k-1][j-1],  
                              web_type[i-1][k-1][j+1],  
                              web_type[i-1][k+1][j-1], 
                              web_type[i-1][k+1][j+1], 
                              web_type[i+1][k-1][j-1], 
                              web_type[i+1][k-1][j+1], 
                              web_type[i+1][k+1][j+1], 
                              web_type[i+1][k+1][j-1],   
                              web_type[i - 1][k-1][j], 
                              web_type[i - 1][k+1][j],
                              web_type[i + 1][k-1][j],
                              web_type[i + 1][k+1][j],
                              web_type[i - 1][k][j-1],
                              web_type[i - 1][k][j+1],
                              web_type[i + 1][k][j-1],
                              web_type[i + 1][k][j+1],
                              web_type[i][k+1][j+1],
                              web_type[i][k+1][j-1],
                              web_type[i][k-1][j-1],
                              web_type[i][k-1][j+1]
                           ])
                           
                          id_grp= [
                             j * len(web_type[i-1][j]) + k,
                             j * len(web_type[i+1][j]) + k,
                             j * len(web_type[i][j]) + (k-1),
                             j * len(web_type[i][j]) + (k+1),
                             j * len(web_type[i][j-1]) + k,
                             j * len(web_type[i][j+1]) + k,
                             (j-1) * len(web_type[i-1][j-1]) + (k-1),
                             (j+1) * len(web_type[i-1][j+1]) + k-1, 
                             (j-1) * len(web_type[i-1][j-1]) + (k+1),
                             (j+1) * len(web_type[i-1][j+1]) + k+1,
                             (j-1) * len(web_type[i+1][j-1]) + k-1,
                             (j-1) * len(web_type[i+1][j+1]) + k-1,
                             (j+1) * len(web_type[i+1][j+1]) + k+1,
                             (j-1) * len(web_type[i+1][j-1]) + k+1,
                             j * len(web_type[i-1][j]) + (k-1),
                             j * len(web_type[i-1][j]) + (k+1),
                             j * len(web_type[i+1][j]) + k-1,
                             j * len(web_type[i+1][j]) + k+1,
                             (j-1) * len(web_type[i-1][j-1]) + k,
                             (j+1) * len(web_type[i-1][j+1]) + k,
                             (j-1) * len(web_type[i+1][j-1]) + k,
                             (j+1) * len(web_type[i+1][j+1]) + k,
                             (j+1) * len(web_type[i][j+1]) + k+1,
                             (j-1) * len(web_type[i][j-1]) + (k+1),
                             (j-1) * len(web_type[i][j-1]) + (k-1),
                             (j+1) * len(web_type[i][j+1]) + (k-1) 
                             ]

                        
                           neighboring_cells.append({"index_central": j * len(web_type[i][j]) + k, "central_cell": cell,"id_cells": id_grp ,"group_list": neu_neighboring_cells})     
                        
      return neighboring_cells 
print(get_cosmic_web_catalogue_3d(web_type))
ll=np.array(get_cosmic_web_catalogue_3d(web_type)) 
l0=np.zeros((128*128*128)) 
for r in range(128*128*128):
             row=ll[r]
             l0[r]=row[0]
l1=l0.reshape(128,128,128) 
l2=l1[:,0,:]

import matplotlib.colors as mc

def addNorm(cmapData):
    cmapData['norm'] = mc.BoundaryNorm(cmapData['bounds'], cmapData['cmap'].N)
    return True
def discretize(cmap, bounds):
    resCmap = {}
    resCmap['cmap'] = mc.ListedColormap( \
        [cmap(i/len(bounds)) for i in range(len(bounds))]
    )
    resCmap['bounds'] = bounds
    addNorm(resCmap)
    return resCmap
    
    
levels = [0,1,2,3,4,5]
cmapData = discretize(plt.cm.jet, bounds=levels)
plt.contourf(l2.T, levels=levels, cmap=cmapData['cmap'], norm=cmapData['norm'])

plt.contourf(l2.T, levels=None, cmap=None, norm=None)
plt.colorbar()
cbar.set_ticks([0,1,2,3])
cbar.set_ticklabels(["0: knot", "1: filament", "2: wall", "3: void"])
plt.title('neighbouring_cells_3d(web_type)') 
plt.show()

#sys.stdout.close() 







