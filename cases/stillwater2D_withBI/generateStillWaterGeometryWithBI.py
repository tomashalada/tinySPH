#!/bin/python

#-------------------------------------------------------------------------------------#
#
# tinySPH: still water case 2D (BT) - generate particles
#
#-------------------------------------------------------------------------------------#
# Parameters of initial state

dp = 0.02                                   #initial particle distance

boxL = 0.8+dp                               #length of box
boxH = 0.6                                  #height of box

fluidL = 0.8-dp                             #length of fluid block
fluidH = 0.5                                #height of fluid block

n_layers = 1                                #number of boundary layers

rho0 = 1000.                                #initial (referential) density
p0 = 0.                                     #initial pressure (recomputed)

#-------------------------------------------------------------------------------------#
# Generate interpolation grid

interpolationFactor = 2                     #number of cell per particle distance
interpolationPlaneL = fluidL                #length of interpolation plane
interpolationPlaneH = fluidH                #heiht of interpolation plane

#-------------------------------------------------------------------------------------#

import numpy as np

#-------------------------------------------------------------------------------------#

boxL_n = int(np.floor(boxL/dp))
boxH_n = int(np.floor(boxH/dp))
fluidL_n = int(np.floor(fluidL/dp))
fluidH_n = int(np.floor(fluidH/dp))

fluid_rx = []; fluid_ry = []; fluid_rz = []
box_rx = []; box_ry = []; box_rz = []
box_nx = []; box_ny = []; box_nz = []

fluid_p = []; box_p = []                    #fluid is initialized in hydrostatic state

#-------------------------------------------------------------------------------------#
# Generate box of fluid

for x in range(fluidL_n):
  for z in range(fluidH_n):
    fluid_rx.append(dp*(x + 1))
    fluid_ry.append(0.)
    fluid_rz.append(dp*(z + 1))

    fluid_p.append(9.81*(fluidH - fluid_rz[-1])*1000)

#-------------------------------------------------------------------------------------#
# Generate walls and corresponging ghost nodes for MDBC


# left wall
for l in range(n_layers):
    for z in range(boxH_n - 1):
        box_rx.append(0. - l*dp)
        box_ry.append(0.) #we use only 2D case
        box_rz.append((z+1)*dp)

        box_nx.append(1.)
        box_ny.append(0.)
        box_nz.append(0.)


# bottom wall
for l in range(n_layers):
    for x in range(boxL_n + (n_layers - 1)*2 - 2):
        box_rx.append((x-(n_layers - 1) + 1)*dp)
        box_ry.append(0.) #we use only 2D case
        box_rz.append(0. - l*dp)
        box_nx.append(0.)
        box_ny.append(0.)
        box_nz.append(1.)

# save x-pos of last particle to know where continue with another wall
x_last = box_rx[-1 -(n_layers - 1)]

# left corner
box_rx.append(0.)
box_ry.append(0.) #we use only 2D case
box_rz.append(0.)
box_nx.append(np.sqrt(0.5))
box_ny.append(0.)
box_nz.append(np.sqrt(0.5))

# right corner
box_rx.append(x_last + dp)
box_ry.append(0.) #we use only 2D case
box_rz.append(0.)
box_nx.append(-np.sqrt(0.5))
box_ny.append(0.)
box_nz.append(np.sqrt(0.5))

x_last = box_rx[-1 -(n_layers - 1)] #due to discretisation, we need to save last value of bottom wall

# right wall - layer 1
for l in range(n_layers):
    for z in range(boxH_n -1):
        box_rx.append(x_last + dp*l)
        box_ry.append(0.) #we use only 2D case
        box_rz.append((z+1)*dp)
        box_nx.append(-1.)
        box_ny.append(0.)
        box_nz.append(0.)

#-------------------------------------------------------------------------------------#

# Write fluid particles
with open("stillwater_fluid.ptcs", "w") as f:
  f.write(str(len(fluid_rx)) + "\n")
  for i in range(len(fluid_rx)):
    f.write(str(round(fluid_rx[i], 5)) + " " + str(round(fluid_ry[i], 5)) + " " + str(round(fluid_rz[i], 5)) + " " + \
            str(0.) + " " + str(0.) + " " + str(0.) + " " + \
            str(rho0) + " " + str(round(fluid_p[i], 5)) + "\n")

# Write boundary particles
with open("stillwater_wall.ptcs", "w") as f:
    f.write(str(len(box_rx)) + "\n")
    for i in range(len(box_rx)):
        f.write(str(round(box_rx[i], 5)) + " " + str(round(box_ry[i], 5)) + " " + str(round(box_rz[i], 5)) + " " + \
                str(0.) + " " + str(0.) + " " + str(0.) + " " + \
                str(rho0) + " " + str(p0) + " " + \
                str(box_nx[i]) + " " + str(box_ny[i]) + " " + str(box_nz[i]) +"\n")

#-------------------------------------------------------------------------------------#

interpolationL_n = int(np.floor(interpolationPlaneL/(dp/interpolationFactor)))
interpolationH_n = int(np.floor(interpolationPlaneH/(dp/interpolationFactor)))

IN_rx = []; IN_ry = []; IN_rz = []

#-------------------------------------------------------------------------------------#

for x in range(interpolationL_n-1):
  for z in range(interpolationH_n-1):
    IN_rx.append((dp/interpolationFactor)*x + dp)
    IN_ry.append(0.)
    IN_rz.append((dp/interpolationFactor)*z + dp)

#-------------------------------------------------------------------------------------#

# Write boundary particles
with open("stillwater_interpolationPlane.ptcs", "w") as f:
  f.write(str(len(IN_rx)) + "\n")
  for i in range(len(IN_rx)):
    f.write(str(round(IN_rx[i], 5)) + " " + str(round(IN_ry[i], 5)) + " " + str(round(IN_rz[i], 5)) + " " + \
            str(0.) + " " + str(0.) + " " + str(0.) + " " + \
            str(0.) + " " + str(0.) + "\n")

#-------------------------------------------------------------------------------------#
