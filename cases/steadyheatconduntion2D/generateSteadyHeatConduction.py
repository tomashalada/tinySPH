#!/bin/python

#-------------------------------------------------------------------------------------#
#
# tinySPH: steady heat conduction 2D - generate particles
#
#-------------------------------------------------------------------------------------#
# Parameters of initial state

dp = 0.0025                                 #initial particle distance

solidL = 0.1 - dp                           #length of solid
solidH = 0.1 - dp                           #height of solid
boundaryL = 0.1 + dp                        #length of boundary
boundaryH = 0.1 + dp                        #height of boundary

n_layers = 3                                #number of boundary layers

rho0 = 8000.                                #initial (referential) density
T0 = 100.                                     #initial temperature
TB = 0.                                   #boundary temperature

#-------------------------------------------------------------------------------------#
# Generate interpolation grid

interpolationFactor = 2                     #number of cell per particle distance
interpolationPlaneL = solidL                #length of interpolation plane
interpolationPlaneH = solidH                #heiht of interpolation plane

#-------------------------------------------------------------------------------------#

import numpy as np

#-------------------------------------------------------------------------------------#

solidL_n = int(np.floor(solidL/dp))
solidH_n = int(np.floor(solidH/dp))

boundaryL_n = int(np.floor(boundaryL/dp))
boundaryH_n = int(np.floor(boundaryH/dp))

solid_rx = []; solid_ry = []; solid_rz = []
boundary_rx = []; boundary_ry = []; boundary_rz = []

#-------------------------------------------------------------------------------------#
# Generate solid

for x in range(solidL_n):
  for z in range(solidH_n):
    solid_rx.append(dp*(x + 1))
    solid_ry.append(0.)
    solid_rz.append(dp*(z + 1))


#-------------------------------------------------------------------------------------#
# Generate walls and corresponging ghost nodes for MDBC

# left wall
for l in range(n_layers):
  for z in range(boundaryH_n - 2):
    boundary_rx.append(0. - l*dp)
    boundary_ry.append(0.)
    boundary_rz.append((z+1)*dp)


# bottom wall
for l in range(n_layers):
  for x in range(boundaryL_n - 2):
    boundary_rx.append((x+1)*dp)
    boundary_ry.append(0.) #we use only 2D case
    boundary_rz.append(0. - l*dp)

# save x-pos of last particle to know where continue with another wall
x_last = boundary_rx[-1] + dp

# right wall
for l in range(n_layers):
  for z in range(boundaryH_n - 2):
    boundary_rx.append(x_last + dp*l)
    boundary_ry.append(0.) #we use only 2D case
    boundary_rz.append((z+1)*dp)

z_last = boundary_rz[-1] + dp

# top wall
for l in range(n_layers):
  for x in range(boundaryL_n - 2):
    boundary_rx.append((x+1)*dp)
    boundary_ry.append(0.) #we use only 2D case
    boundary_rz.append(z_last + l*dp)

# generate the corners
def generate90degCorner(x, z, dirx, dirz):
  for l in range(n_layers):
    for k in range(n_layers):
      boundary_rx.append(x + k*dp*dirx)
      boundary_ry.append(0.)
      boundary_rz.append(z + l*dp*dirz)

generate90degCorner(0, 0., -1, -1)
generate90degCorner(0, z_last, -1, +1)
generate90degCorner(x_last, 0., +1, -1)
generate90degCorner(x_last, z_last, +1, +1)

#-------------------------------------------------------------------------------------#

# Write solid particles
with open("heatconduction_solid.ptcs", "w") as f:
  f.write(str(len(solid_rx)) + "\n")
  for i in range(len(solid_rx)):
    f.write(str(round(solid_rx[i], 5)) + " " + str(round(solid_ry[i], 5)) + " " + str(round(solid_rz[i], 5)) + " " + \
            str(round(rho0, 5)) + " " + str(round(T0)) + "\n")

# Write boundary particles
with open("heatconduction_boundary.ptcs", "w") as f:
  f.write(str(len(boundary_rx)) + "\n")
  for i in range(len(boundary_rx)):
    f.write(str(round(boundary_rx[i], 5)) + " " + str(round(boundary_ry[i], 5)) + " " + str(round(boundary_rz[i], 5)) + " " + \
            str(round(rho0, 5)) + " " + str(round(TB)) + "\n")

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
    f.write(str(IN_rx[i]) + " " + str(IN_ry[i]) + " " + str(IN_rz[i]) + " " + \
            str(0.) + " " + str(0.) + " " + str(0.) + " " + \
            str(0.) + " " + str(0.) + "\n")

#-------------------------------------------------------------------------------------#
