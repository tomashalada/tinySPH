#!/bin/python

#-------------------------------------------------------------------------------------#
#
# tinySPH: still water case 2D (MDBC) - generate particles
#
#-------------------------------------------------------------------------------------#
# Parameters of initial state

dp = 0.02                                   #initial particle distance

boxL = 0.8+dp                               #length of box
boxH = 0.6                                  #height of box

fluidL = 0.8-dp                             #length of fluid block
fluidH = 0.5                                #height of fluid block

n_layers = 3                                #number of boundary layers

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
ghost_rx = []; ghost_ry = []; ghost_rz = []

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
    box_ry.append(0.)
    box_rz.append((z+1)*dp)

    ghost_rx.append(0. + dp*(l+1))
    ghost_ry.append(0.)
    ghost_rz.append((z+1)*dp)


# bottom wall
for l in range(n_layers):
  for x in range(boxL_n - 2):
    box_rx.append((x+1)*dp)
    box_ry.append(0.) #we use only 2D case
    box_rz.append(0. - l*dp)

    ghost_rx.append((x+1)*dp)
    ghost_ry.append(0.)
    ghost_rz.append(0. + dp*(l+1))

# save x-pos of last particle to know where continue with another wall
x_last = box_rx[-1] + dp

# right wall
for l in range(n_layers):
  for z in range(boxH_n - 1):
    box_rx.append(x_last + dp*l)
    box_ry.append(0.) #we use only 2D case
    box_rz.append((z+1)*dp)
    ghost_rx.append(x_last - dp*(l+1))
    ghost_ry.append(0.)
    ghost_rz.append((z+1)*dp)

# generate the corners
def generate90degCorner(x, z, dirx, dirz):
  for l in range(n_layers):
    for k in range(n_layers):
      box_rx.append(x + k*dp*dirx)
      box_ry.append(0.)
      box_rz.append(z + l*dp*dirz)

      ghost_rx.append(x + (k+1)*dp*dirx*(-1))
      ghost_ry.append(0.)
      ghost_rz.append(z + (l+1)*dp*dirz*(-1))

generate90degCorner(0, 0., -1, -1)
generate90degCorner(x_last, 0., +1, -1)

#-------------------------------------------------------------------------------------#

# Write fluid particles
with open("stillwater_fluid.ptcs", "w") as f:
  f.write(str(len(fluid_rx)) + "\n")
  for i in range(len(fluid_rx)):
    f.write(str(fluid_rx[i]) + " " + str(fluid_ry[i]) + " " + str(fluid_rz[i]) + " " + \
            str(0.) + " " + str(0.) + " " + str(0.) + " " + \
            str(rho0) + " " + str(fluid_p[i]) + "\n")

# Write boundary particles
with open("stillwater_wall.ptcs", "w") as f:
  f.write(str(len(box_rx)) + "\n")
  for i in range(len(box_rx)):
    f.write(str(box_rx[i]) + " " + str(box_ry[i]) + " " + str(box_rz[i]) + " " + \
            str(0.) + " " + str(0.) + " " + str(0.) + " " + \
            str(rho0) + " " + str(p0) + " " + \
            str(ghost_rx[i]) + " " + str(ghost_ry[i]) + " " + str(ghost_rz[i]) + "\n")

# Write ghost node
with open("stillwater_ghostNodes.ptcs", "w") as f:
  f.write(str(len(ghost_rx)) + "\n")
  for i in range(len(ghost_rx)):
    f.write(str(ghost_rx[i]) + " " + str(ghost_ry[i]) + " " + str(ghost_rz[i]) + " " + \
            str(0.) + " " + str(0.) + " " + str(0.) + " " + \
            str(rho0) + " " + str(p0) + "\n")

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
