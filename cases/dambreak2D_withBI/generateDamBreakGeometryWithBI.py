"""
TestCase Dambreak

"""

boxL = 1.61
boxH = 0.8

fluidL = 0.6
fluidH = 0.3

#dp = 0.005
dp = 0.005
#dp = 0.1
rho0 = 1000.
p0 = 0.
n_layers = 1

import numpy as np

boxL_n = round(boxL/dp)
boxH_n = round(boxH/dp)

fluidL_n = round(fluidL/dp)
fluidH_n = round(fluidH/dp)


### Generate fluid particles
fluid_rx = np.zeros(fluidL_n*fluidH_n)
fluid_ry = np.zeros(fluidL_n*fluidH_n)
fluid_rz = np.zeros(fluidL_n*fluidH_n)


for x in range(fluidL_n):
    for z in range(fluidH_n):
        fluid_rx[x*fluidH_n + z] = dp*(x + 1)
        fluid_ry[x*fluidH_n + z] = 0. #we use only 2D case
        fluid_rz[x*fluidH_n + z] = dp*(z + 1)


### Generate boundary particles
box_rx = []
box_ry = []
box_rz = []
#box normals
box_nx = []
box_ny = []
box_nz = []

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

x_last = box_rx[-1 -(n_layers - 1)] #due to discretisation, we need to save last value of bottom wall

# corners
box_rx.append(0.)
box_ry.append(0.) #we use only 2D case
box_rz.append(0.)
box_nx.append(np.sqrt(0.5))
box_ny.append(0.)
box_nz.append(np.sqrt(0.5))

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

### Write fluid particles
with open("dambreak_fluid.ptcs", "w") as f:
    f.write(str(len(fluid_rx)) + "\n")
    for i in range(len(fluid_rx)):
        f.write(str(fluid_rx[i]) + " " + str(fluid_ry[i]) + " " + str(fluid_rz[i]) + " " + \
                str(0.) + " " + str(0.) + " " + str(0.) + " " + \
                str(rho0) + " " + str(p0) + "\n")

### Write boundary particles
with open("dambreak_wall.ptcs", "w") as f:
    f.write(str(len(box_rx)) + "\n")
    for i in range(len(box_rx)):
        f.write(str(box_rx[i]) + " " + str(box_ry[i]) + " " + str(box_rz[i]) + " " + \
                str(0.) + " " + str(0.) + " " + str(0.) + " " + \
                str(rho0) + " " + str(p0) + " " + \
                str(box_nx[i]) + " " + str(box_ny[i]) + " " + str(box_nz[i]) +"\n")

### Generate interpolation nodes
interpolationFactor = 1
interpolationPlaneL = boxL
interpolationPlaneH = fluidH

interpolationL_n = round(interpolationPlaneL/(dp/interpolationFactor))
interpolationH_n = round(interpolationPlaneH/(dp/interpolationFactor))

IN_rx = []
IN_ry = []
IN_rz = []

for x in range(interpolationL_n):
    for z in range(interpolationH_n):
        IN_rx.append((dp/interpolationFactor)*(x + 1))
        IN_ry.append(0.) #we use only 2D case
        IN_rz.append((dp/interpolationFactor)*(z + 1))


### Write boundary particles
with open("dambreak_interpolationPlane.ptcs", "w") as f:
    f.write(str(len(IN_rx)) + "\n")
    for i in range(len(IN_rx)):
        f.write(str(IN_rx[i]) + " " + str(IN_ry[i]) + " " + str(IN_rz[i]) + " " + \
                str(0.) + " " + str(0.) + " " + str(0.) + " " + \
                str(0.) + " " + str(0.) + "\n")

