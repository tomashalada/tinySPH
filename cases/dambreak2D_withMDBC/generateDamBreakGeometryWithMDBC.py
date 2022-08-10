"""
TestCase Dambreak

"""

boxL = 1.61
boxH = 0.8

fluidL = 0.6
fluidH = 0.3

dp = 0.005
#dp = 0.1
rho0 = 1000.
p0 = 0.
n_layers = 3

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

ghost_rx = []
ghost_ry = []
ghost_rz = []

# left wall
for l in range(n_layers):
    for z in range(boxH_n - 1):
        box_rx.append(0. - l*dp)
        box_ry.append(0.) #we use only 2D case
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

x_last = box_rx[-1] + dp #due to discretisation, we need to save last value of bottom wall

# left corner - 1st layer
box_rx.append(0.); box_ry.append(0.); box_rz.append(0.)
ghost_rx.append(dp); ghost_ry.append(0.); ghost_rz.append(dp)
# left corner - 2nd layer
box_rx.append(-dp); box_ry.append(0.); box_rz.append(0.)
ghost_rx.append(2*dp); ghost_ry.append(0.); ghost_rz.append(dp)
box_rx.append(-dp); box_ry.append(0.); box_rz.append(-dp)
ghost_rx.append(2*dp); ghost_ry.append(0.); ghost_rz.append(2*dp)
box_rx.append(0.); box_ry.append(0.); box_rz.append(-dp)
ghost_rx.append(dp); ghost_ry.append(0.); ghost_rz.append(2*dp)
# left corner -3nd layer
box_rx.append(-dp*2); box_ry.append(0.); box_rz.append(0.)
ghost_rx.append(3*dp); ghost_ry.append(0.); ghost_rz.append(dp)
box_rx.append(-dp*2); box_ry.append(0.); box_rz.append(-dp)
ghost_rx.append(3*dp); ghost_ry.append(0.); ghost_rz.append(2*dp)
box_rx.append(-dp*2); box_ry.append(0.); box_rz.append(-dp*2)
ghost_rx.append(3*dp); ghost_ry.append(0.); ghost_rz.append(3*dp)
box_rx.append(-dp); box_ry.append(0.); box_rz.append(-dp*2)
ghost_rx.append(2*dp); ghost_ry.append(0.); ghost_rz.append(3*dp)
box_rx.append(0.); box_ry.append(0.); box_rz.append(-dp*2)
ghost_rx.append(dp); ghost_ry.append(0.); ghost_rz.append(3*dp)

# right corner
box_rx.append(x_last); box_ry.append(0.); box_rz.append(0.)
ghost_rx.append(x_last-dp); ghost_ry.append(0.); ghost_rz.append(dp)
# left corner - 2nd layer
box_rx.append(x_last+dp); box_ry.append(0.); box_rz.append(0.)
ghost_rx.append(x_last-2*dp); ghost_ry.append(0.); ghost_rz.append(dp)
box_rx.append(x_last+dp); box_ry.append(0.); box_rz.append(-dp)
ghost_rx.append(x_last-2*dp); ghost_ry.append(0.); ghost_rz.append(2*dp)
box_rx.append(x_last); box_ry.append(0.); box_rz.append(-dp)
ghost_rx.append(x_last-dp); ghost_ry.append(0.); ghost_rz.append(2*dp)
# left corner -3nd layer
box_rx.append(x_last+dp*2); box_ry.append(0.); box_rz.append(0.)
ghost_rx.append(x_last-dp*3); ghost_ry.append(0.); ghost_rz.append(dp)
box_rx.append(x_last+dp*2); box_ry.append(0.); box_rz.append(-dp)
ghost_rx.append(x_last-dp*3); ghost_ry.append(0.); ghost_rz.append(2*dp)
box_rx.append(x_last+dp*2); box_ry.append(0.); box_rz.append(-dp*2)
ghost_rx.append(x_last-dp*3); ghost_ry.append(0.); ghost_rz.append(3*dp)
box_rx.append(x_last+dp); box_ry.append(0.); box_rz.append(-dp*2)
ghost_rx.append(x_last-dp*2); ghost_ry.append(0.); ghost_rz.append(3*dp)
box_rx.append(x_last); box_ry.append(0.); box_rz.append(-dp*2)
ghost_rx.append(x_last-dp); ghost_ry.append(0.); ghost_rz.append(3*dp)

# right wall - layer 1
for l in range(n_layers):
    for z in range(boxH_n - 1):
        box_rx.append(x_last + dp*l)
        box_ry.append(0.) #we use only 2D case
        box_rz.append((z+1)*dp)
        ghost_rx.append(x_last - dp*(l+1))
        ghost_ry.append(0.)
        ghost_rz.append((z+1)*dp)

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
                str(ghost_rx[i]) + " " + str(ghost_ry[i]) + " " + str(ghost_rz[i]) + "\n")

### Write ghost node
with open("dambreak_ghostNodes.ptcs", "w") as f:
    f.write(str(len(ghost_rx)) + "\n")
    for i in range(len(ghost_rx)):
        f.write(str(ghost_rx[i]) + " " + str(ghost_ry[i]) + " " + str(ghost_rz[i]) + " " + \
                str(0.) + " " + str(0.) + " " + str(0.) + " " + \
                str(rho0) + " " + str(p0) + "\n")


