import numpy as np
import vtk

import sys
import argparse
from contextlib import redirect_stdout


# ------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument("-input")
parser.add_argument("-output")

args = parser.parse_args()
inputFileName = args.input
outputFileName = args.output

if args.input == None:
    print("Error: Please insert name of input .ptcs file to read!")
    exit()

if args.output == None:
    print("Error: Please insert name of output .vtk file to write!")
    exit()

print("Input file: ", inputFileName)
print("Output file: ", outputFileName)

# ------------------------------------------------

with open(inputFileName, "r") as f:
    #DATA = f.read()
    LINES =  f.readlines()

#print(LINES)
#print(type(LINES))

#ptcs = LINES[1:]
ptcs = LINES
#print(ptcs)

DATA = np.genfromtxt(ptcs, delimiter=" ")
#print(DATA)
#print(type(DATA))
#print(DATA.shape)

r = DATA[:,:3]
rho = DATA[:,3]
T = DATA[:,4]

xdim = np.unique(r[:,0]).size
ydim = np.unique(r[:,1]).size
zdim = np.unique(r[:,2]).size

print("r_size: ", r[:,0].size)
print("x_dim: ", xdim, " y_dim: ", ydim, " z_dim: ", zdim)

r_resampled = r.copy()
rho_resampled = rho.copy()
T_resampled = T.copy()

def resample():
    for z in range(zdim):
        for x in range(xdim):
            r_resampled[xdim*z + x,0] = r[x*zdim + z,0]
            r_resampled[xdim*z + x,1] = r[x*zdim + z,1]
            r_resampled[xdim*z + x,2] = r[x*zdim + z,2]

            rho_resampled[xdim*z + x] = rho[x*zdim + z]
            T_resampled[xdim*z + x] = T[x*zdim + z]


def myGrid():
    with open(outputFileName, 'w') as f:
        with redirect_stdout(f):
            print("# vtk DataFile Version 3.0")
            print("vtk output")
            print("ASCII")
            print("DATASET STRUCTURED_GRID")
            print("DIMENSIONS ", xdim, " ", ydim, " ", zdim)
            print("POINTS ", xdim*ydim*zdim, " float")
            for i in range(xdim*ydim*zdim):
                print("", r_resampled[i,0], " ", r_resampled[i,1], " ", r_resampled[i,2])
            print("POINT_DATA ", xdim*ydim*zdim)
            print("FIELD FieldData 2")
            print("DENSITY 1", xdim*ydim*zdim, " float")
            for i in range(xdim*ydim*zdim):
                print("", rho_resampled[i])
            print("TEMPERATURE 1", xdim*ydim*zdim, " float")
            for i in range(xdim*ydim*zdim):
                print("", T_resampled[i])

resample()
myGrid()

