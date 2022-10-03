#!/bin/bash

#-------------------------------------------------------------------------------------#
#
# tinySPH: steady heat conduction case
#
#-------------------------------------------------------------------------------------#
### Case settings

SaveResults=/home/tomas/Documents/temp/devel/heat
CaseName="CASE_STEADY_HEAT_CONDUCTION_2D"

CaseFolder=$(pwd)

#-------------------------------------------------------------------------------------#

echo --- CASE_DAMBREAK_2D - TinySPH v0.0-9-2022 ---
echo - Generate output folders, set parameters.

### Prepare folders, copy data

rm -r $SaveResults/OUTPUT

mkdir $SaveResults/OUTPUT
mkdir $SaveResults/OUTPUT/SOLID
mkdir $SaveResults/OUTPUT/BOUND
mkdir $SaveResults/OUTPUT/INTERPOLATION

./../../tools/parseCaseResults.sh $CaseFolder parameters.hpp
./../../tools/parseCasePath.sh $SaveResults parameters.hpp

echo ... DONE.

#-------------------------------------------------------------------------------------#
set -e

echo - Generate geometry and particles.

### Generate
python generateSteadyHeatConduction.py

### Initial condition to VTK
#Fluid
python ../../tools/PTCStoVTKheatSimple.py -input="heatconduction_solid.ptcs" -output="heatconduction_solid.vtk"
#Bound
python ../../tools/PTCStoVTKheatSimple.py -input="heatconduction_boundary.ptcs" -output="heatconduction_boundary.vtk"

echo ... DONE.

#-------------------------------------------------------------------------------------#

echo - Run TinySPH $CASE_DAMBREAK_2D.

### Run the case.

cd ../../build/heatequation
make
time ./../../build/heatequation/heatequation

echo ... DONE.

#-------------------------------------------------------------------------------------#

echo - Postprocess .ptcs files.

cd $CaseFolder

#fluid
for f in $SaveResults/OUTPUT/SOLID/*
do
	echo $f
	fileIn=${f%.*}
	fileVTK=$fileIn.vtk
	python ../../tools/PTCStoVTKheatSimple.py -input="$f" -output="$fileVTK"
done

#Bound
for f in $SaveResults/OUTPUT/BOUND/*
do
	echo $f
	fileIn=${f%.*}
	fileVTK=$fileIn.vtk
	python ../../tools/PTCStoVTKheatSimple.py -input="$f" -output="$fileVTK"
done

#:#Interpolation
#:for f in $SaveResults/OUTPUT/INTERPOLATION/*
#:do
#:	echo $f
#:	fileIn=${f%.*}
#:	fileVTK=$fileIn.vtk
#:	#python ../../tools/PTCStoVTKstructuredGrid.py -input="$f" -output="$fileVTK"
#:	python ../../tools/structuredGridSimple.py -input="$f" -output="$fileVTK"
#:done
#:
#:echo ... Postprocessing DONE.
#:
#:#-------------------------------------------------------------------------------------#

