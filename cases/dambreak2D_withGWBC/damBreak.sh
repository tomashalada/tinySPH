#!/bin/bash

#-------------------------------------------------------------------------------------#
#
# tinySPH: dambreak case
#
#-------------------------------------------------------------------------------------#
### Case settings

SaveResults=/home/tomas/Documents/temp/tinySPH/dambreak2D_withGWBC
CaseName="CASE_DAMBREAK_2D"

CaseFolder=$(pwd)

#-------------------------------------------------------------------------------------#

echo --- CASE_DAMBREAK_2D - TinySPH v0.0-7-2022 ---
echo - Generate output folders, set parameters.

### Prepare folders, copy data

rm -r $SaveResults/OUTPUT

mkdir $SaveResults/OUTPUT
mkdir $SaveResults/OUTPUT/FLUID
mkdir $SaveResults/OUTPUT/BOUND

./../../tools/parseCaseResults.sh $CaseFolder parameters.hpp
./../../tools/parseCasePath.sh $SaveResults parameters.hpp

echo ... DONE.

#-------------------------------------------------------------------------------------#
set -e

echo - Generate geometry and particles.

### Generate
python generateDamBreakGeometryWithDBC.py

### Initial condition to VTK
#Fluid
python ../../tools/PTCStoVTK.py -input="dambreak_fluid.ptcs" -output="dambreak_fluid.vtk"
#Bound
python ../../tools/PTCStoVTK.py -input="dambreak_wall.ptcs" -output="dambreak_wall.vtk"

echo ... DONE.

#-------------------------------------------------------------------------------------#

echo - Run TinySPH $CASE_DAMBREAK_2D.

### Run the case.

cd ../../build/dambreak2D_withGWBC
make
time ./../../build/dambreak2D_withGWBC/dambreak2D_GWBC

echo ... DONE.

#-------------------------------------------------------------------------------------#

echo - Postprocess .ptcs files.

cd $CaseFolder

#Fluid
for f in $SaveResults/OUTPUT/FLUID/*
do
	echo $f
	fileIn=${f%.*}
	fileVTK=$fileIn.vtk
	python ../../tools/PTCStoVTK.py -input="$f" -output="$fileVTK"
done

#Bound
for f in $SaveResults/OUTPUT/BOUND/*
do
	echo $f
	fileIn=${f%.*}
	fileVTK=$fileIn.vtk
	python ../../tools/PTCStoVTK.py -input="$f" -output="$fileVTK"
done

echo ... Fluid postprocessing DONE.

#-------------------------------------------------------------------------------------#

