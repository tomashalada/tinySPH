#!/bin/bash

#-------------------------------------------------------------------------------------#
#
# tinySPH: dambreak case
#
#-------------------------------------------------------------------------------------#
### Case settings

SaveResults=/home/tomas/Documents/temp/tinySPH_dev/dambreak2D_withBI
CaseName="CASE_DAMBREAK_2D"

CaseFolder=$(pwd)

#:#-------------------------------------------------------------------------------------#
#:
#:echo --- CASE_DAMBREAK_2D - TinySPH v0.0-7-2022 ---
#:echo - Generate output folders, set parameters.
#:
#:### Prepare folders, copy data
#:
#:rm -r $SaveResults/OUTPUT
#:
#:mkdir $SaveResults/OUTPUT
#:mkdir $SaveResults/OUTPUT/FLUID
#:mkdir $SaveResults/OUTPUT/BOUND
#:mkdir $SaveResults/OUTPUT/INTERPOLATION
#:
#:./../../tools/parseCaseResults.sh $CaseFolder parameters.hpp
#:./../../tools/parseCasePath.sh $SaveResults parameters.hpp
#:
#:echo ... DONE.
#:
#:#-------------------------------------------------------------------------------------#
set -e

echo - Generate geometry and particles.

### Generate
python generateStillWaterGeometryWithBI.py

### Initial condition to VTK
#Fluid
python ../../tools/PTCStoVTK.py -input="stillwater_fluid.ptcs" -output="stillwater_fluid.vtk"
#Bound
python ../../tools/PTCStoVTK.py -input="stillwater_wall.ptcs" -output="stillwater_wall.vtk"

echo ... DONE.

#:#-------------------------------------------------------------------------------------#
#:
#:echo - Run TinySPH $CASE_DAMBREAK_2D.
#:
#:### Run the case.
#:
#:cd ../../build/dambreak2D_withBI
#:make
#:time ./../../build/dambreak2D_withBI/dambreak2D_BI
#:
#:echo ... DONE.
#:
#:#-------------------------------------------------------------------------------------#
#:
#:echo - Postprocess .ptcs files.
#:
#:cd $CaseFolder
#:
#:#fluid
#:for f in $SaveResults/OUTPUT/FLUID/*
#:do
#:	echo $f
#:	fileIn=${f%.*}
#:	fileVTK=$fileIn.vtk
#:	python ../../tools/PTCStoVTK.py -input="$f" -output="$fileVTK"
#:done
#:
#:#Bound
#:for f in $SaveResults/OUTPUT/BOUND/*
#:do
#:	echo $f
#:	fileIn=${f%.*}
#:	fileVTK=$fileIn.vtk
#:	python ../../tools/PTCStoVTK.py -input="$f" -output="$fileVTK"
#:done
#:
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

