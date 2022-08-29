#!/bin/bash

#-------------------------------------------------------------------------------------#
#
# tinySPH: Still water testcases
#
#-------------------------------------------------------------------------------------#
# Still water with MDBC:

SaveResults=/home/tomas/Documents/temp/tinySPH_glob/MDBC
CaseFolder=/home/tomas/Documents/testovaci/cpp/TINYSPH_lib/cases/stillwater2D_withMDBC

geometries="geometryStraightCutHR
geometryStraightCutLR
geometryStraightHR
geometryStraightLR
geometryWithWedgeHR
geometryWithWedgeLR"

for i in $geometries; do

  echo Computing geometry - $i:
  echo $CaseFolder/$i
  echo $CaseFolder/$i/parameters.hpp
  echo $CaseFolder/$i/stillwater_fluid.ptcs
  echo $CaseFolder/$i/stillwater_wall.ptcs
  echo $CaseFolder/$i/stillwater_ghostNodes.ptcs
  echo $CaseFolder/$i/stillwater_interpolationPlane.ptcs

  echo Prepare case settings.
  cp $CaseFolder/$i/parameters.hpp $CaseFolder
  cp $CaseFolder/$i/stillwater_fluid.ptcs $CaseFolder
  cp $CaseFolder/$i/stillwater_wall.ptcs $CaseFolder
  cp $CaseFolder/$i/stillwater_ghostNodes.ptcs $CaseFolder
  cp $CaseFolder/$i/stillwater_interpolationPlane.ptcs $CaseFolder
  echo ... OK

  echo Create results folder.
  rm -r $SaveResults/OUTPUT_$i

  mkdir $SaveResults/OUTPUT_$i
  mkdir $SaveResults/OUTPUT_$i/FLUID
  mkdir $SaveResults/OUTPUT_$i/BOUND
  mkdir $SaveResults/OUTPUT_$i/INTERPOLATION
  echo ... OK

  echo Setup paths.
  cd $CaseFolder
  ./../../tools/parseCaseResults.sh $CaseFolder parameters.hpp
  ./../../tools/parseCasePath.sh $SaveResults/OUTPUT_$i parameters.hpp
  echo ... OK

  echo Compile and run the case.
  cd ../../build/stillwater2D_withMDBC
  make
  time ./../../build/stillwater2D_withMDBC/stillwater2D_mDBC
  echo ... OK

  echo Remove current settings.
  rm $CaseFolder/parameters.hpp
  rm $CaseFolder/stillwater_fluid.ptcs
  rm $CaseFolder/stillwater_wall.ptcs
  rm $CaseFolder/stillwater_ghostNodes.ptcs
  rm $CaseFolder/stillwater_interpolationPlane.ptcs
  echo ... OK
  echo ''

done

#-------------------------------------------------------------------------------------#
