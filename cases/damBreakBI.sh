#!/bin/bash

#-------------------------------------------------------------------------------------#
#
# tinySPH: DamBreak testcases
#
#-------------------------------------------------------------------------------------#
# Still water with BI:

#SaveResults=/home/tomas/Documents/temp/tinySPH_glob/BI
SaveResults=/home/tomas/Documents/temp/devel
CaseFolder=/home/tomas/Documents/testovaci/cpp/TINYSPH_lib/cases/dambreak2D_withBI

geometries="geometryLR"

for i in $geometries; do

  echo Computing geometry - $i:
  echo $CaseFolder/$i
  echo $CaseFolder/$i/parameters.hpp
  echo $CaseFolder/$i/dambreak_fluid.ptcs
  echo $CaseFolder/$i/dambreak_wall.ptcs
  echo $CaseFolder/$i/dambreak_interpolationPlane.ptcs

  echo Prepare case settings.
  cp $CaseFolder/$i/parameters.hpp $CaseFolder
  cp $CaseFolder/$i/dambreak_fluid.ptcs $CaseFolder
  cp $CaseFolder/$i/dambreak_wall.ptcs $CaseFolder
  cp $CaseFolder/$i/dambreak_interpolationPlane.ptcs $CaseFolder
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
  cd ../../build/dambreak2D_withBI
  make
  time ./../../build/dambreak2D_withBI/dambreak2D_BI
  echo ... OK

  echo Remove current settings.
  rm $CaseFolder/parameters.hpp
  rm $CaseFolder/dambreak_fluid.ptcs
  rm $CaseFolder/dambreak_wall.ptcs
  rm $CaseFolder/dambreak_interpolationPlane.ptcs
  echo ... OK
  echo ''

done

#-------------------------------------------------------------------------------------#
