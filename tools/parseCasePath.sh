#!/bin/bash

#-------------------------------------------------------------------------------------#
#
# tinySPH: dambreak case
#
#-------------------------------------------------------------------------------------#

#path of case folder
SaveResults=$1
ParamFile=$2

#read config file and create new path string
CasePathLine=$(grep 'caseResults' $ParamFile)
CasePathLineArray=($CasePathLine)
CasePathLineArray[3]=\"$SaveResults\"
NewCasePathLine=""

for i in "${CasePathLineArray[@]}"
do
  NewCasePathLine=$NewCasePathLine" "$i
done
NewCasePathLine=$NewCasePathLine\;

#replace path in file with parameters
sed -i "s|${CasePathLine}|${NewCasePathLine}|g" $ParamFile

#-------------------------------------------------------------------------------------#


