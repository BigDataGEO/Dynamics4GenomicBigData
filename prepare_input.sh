#!/bin/bash

# Usage example:
# ./prepare_input.sh GSE52428 input.csv

series=$1

outputFile=$series"_-_$2"

conditionsFolder="Output/$series/Conditions/"

# echo "condition" > $outputFile

for condition in "$conditionsFolder"/*; do
  echo "`basename "$condition"`" >> "$outputFile"
done

mv $outputFile Input