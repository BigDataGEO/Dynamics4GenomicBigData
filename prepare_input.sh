#!/bin/bash

# Usage example:
# ./prepare_input.sh GSE52428 input.csv

series=$1

outputFile=$2

conditionsFolder="Results/$series/Conditions/"

echo "series,condition" > "$outputFile"

for condition in "$conditionsFolder"/*; do
  echo "$series,`basename "$condition"`" >> "$outputFile"
done

mv $outputFile Input