SUMMARY

* The pipeline analysis is performed individually for each experimental condition.

* An experimental condition is composed of several GEO samples, all associated to a GEO series.

* The information of each sample must be provided in an input file located in folder Input.

* Several condition files can be deposited in folder Input and the pipeline analysis will be performed on all the conditions.

* The analysis results for all condition files given as input are output in folder Results.

* Section 1 explains how to prepare the condition files that will be given as input.

* Section 2 explains how to run the analysis on the condition files created in Section 1.



Section 1. Constructing the input files for the pipeline analysis of one or more experimental conditions.

The information of each experimental condition must be provided in a text file. One file for each condition.

The condition files must be stored in folder Input.

Each condition file must be named with the format below

  <GEO SERIES NUMBER>_-_<NAME OF EXPERIMENTAL CONDITION>_-_<NUMBER OF TOP DRGs FOR CLUSTERING>.txt

All the samples associated with the experimental condition must be provided along with the time points using the format described as follows.

Each line of the file must consist of the accession number one sample followed by a comma (,) and then followed by the time point, including the time unit (e.g., 2 hours).

Example: The file for condition "D10" of GEO series GSE59015 should be named GSE59015_-_D10_-_3000.csv and the contents must be as follows.

GSM1424453,0 hours
GSM1424454,6 hours
GSM1424455,12 hours
GSM1424456,18 hours
GSM1424457,24 hours
GSM1424458,30 hours
GSM1424459,36 hours
GSM1424460,42 hours

The condition files can be written manually, following strictly the format described above. Optionally, this task can be carried out more easily by using script create_input_files.m. In order to do this, create_input_files.m must be run and the instructions provided thereafter by the program must be followed.

Section 2.

Once all condition files are located in folder Input, then the analysis can be started by running pipeline.m. The script will read ALL the condition files in folder Input and run the analysis for each one of them.

Section 3.

Extracting additional comparisons from a group of condition analyses already finished.

The conditions that need to be compared must be provided in a text file that must be located in folder Input.

This file must be named using the format below.

  <GEO SERIES NUMBER>_-_<NAME OF MACRO CONDITION>.txt
  
Example: The following are the contents of a macrocondition of five H3N1 subjects from GSE52428.

H3N1_001
H3N1_002
H3N1_003
H3N1_004
H3N1_005

This file can be constructed manually, or with BASH script prepare_input.sh.

The latter option requires running the script with the following syntax.

./prepare_input.sh <GEO SERIES> <NAME OF MACRO CONDITION>.txt

Example

./prepare_input.sh GSE52428 H1N1.txt

This will read ALL the conditions whose analyses have been completed for the GEO series indicated and write the file with the format described earlier.



