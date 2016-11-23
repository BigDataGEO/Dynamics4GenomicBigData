1. Running the pipeline analysis for one or more experimental conditions.

The information of each experimental condition must be provided in a text file. One file for each condition.

The condition files must be stored in folder Input.

Each condition file must be named with the format below

  <GEO SERIES NUMBER>_<NAME OF EXPERIMENTAL CONDITION>_<NUMBER OF TOP DRGs FOR CLUSTERING>.txt

All the samples associated with the experimental condition must be provided along with the time points using the format described as follows.

Each line of the file must consist of the accession number one sample followed by a comma (,) and then followed by the time point, including the time unit (e.g., 2 hours).

Example: The file for condition "D10" of GEO series GSE59015 should be named GSE59015_D10_3000.csv and the contents must be as follows.

GSM1424453,0 hours
GSM1424454,6 hours
GSM1424455,12 hours
GSM1424456,18 hours
GSM1424457,24 hours
GSM1424458,30 hours
GSM1424459,36 hours
GSM1424460,42 hours

Once all condition files are located in folder Input, then the analysis can be started by running pipeline.m. The script will read all the condition files in folder Input and run the analysis for each one of them.