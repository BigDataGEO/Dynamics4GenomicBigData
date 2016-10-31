1. Running the pipeline analysis.


2. Performing an integrated analysis

2.1 Open a Linux console and change to the main Pipeline directory. Then enter the following.

./prepare_input.sh [GSE] [Outputfile]

The above will look for the results of the analysis carried out on dataset [GSE] and will load all the analyzed conditions into [Outputfile].

2.2 From the Matlab console enter

integrated_analysis

The program will request the user to enter an input file name. This file must be the same from step 2.1.
