1. Running the pipeline analysis.

1.1 Execute main.m and follow instructions.

2. Performing an integrated analysis

2.1 Open a Linux console and change to the main Pipeline directory. Then enter the following.

./prepare_input.sh [GSE] [Outputfile]

The above will look for the results of the analysis carried out on dataset [GSE] and will load all the analyzed conditions into [Outputfile].

2.2 From the Matlab console enter

integrated_analysis

The program will request the user to enter an input file name. This file must be the same from step 2.1.

3. Compare conditions for which the pipeline analysis has been completed.

3.1 Open a Linux console and change to the main Pipeline directory. Then enter the following.

./prepare_input.sh [GSE] [Outputfile]

The above will look for the results of the analysis carried out on dataset [GSE] and will load all the analyzed conditions into [Outputfile].

3.2 From the Matlab console enter

compare

The program will request the user to enter an input file name. This file must be the same from step 3.1.


4. Integrated report of the pipeline analysis for all conditions.

4.1 Open a Linux console and change to the main Pipeline directory. Then enter the following.

./prepare_input.sh [GSE] [Outputfile]

The above will look for the results of the analysis carried out on dataset [GSE] and will load all the analyzed conditions into [Outputfile].

4.2 From the Matlab console enter

compare

The program will request the user to enter an input file name. This file must be the same from step 3.1.