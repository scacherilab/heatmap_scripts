IMPORTANT: 
--DO NOT rename or delete the Heatmap_Script_Files directory since it will be accessed by the Heatmap script.
--Individual chromosomes must be in SEPARATE, SORTED SGR files
--The SGR files and the peak file should be readable by the program. (chmod 444)
--The tab-delimited peaks file cannot have any special characters and must have UNIX-style line separators. Make sure to do "wc peaks_file.txt" before running the script to ensure UNIX can parse the file properly. If "wc peaks_file.txt" returns 0, then the file most likely contains Mac-style line separators. You can use Text Wrangler or the UNIX tr command to replace the MAC new line characters with the UNIX new line characters

1. Update the parameters file.

EXAMPLE PARAMETERS FILE:

chromosomes:1,2,3           <----List the chromosomes to run the script on separated by commas (separated by commas)
absolute_path_to_peak_file:/data/home/pxs183/HEATMAP_SCRIPT/peaks.txt       <---- tab-delimited 2-column peak file in this format <chr><peak midpoint> or tab delimited 3-column peak file in this format <chr> <start> <stop>
absolute_path_to_sample_dirs:/data/home/HEATMAP/SGR_FILES/           <----absolute path to the directory containing the sample directories.
sample_dirs:DLD1,C101       <---the individual samples located in <absolute_path_to_sample_dires (separated by commas)
calculate_peak_center?:no       <---"no" if the peak file is a 2 column file, "yes" if the peak file is a 3 column file
num_windows:50           <---------number of windows for each peak
each_dir:5000       <------number of bases in each direction from the center of each peak. Hence, the size of each window in this case is 2*<each_dir>/<num_windows> or 200 bp.


2. In command prompt, enter this:

<path to runHeatmapWindowMedians>/./runHeatmapWindowMedians.pl params.txt 4

or, if you are already in the main script directory, enter this:

perl ./runHeatmapWindowMedians.pl params.txt 4


NOTE: At this time, you can only run this on up to 4 processors

----------------------------------------------------------------------------------------------------
OUTPUT:

An output folder will be created for each sample provided on the sample_dirs line. 
Each output folder will be named <sample-name>_OUTPUT and will be created in the main heatmap script directory
Each <chr>_OUT_NO_DUPS lists the median SGR signal in a peak window on every line
FINAL_ALL_CHR.txt will have all of the individual FINAL_<chromosome name>_NO_DUPS  files appended together. 
The z-scores output file for a sample is z_scores_final.txt

