IMPORTANT: 
--DO NOT rename or delete the Heatmap_Script_Files directory since it will be accessed by the Heatmap script.
--Individual chromosomes must be in SEPARATE, SORTED SGR files
--The SGR files should be readable by the program. (chmod 444)
--The tab-delimited peaks file cannot have any special characters and must have UNIX-style line separators. Make sure to do "wc peaks_file.txt" before running the script to ensure UNIX can parse the file properly. If "wc peaks_file.txt" returns 0, then the file most likely contains Mac-style line separators. You can use Text Wrangler or the UNIX tr command to replace the MAC new line characters with the UNIX new line characters. 



1. Update the parameters file.

EXAMPLE PARAMETERS FILE:

chromosomes:1,2,3           <----List the chromosomes to run the script on separated by commas (separated by commas)
absolute_path_to_peak_file:/data/home/pxs183/HEATMAP_SCRIPT/peaks.txt       <----  tab delimited 3-column peak file in this format <chr> <start> <stop>
absolute_path_to_sample_dirs:/data/home/HEATMAP/SGR_FILES/           <----absolute path to the directory containing the sample directories.
sample_dirs:DLD1,C101       <---the individual samples located in <absolute_path_to_sample_dires (separated by commas)

2. In command prompt, enter this:

perl <path to getMinMedianMax.pl>getMinMedianMax.pl <parameters file> <number of processors>


or, if you are already in the main script folder,  just this:

perl getMinMedianMax.pl <name of the parameters file> <number of processors>


NOTE: Currently can only use up to 4 processors per sample

----------------------------------------------------------------------------------------------------
OUTPUT:

An output folder will be created for each sample provided on the sample_dirs line. 
Each output folder will be named <sample-name>_OUTPUT and will be created in the main script directory
Each <chr>_OUT_NO_DUPS lists the median SGR signal in a peak window on every line
FINAL_ALL_CHR.txt will have all the individual FINAL_<chromosome name>_NO_DUPS  files appended together. 

The output format is as follows:
chromosome start peak coordinate end peak coordinate min median max
