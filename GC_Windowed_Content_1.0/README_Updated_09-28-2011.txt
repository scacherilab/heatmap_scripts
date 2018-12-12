IMPORTANT: 
--DO NOT rename or delete the Helper_Scripts directory since it will be accessed by the Heatmap script.
--The fasta files and the regions file should be readable by the program. (chmod 444)
--The tab-delimited regions file cannot have any special characters and must have UNIX-style line separators. Make sure to do "wc regions_file.txt" before running the script to ensure UNIX can parse the file properly. If "wc regions_file.txt" returns 0, then the file most likely contains Mac-style line separators. You can use Text Wrangler or the UNIX tr command to replace the MAC new line characters with the UNIX new line characters

1. Update the parameters file.

EXAMPLE PARAMETERS FILE:
chromosomes:1,2,3           <----List the chromosomes to run the script on (separated by commas)
absolute_path_to_regions_file:/data/home/regions.txt       <---- tab-delimited 2-column regions file in this format <chr><regions midpoint> or tab delimited 3-column regions file in this format <chr> <start> <stop>
absolute_path_to_fasta_files:/data/home/HEATMAP/FASTA_FILES/           <----absolute path to the directory containing the fasta files in this format *chr(num).fa
calculate_region_center?:no       <---"no" if the regions file is a 2 column file, "yes" if the regions file is a 3 column file
num_windows:50           <---------number of windows for each region
each_dir:5000       <------number of bases in each direction from the center of each region. Hence, the size of each window in this case is 2*<each_dir>/<num_windows> or 200 bp.


2. In command prompt, enter this:

./runGCWindowedContent.pl params.txt 4

NOTE: At this time, you can only run this on up to 4 processors

----------------------------------------------------------------------------------------------------
OUTPUT:

The output folder, OUTPUT, will be created in the main heatmap script directory and will contain all output files.
z_scores_final.txt is the z_scored output for every region.

