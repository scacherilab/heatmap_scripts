#!/usr/bin/perl -w
#Author: Alina Saiakhova, Scacheri Lab, Case Western Reserve University
BEGIN { $^W = 0 }
use Data::Dumper;
use Cwd;
my $working_dir = getcwd; #directory from which the script will be executed
my $SCRIPTS_PATH="$working_dir/Heatmap_Script_Files";
my @chromosomes=@sample_chromosomes=();
my @wig_subfolders=();
my $path_to_peak_file=$path_to_wig_files=$params_file=$all_chroms=$calc_peak_center="";
my $num_cores=$num_windows=$dist_each_dir=0;


getParams(@ARGV); #populates @chromosomes, @wig_subfolders, $path_to_peak_file and $path_to_wig_files



print "-----------Main program started on: ".(localtime)."-------------\n";
my $all_chroms_old=$all_chroms;


foreach(@wig_subfolders){
    my $cur_wig_subfolder=$_;
    chdir "$path_to_wig_files/$cur_wig_subfolder";
    my @wig_files=<*.wig>;
    my @temp_chromosomes=@chromosomes;
@sample_chromosomes=();



foreach(my $i=0; $i<scalar @temp_chromosomes;$i++){
   my $chr=$temp_chromosomes[$i];
my $valid=0;
   foreach(@wig_files){
      my $wig_file=$_;
      if($wig_file=~/chr$chr\.wig|chr$chr\_/){
         $valid=1;
         if  (! -r $wig_file){ #if read permissions are set
                     print "NOTE:$wig_file does not have read permissions. Chromosome $chr will be skipped!\n";
                    
                    
                    
                }else{
            push @sample_chromosomes, $chr;
            
         }
      }
   }
   print "WARNING: No wig file for chr$chr found. The chromosome will be skipped\n" if $valid==0;
}
  
    if (scalar @sample_chromosomes > 0){
        
    chdir $working_dir;
   my $out_sample_dir="$working_dir/$cur_wig_subfolder"."_OUTPUT";
if (  -d "$out_sample_dir" ){
   deldir($out_sample_dir);
print "$out_sample_dir will be overwritten! \n";
    
    
}
mkdir $out_sample_dir;

if ($calc_peak_center eq "yes"){
  print "Calculating peak medians and splitting $path_to_peak_file into individual chromosomes...";  
    calcMediansAndPickOutChromosomes($path_to_peak_file, $out_sample_dir);
}
else {
    print "Splitting $path_to_peak_file into individual chromosomes...";
    pickOutChromosomes($path_to_peak_file, $out_sample_dir);
}
print "Done\n";


print "Scheduling parallel processing jobs and running HeatmapWindowMedians on $cur_wig_subfolder...\n";

$all_chroms=getNewAllChroms(@sample_chromosomes) if scalar @sample_chromosomes < scalar @chromosomes;
$all_chroms=$all_chroms_old if scalar @sample_chromosomes == scalar @chromosomes;
 

system("perl $SCRIPTS_PATH/processManager.pl $all_chroms $cur_wig_subfolder $path_to_wig_files $out_sample_dir  $num_windows $dist_each_dir $SCRIPTS_PATH $num_cores");


#######Convert to right format########
chdir $out_sample_dir;
my @temp_out_files=<*_OUT>;
foreach(@temp_out_files){
    my $cur_file=$_;
    print "Converting $cur_file to the right format...";
    removeDuplicates($cur_file, $out_sample_dir);
    convertToRightFormat($cur_file, $out_sample_dir, $num_windows);
    print "Done\n";
   
    
}


`cat $out_sample_dir/FINAL*_NO_DUPS > $out_sample_dir/FINAL_ALL_CHR.txt`;
`rm $out_sample_dir/FINAL*_NO_DUPS`;
print "Z-scoring...";
z_score("$out_sample_dir/FINAL_ALL_CHR.txt", $cur_wig_subfolder, $out_sample_dir);
print "Done\n";
chdir $working_dir;
##Clean up###

    
    }
    
}
@time_data=localtime(time);
print "-----------Main program finished on: ".(localtime)."-------------\n";

sub z_score{
    my $file_to_z_score=$_[0];
    my $sample=$_[1];
    my $output_dir=$_[2];
    open (OUTPUT, ">$output_dir/z_scores_final.txt");
    $max=$num_windows+1;
    
    `R --vanilla $file_to_z_score < $SCRIPTS_PATH/z_score.R`;
    print OUTPUT "chr_coord\t";
    
    for(my $i=1; $i<$num_windows; $i++){
        print OUTPUT "$sample"."_$i\t";
        
    }
    print OUTPUT "$sample"."_$num_windows\n";
    system("cut -f1-$max $output_dir/temp.txt >> $output_dir/z_scores_final.txt");
    
    
    
    
}
sub convertToRightFormat{
    
  my $input_file=$_[0]."_NO_DUPS";
  my $output_dir=$_[1];
  my $num_windows=$_[2];
  my $num_wi=$num_windows;
 
open(INPUT, "<$output_dir/$input_file");
open (OUTPUT, ">$output_dir/FINAL_$input_file");


my %peaks;
my %sizes;

@name=split(/\_/,$input_file);
$chr_name=shift @name;


while (<INPUT>){
   
    my($line)=$_;
   @columns=split (/\t|\s+/, $line);
   if(!exists($peaks{$columns[0]})){
    
    
    $sizes{$columns[0]}='1';
    
   }
    else{
        $sizes{$columns[0]}=$sizes{$columns[0]}+1;
        
    }
   
   $size=$sizes{$columns[0]};
    if($size<=$num_windows ){
        
    push @{$peaks{$columns[0]}}, $columns[1];


    
    }
  
       
   if($size==$num_windows){
    
    @arr=@{$peaks{$columns[0]}};
    print OUTPUT "$chr_name\t$columns[0]\t";
    foreach(@arr)
    {
        print OUTPUT  "$_\t";
    }
    print OUTPUT "\n";
    #new line of code
   delete  $peaks{$columns[0]};
   }
    
    
}
foreach my $key (keys %peaks) {
    
    @arr=@{$peaks{$key}};
    print OUTPUT "$chr_name\t$key\t";
    foreach(@arr)
    {
       
        if($counter<=$num_wi){
        print OUTPUT  "$_\t";}
    }
    print OUTPUT "\n";
}




 close INPUT;
 close OUTPUT;
 
    
    
}
sub removeDuplicates{
    my $file_name=$_[0];
    my $output_dir_name=$_[1];
    
open(INPUT, "<$output_dir_name/$file_name");
open(OUTPUT, ">$output_dir_name/$file_name"."_NO_DUPS");
%repeats=();
while(<INPUT>){
    
    @liner=split(/\t|\s+/);
    if(!exists($repeats{"$liner[0] and $liner[3] and $liner[5]"})){
        $repeats{"$liner[0] and $liner[3] and $liner[5]"}="";
        
        foreach(@liner){
            print OUTPUT "$_\t";
        }
        print OUTPUT "\n";
    }
    
}
    `rm $output_dir_name/$file_name`;
}

sub pickOutChromosomes{
    my $main_peak_file=$_[0];
    my $output_dir=$_[1];
    
    open(PEAK_FILE, "<$main_peak_file");
    while(my $line=<PEAK_FILE>){
        
          if($line!~m/chr[^\t\s]+[\t\s]+\d+.*/){
            
            die ("The peak file is improperly formatted.The peak file must be in this format: chromosome< TAB or space>< midpoint coord>... . The program will exit.\n")
            
        }
    }
    
    foreach(@sample_chromosomes){
       my $chr=$_;
      
    system("awk '\$1==\"chr$chr\"' $main_peak_file > $output_dir/chr$chr.txt") ;#or die "Can't extract chr$chr peaks from $main_peak_file\n";
       system("sort  $output_dir/chr$chr.txt -n -k 2 -u > $output_dir/tmp$chr.txt");
       system("mv $output_dir/tmp$chr.txt $output_dir/chr$chr.txt");
        
        
    }
    
    
    
    
    
}

sub calcMediansAndPickOutChromosomes{
    
    my $main_peak_file=$_[0];
    my $output_dir=$_[1];
    open(INPUT, "<$main_peak_file" ) || die "Problem calculating peak centers\n";
    open(OUTPUT, ">$output_dir/main_peak_file.txt") || die "Problem calculating peak centers\n";
   
    while ($line=<INPUT>){
    if($line!~m/chr[^\t\s]+[\t\s]+\d+[\t\s]+\d+.*/){
            print "$line\n";
            die ("The peak file is improperly formatted. To calculate peak medians, the peak file must be in this format: chromosome< TAB or space>< start coord><TAB or space><stop coord>... . The program will exit.\n")
            
        }
        my @liner=split(/\t|\s+/, $line);
        my $peak_center=eval{($liner[1]+$liner[2])/2};
        print OUTPUT "$liner[0]\t$peak_center\n" ;
        
        
    }
    foreach(@sample_chromosomes){
       my $chr=$_;
      
    `awk '\$1==\"chr$chr\"' $output_dir/main_peak_file.txt > $output_dir/chr$chr.txt` ;#or die "Can't extract chr$chr peaks from $main_peak_file\n";
       system("sort  $output_dir/chr$chr.txt -n -k 2 > $output_dir/tmp$chr.txt");
       system("mv $output_dir/tmp$chr.txt $output_dir/chr$chr.txt");
      
        
    
    
    
}
     system("rm $output_dir/main_peak_file.txt");
}
sub getParams{
    
 ##read commandline parameters   
if($#_<1){die "Usage: ./getPeakWindowMedians params_file.txt num_cores\n";}
$params_file=$_[0];
open(PARAMS_FILE, "<$params_file") || die "Invalid parameters file\n";

 $num_cores=$_[1];
if ($num_cores<0 or $num_cores>4) { die "You must specify between 1 and 4 cores\n";}


##read params file parameters
while($line=<PARAMS_FILE>){
    chomp $line;
    
    if($line =~m/^chromosomes:([\-\w]+)$/){
      
        push @chromosomes, $1;
        $all_chroms=$1;
       
    }
    elsif($line=~m/^chromosomes:(([\-\w]+\,)+[\-\w]+)$/){
      
 @chromosomes=split(/,/, $1);
 $all_chroms=$1;
    }
    elsif($line=~m/^absolute_path_to_peak_file:(.+)$/){
        
        if (! -e $1){
            
            die "peak file $1 doesn't exist\n";
        }
        
        $path_to_peak_file=$1;
        
    }
    elsif($line=~m/^absolute_path_to_sample_dirs:(.+)$/){
        if( ! -d $1){
            
            die "wig directory $1 doesn't exist\n";
        }
        
        my @path_to_arr=split(//, $1);
        $path_to_wig_files="$1/";
      $path_to_wig_files=$1  if ($path_to_arr[$#path_to_arr] eq "/");
     
        
        
        
    }
    elsif ($line=~m/^sample_dirs:([\+\-\w]+)$/) {
        
        if( ! -d "$path_to_wig_files$1"){
            
            die "sample directory $path_to_wig_files$1 doesn't exist\n";
        }
        push @wig_subfolders, $1;
        
        
    }
    elsif($line=~m/^sample_dirs:(([\+\-\w]+\,)+[\+\-\w]+)$/){
       @wig_subfolders=split(/,/, $1);
       my @new_wig_subfolders=();
       
       for( my $i=0; $i<=$#wig_subfolders; $i++){
        
        if(! -d "$path_to_wig_files$wig_subfolders[$i]"){
            
            print "WARNING: sample $wig_subfolders[$i] doesn't exist and will be ignored\n";
            
            
        }
        else {
            
            push @new_wig_subfolders, $wig_subfolders[$i];
        }
        
       }
       @wig_subfolders=@new_wig_subfolders;
       
       
       
       
        
        
    }
    elsif($line=~m/^calculate_peak_center\?:(yes|no)$/){
        
        $calc_peak_center=$1;
       
        
    }
    elsif($line=~m/^num_windows:([\d]+)$/){
       
        
        $num_windows=$1;
        
        
    }
    elsif($line=~m/^dist_from_peak_in_each_direction:(\d+)$/){
        
        $dist_each_dir=$1;
        
       
    }
    else {print "Not recognized as a parameter line: $line\n";}
   
    
}
    scalar @chromosomes>0 || die "You must specify at least 1 chromosome\n";
    $path_to_peak_file ne "" || die "You must provide the full path to a valid regions file\n";
    $path_to_wig_files ne ""|| die "You must provide the full path to sample directories (including all forward slashes)\n";
    scalar @wig_subfolders >0 || die "You must specify at least one sample directory (no forward slashes)\n";
    $calc_peak_center ne "" || die "You must specify if peak medians will need to be calculated\n";
    $num_windows >0 || die "You must provide number of peak windows\n";
    $dist_each_dir >0 || die "You must specify the distance in each direction of the peak the median signal will need to be calculated for\n";
    print "Finished params validation\n";

}


sub getNewAllChroms{
    my $all_chroms;
    
    my @new_chroms=@_;
    
    if (scalar @new_chroms == 1){
        
        $all_chroms="$new_chroms[0]";
    }
    else{
       $all_chroms="$new_chroms[0]";
    
    for(my $i=1; $i<=$#new_chroms;$i++){
        $all_chroms=$all_chroms.",$new_chroms[$i]";
        
        
    }
   
        
        
    }
    
    
    return $all_chroms;
}
sub deldir {
  my $dirtodel = pop;
  my $sep = '/';
  opendir(DIR, $dirtodel);
  my @files = readdir(DIR);
  closedir(DIR);
 
  @files = grep { !/^\.{1,2}/ } @files;
  @files = map { $_ = "$dirtodel$sep$_"} @files;
  @files = map { (-d $_)?deldir($_):unlink($_) } @files;
 
  rmdir($dirtodel);
}

