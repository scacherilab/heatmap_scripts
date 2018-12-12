#!/usr/bin/perl -w
BEGIN { $^W=0}


 $time{$format};

 use Cwd;
my $working_dir = getcwd; #directory from which the script will be executed
my $SCRIPTS_PATH="$working_dir/Heatmap_Script_Files";
my @chromosomes=@sample_chromosomes=();
my @wig_subfolders=();
my $path_to_peak_file=$path_to_wig_files=$params_file=$all_chroms=$has_overlapping_peaks=$single_wig_file_name="";
my $num_cores=0;






getParams(@ARGV); #populates @chromosomes, @wig_subfolders, $path_to_peak_file and $path_to_wig_files




print "-----------Main program started on: ".(localtime)."-------------\n";


#if single wig file, then run splitWigIntoChr.pl
my $all_chroms_old=$all_chroms;
foreach(@wig_subfolders){
    my $cur_wig_subfolder=$_;
    chdir "$path_to_wig_files/$cur_wig_subfolder";
    my @wig_files=<*.wig>;
     @temp_chromosomes=@chromosomes;
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
   print "NOTE: No wig file for chr$chr found. The chromosome will be skipped\n" if $valid==0;
}
  
    if (scalar @sample_chromosomes > 0){
        
    chdir $working_dir;
   my $out_sample_dir="$working_dir/$cur_wig_subfolder"."_OUTPUT";
if (  -d "$out_sample_dir" ){
 deldir($out_sample_dir);
    print "$out_sample_dir will be overwritten! \n";
    
    
}
mkdir $out_sample_dir;
print "Splitting $path_to_peak_file into individual chromosomes...";
pickOutChromosomes($path_to_peak_file, $out_sample_dir);
print "Done\n";
print "Checking $path_to_peak_file for overlapping peaks... ";
my $contains_overlapping=checkOverlappingPeaks($path_to_peak_file);
if($contains_overlapping){
    $has_overlapping_peaks="yes";
    print "found overlapping peaks "
    
}
else{
    $has_overlapping_peaks="no";
}
print "\n";
print "Scheduling parallel processing jobs and running MinMedianMax on $cur_wig_subfolder...\n";

$all_chroms=getNewAllChroms(@sample_chromosomes) if scalar @sample_chromosomes < scalar @chromosomes;
$all_chroms=$all_chroms_old if scalar @sample_chromosomes == scalar @chromosomes;

system("perl $SCRIPTS_PATH/processManager.pl $all_chroms $cur_wig_subfolder $path_to_wig_files $out_sample_dir $SCRIPTS_PATH $num_cores $has_overlapping_peaks");


`cat $out_sample_dir/*_OUT > $out_sample_dir/OUTPUT_ALL_CHR.txt`;
#$result2=`rm $out_sample_dir/chr*.txt`;








    
    
    }
    
}
@time_data=localtime(time);
print "-----------Main program finished on: ".(localtime)."-------------\n";

sub checkOverlappingPeaks{
    
    my $peak_file=$_[0];
    
    my $num_un_merged=  `cat $peak_file | awk '{printf "%s\t%d\t%d\\n", \$1, \$2, \$3}' | sort -u  |wc |awk '{print \$1}'`;
    chomp $num_un_merged;
    
    my $num_merged=`cat $peak_file | awk '{printf "%s\t%d\t%d\\n", \$1, \$2, \$3}' | sort -u |mergeBed -i stdin | wc|awk '{print \$1}'`;
    chomp $num_merged;
    
return 1    if($num_merged< $num_un_merged); # contains overlapping peaks

return 0;
}
sub getNewAllChroms{
    my $all_chroms;
    
    my @new_chroms=@_;
    
    if (scalar @new_chroms == 1){
        
        $all_chroms="$new_chroms[0]";
    }
    else{
       $all_chroms=$new_chroms[0];
    
    for( my $i=1;$i<=$#new_chroms;$i++){
        $all_chroms=$all_chroms.",$new_chroms[$i]";
        
        
    }
        
        
    }
    
    
    return $all_chroms;
}


sub pickOutChromosomes{
    my $main_peak_file=$_[0];
    my $output_dir=$_[1];
    
    
    open(PEAK_FILE, "<$main_peak_file");
    
    while(my $line=<PEAK_FILE>){
        if($line!~m/chr[^\t\s]+[\t\s]+\d+[\t\s]+\d+.*/){
            
            die ("The peak file is improperly formatted. The peak file must be in this format: chromosome< TAB or space>< start coord><TAB or space><stop coord>... . The program will exit.\n")
            
        }
        
    }
    
    foreach(@sample_chromosomes){
       my $chr=$_;
      
    system("awk '\$1==\"chr$chr\"' $main_peak_file > $output_dir/chr$chr.txt") ;#or die "Can't extract chr$chr peaks from $main_peak_file\n";
       system("sort  $output_dir/chr$chr.txt -n -k 2 > $output_dir/tmp$chr.txt");
       system("mv $output_dir/tmp$chr.txt $output_dir/chr$chr.txt");
        
        
    }
    
    
    
    
    
}
sub getParams{
    
if($#_<1){die "Usage: ./getMinMedianMax.pl params.txt num_cores\n";}
    open(PARAMS_FILE, "$_[0]") || die "Invalid parameters file\n";
    $params_file=$_[0];
    if ($_[1]<0 or $_[1]>4) {die "You must specify between 1 and 4 cores\n";}
$num_cores=$_[1];


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
    elsif ($line=~m/^sample_dirs:([\-\w\+]+)$/) {
        
        if( ! -d "$path_to_wig_files$1"){
            
            die "sample directory $path_to_wig_files$1 doesn't exist\n";
        }
        push @wig_subfolders, $1;
        
        
    }
    elsif($line=~m/^sample_dirs:(([\-\w\+]+\,)+[\-\w\+]+)$/){
       @wig_subfolders=split(/,/, $1);
       my @new_wig_subfolders=();
       
       for( my $i=0; $i<=$#wig_subfolders; $i++){
        
        if(! -d "$path_to_wig_files$_"){
            
            print "sample $_ doesn't exist and will be ignored\n";
            
            
        }
        else {
            
            push @new_wig_subfolders, $wig_subfolders[$i];
        }
        
       }
       @wig_subfolders=@new_wig_subfolders;
       
       
       
       
        
        
    }
    
    elsif($line!~m/^#/ and $line ne "") {print "Not recognized as a parameter line: $line\n";}
   
    
}
$path_to_wig_files ne ""|| die "You must provide the full path to sample directories (including all forward slashes)\n";
   scalar @wig_subfolders >0 || die "You must specify at least one sample directory (no forward slashes)\n";
   scalar @chromosomes>0 || die "You must specify at least 1 chromosome\n";
   $path_to_peak_file ne "" || die "You must provide the absolute path to a valid regions file\n";
   print "All input parameters seem valid\n";

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