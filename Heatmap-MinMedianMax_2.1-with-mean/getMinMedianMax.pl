#!/usr/bin/perl -w
BEGIN { $^W=0}


 $time{$format};

 use Cwd;
my $working_dir = getcwd; #directory from which the script will be executed
my $SCRIPTS_PATH="$working_dir/Heatmap_Script_Files";
my @chromosomes=();
my @sgr_subfolders=();
my $path_to_peak_file=$path_to_sgr_files=$params_file=$all_chroms="";
my $num_cores=0;






getParams(@ARGV); #populates @chromosomes, @sgr_subfolders, $path_to_peak_file and $path_to_sgr_files




print "-----------Main program started on: ".(localtime)."-------------\n";
my $all_chroms_old=$all_chroms;
foreach(@sgr_subfolders){
    my $cur_sgr_subfolder=$_;
    chdir "$path_to_sgr_files/$cur_sgr_subfolder";
    my @sgr_files=<*.sgr>;
     @sample_chromosomes=@chromosomes;
     
    foreach(@sgr_files){
        my $sgr_file=$_;
       my $i=-1;
       my $invalid_chr=-1;
      
        foreach(@sample_chromosomes){
                $i++;
            my $chr=$_;
            
            if($sgr_file=~/chr$chr\.sgr|chr$chr\_/){ #this is the SGR file we will be working with
            
                if  (! -r $sgr_file){ #if read permissions are set
                     print "NOTE:$sgr_file does not have read permissions. $cur_sgr_subfolder chromosome $chr will be skipped!\n";
                    $invalid_chr=$i;
                    
                    
                }
                else {
                    
                    $invalid_chr=-1;
                }
             
                
                
                
            }
            
            
            
        }
        
        delete $sample_chromosomes[$invalid_chr] if $invalid_chr >=0;
        
    } #done looping through sgr files
  
    if (scalar @sample_chromosomes > 0){
        
    chdir $working_dir;
   my $out_sample_dir="$working_dir/$cur_sgr_subfolder"."_OUTPUT";
if (  -d "$out_sample_dir" ){
 deldir($out_sample_dir);
    print "$out_sample_dir will be overwritten! \n";
    
    
}
mkdir $out_sample_dir;
print "Splitting $path_to_peak_file into individual chromosomes...";
pickOutChromosomes($path_to_peak_file, $out_sample_dir);
print "Done\n";
print "Scheduling parallel processing jobs and running MinMedianMax on $cur_sgr_subfolder...\n";

$all_chroms=getNewAllChroms(@sample_chromosomes) if scalar @sample_chromosomes < scalar @chromosomes;
$all_chroms=$all_chroms_old if scalar @sample_chromosomes == scalar @chromosomes;

system("perl $SCRIPTS_PATH/processManager.pl $all_chroms $cur_sgr_subfolder $path_to_sgr_files $out_sample_dir $SCRIPTS_PATH $num_cores");


`cat $out_sample_dir/*_OUT > $out_sample_dir/OUTPUT_ALL_CHR.txt`;
#$result2=`rm $out_sample_dir/chr*.txt`;








    
    
    }
    
}
@time_data=localtime(time);
print "-----------Main program finished on: ".(localtime)."-------------\n";

sub getNewAllChroms{
    my $all_chroms;
    
    my @new_chroms=@_;
    
    if (scalar @new_chroms == 1){
        
        $all_chroms="$new_chroms[0]";
    }
    else{
        my $last_chrom = pop @new_chroms;
    
    foreach(@new_chroms){
        $all_chroms=$all_chroms.",$_";
        
        
    }
    $all_chroms="$all_chroms,$last_chrom";
        
        
    }
    
    
    return $all_chroms;
}


sub pickOutChromosomes{
    my $main_peak_file=$_[0];
    my $output_dir=$_[1];
    
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
            
            die "sgr directory $1 doesn't exist\n";
        }
        
        my @path_to_arr=split(//, $1);
        $path_to_sgr_files="$1/";
      $path_to_sgr_files=$1  if ($path_to_arr[$#path_to_arr] eq "/");
        
    
        
        
    }
    elsif ($line=~m/^sample_dirs:([\-\w]+)$/) {
        
        if( ! -d "$path_to_sgr_files$1"){
            
            die "sample directory $path_to_sgr_files$1 doesn't exist\n";
        }
        push @sgr_subfolders, $1;
        
        
    }
    elsif($line=~m/^sample_dirs:(([\-\w]+\,)+[\-\w]+)$/){
       @sgr_subfolders=split(/,/, $1);
       my @new_sgr_subfolders=();
       
       for( my $i=0; $i<=$#sgr_subfolders; $i++){
        
        if(! -d "$path_to_sgr_files$_"){
            
            print "sample $_ doesn't exist and will be ignored\n";
            
            
        }
        else {
            
            push @new_sgr_subfolders, $sgr_subfolders[$i];
        }
        
       }
       @sgr_subfolders=@new_sgr_subfolders;
       
       
       
       
        
        
    }
    else {print "Not recognized as a parameter line: $line\n";}
   
    
}
$path_to_sgr_files ne ""|| die "You must provide the full path to sample directories (including all forward slashes)\n";
   scalar @sgr_subfolders >0 || die "You must specify at least one sample directory (no forward slashes)\n";
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