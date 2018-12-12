#!/usr/bin/perl -w
#Author: Alina Saiakhova, Scacheri Lab, Case Western Reserve University
BEGIN { $^W = 0 }
use Data::Dumper;
use Cwd;
my $working_dir = getcwd; #directory from which the script will be executed
my $SCRIPTS_PATH="$working_dir/Helper_Scripts";
my $path_to_chrom_fasta_files="";
my @chromosomes=();
my @sample_chromosomes=();
my $path_to_regions_file=$params_file=$all_chroms=$calc_region_center="";
my $num_cores=$num_windows=$dist_each_dir=0;

print "-----------Main program started on: ".(localtime)."-------------\n";

getParams(@ARGV); #populates @chromosomes, $path_to_regions_file and $path_to_fasta_files







    chdir "$path_to_fasta_files/";
    
 
 foreach(my $i=0; $i<=$#chromosomes;$i++){
                
    my $chr=$chromosomes[$i];
    
    my $fasta_file=`ls $path_to_fasta_files/| grep 'chr$chr.fa' | head -1 `;
    chomp $fasta_file;
    
    if($fasta_file eq "" and $chr ne ""){
        
        print "WARNING: No fasta file (*chr$chr.fa) for chromosome $chr was found in $path_to_fasta_files! The chromosome will be skipped\n";
         $chromosomes[$i]=undef;
    }
    elsif( $fasta_file ne ""){
        
          if  (! -r $fasta_file){ #if read permissions are set
                     print "WARNING:$fasta_file does not have read permissions. The chromosome will be skipped!\n";
                   $chromosomes[$i]=undef;
                    
                    
                }
    }

   }
        
    foreach(@chromosomes){
        
        if (defined($_)){
            
            push @sample_chromosomes, $_;
        }
        
    }
    
  
    if (scalar @sample_chromosomes > 0){
        
    chdir $working_dir;
   my $out_dir="$working_dir/"."OUTPUT";
if (  -d "$out_dir" ){
   deldir($out_dir);
    
    
}
mkdir $out_dir;

if ($calc_region_center eq "yes"){
  print "Calculating regions midpoints and splitting $path_to_regions_file into individual chromosomes...";  
    calcMediansAndPickOutChromosomes($path_to_regions_file, $out_dir);
}
else {
    print "Splitting $path_to_regions_file into individual chromosomes...";
    pickOutChromosomes($path_to_regions_file, $out_dir);
}
print "Done\n";


print "Scheduling parallel processing jobs and running main script...\n";

$all_chroms=getNewAllChroms(@sample_chromosomes) if scalar @sample_chromosomes < scalar @chromosomes;

system("perl $SCRIPTS_PATH/processManager.pl $all_chroms $path_to_fasta_files $out_dir  $num_windows $dist_each_dir $SCRIPTS_PATH $num_cores");


#######Convert to right format########
chdir $out_dir;
my @temp_out_files=<*_OUT>;
foreach(@temp_out_files){
    my $cur_file=$_;
    print "Converting $cur_file to the right format...";
    removeDuplicates($cur_file, $out_dir);
    convertToRightFormat($cur_file, $out_dir, $num_windows);
    print "Done\n";
   
    
}


`cat $out_dir/FINAL*_NO_DUPS > $out_dir/FINAL_ALL_CHR.txt`;
`rm $out_dir/FINAL*_NO_DUPS`;
print "Z-scoring...";
z_score("$out_dir/FINAL_ALL_CHR.txt", $out_dir);
print "Done\n";
chdir $working_dir;
##Clean up###

    
    }
    

@time_data=localtime(time);
print "-----------Main program finished on: ".(localtime)."-------------\n";

sub z_score{
    my $file_to_z_score=$_[0];
    my $output_dir=$_[1];
    open (OUTPUT, ">$output_dir/z_scores_final.txt");
    $max=$num_windows+1;
    
    `R --vanilla $file_to_z_score < $SCRIPTS_PATH/z_score.R`;
    print OUTPUT "chr_coord\t";
    
    for(my $i=1; $i<$num_windows; $i++){
        print OUTPUT "$i\t";
        
    }
    print OUTPUT "$num_windows\n";
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
   my @columns=split (/\t|\s+/, $line);
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
    my $num_dup_regions=0;
    foreach(@sample_chromosomes){
       my $chr=$_;
      
    system("awk '\$1==\"chr$chr\"' $main_peak_file > $output_dir/chr$chr.txt") ;#or die "Can't extract chr$chr peaks from $main_peak_file\n";
    my $before=`wc -l $output_dir/chr$chr.txt| awk '{print \$1}'`;
    
       system("sort  $output_dir/chr$chr.txt -n -k 2 -u > $output_dir/tmp$chr.txt");
       system("mv $output_dir/tmp$chr.txt $output_dir/chr$chr.txt");
     my $after=  `wc -l $output_dir/chr$chr.txt| awk '{print \$1}'`;
        
        my $diff=eval{$before - $after};
        
        $num_dup_regions=$num_dup_regions+$diff;
    }
    
    print "removed $num_dup_regions duplicate regions..." if $num_dup_regions>0;
    
    
    
}

sub calcMediansAndPickOutChromosomes{
    
    my $main_peak_file=$_[0];
    my $output_dir=$_[1];
    open(INPUT, "<$main_peak_file" ) || die "Problem calculating peak centers\n";
    open(OUTPUT, ">$output_dir/main_peak_file.txt") || die "Problem calculating peak centers\n";
   
    while ($line=<INPUT>){
        my @liner=split(/\t|\s+/, $line);
        my $peak_center=eval{($liner[1]+$liner[2])/2};
        print OUTPUT "$liner[0]\t$peak_center\n" ;
        
        
    }
    
    foreach(@sample_chromosomes){
       my $chr=$_;
      
    `awk '\$1==\"chr$chr\"' $output_dir/main_peak_file.txt > $output_dir/chr$chr.txt` ;#or die "Can't extract chr$chr peaks from $main_peak_file\n";
       system("sort  $output_dir/chr$chr.txt -n -k 2 -u > $output_dir/tmp$chr.txt");
       system("mv $output_dir/tmp$chr.txt $output_dir/chr$chr.txt");
      
        
    
    
    
}
     system("rm $output_dir/main_peak_file.txt");
}
sub getParams{
    
 ##read commandline parameters   
if($#_<1){die "Usage: ./getGCWindowedContent params_file.txt num_cores\n";}
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
    elsif($line=~m/^absolute_path_to_regions_file:(.+)$/){
        
        if (! -e $1){
            
            die "regions file $1 doesn't exist\n";
        }
        if(`wc -l $1 | awk '{print \$1}'`==0){
            
            die "Regions file $1 appears to be empty. It may contain MAC-style newline characters.\n";
        }
        
        $path_to_regions_file=$1;
        
    }
    elsif ($line=~m/^absolute_path_to_fasta_files:(.+)$/) {
        
        if( ! -d "$1"){
            
            die "Invalid path to fasta files: $1 \n";
        }
          
       my @path_to_arr=split(//, $1);
        $path_to_fasta_files="$1/";
      $path_to_fasta_files=$1  if ($path_to_arr[$#path_to_arr] eq "/");
        
      
    }
  
    elsif($line=~m/^calculate_region_center\?:(yes|no)$/){
        
        $calc_region_center=$1;
       
        
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
    $path_to_regions_file ne "" || die "You must provide the full path to a valid regions file\n";
    $path_to_fasta_files ne ""|| die "You must provide the full path to the fasta files\n";
    $calc_region_center ne "" || die "You must specify if peak medians will need to be calculated\n";
    $num_windows >0 || die "You must provide number of peak windows\n";
    $dist_each_dir >0 || die "You must specify the distance in each direction of the peak the GC content will need to be calculated for\n";
 

}


sub getNewAllChroms{
    my $all_chroms;
    
    my @new_chroms=@_;
    
  
    
    if (scalar @new_chroms == 1){
        
        $all_chroms="$new_chroms[0]";
    }
    else{
        my $first_chrom = shift @new_chroms;
        $all_chroms=$first_chrom;
    
    foreach(@new_chroms){
        $all_chroms=$all_chroms.",$_";
        
        
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

