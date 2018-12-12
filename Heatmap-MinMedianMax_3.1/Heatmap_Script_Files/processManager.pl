#!/usr/local/roadm/bin/perl
BEGIN { $^W = 0 }
#Author: Alina Saiakhova, Scacheri Lab, Case Western Reserve University

use POSIX;
use List::Util qw[min max];
use Data::Dumper;



my $all_chroms=shift @ARGV;
my $cur_wig_subfolder=shift @ARGV;

my $path_to_wig_files=shift @ARGV;
my $output_dir=shift @ARGV;
my $SCRIPTS_PATH=shift @ARGV;
my $num_cores=shift @ARGV;
my $has_overlapping=shift @ARGV;
my @children;
my %processors;
my %sizes;
my %file_sizes;
my @chroms_to_run;
my $chr_count=0;

chdir $output_dir;

#get chromosome names
if($all_chroms =~/,/){
        
     @chroms_to_run=split(/,/, $all_chroms);
}
else{
        
        push @chroms_to_run, $all_chroms;
}

$chr_count=scalar @chroms_to_run;


my $num_processors=min($num_cores,4, $#chroms_to_run+1);

foreach(@chroms_to_run){
       my( $chrom)=$_;
        $file_sizes{$chrom}= -s "$output_dir/chr$chrom.txt"; #get file sizes
       
        if (!defined($file_sizes{$chrom})){
                
                $file_sizes{$chrom}=100;
        }
        elsif($file_sizes{$chrom}==0){delete $file_sizes{$chrom}; print "WARNING: The peak file for chr$chrom appears to be empty. The chromosome will not be processed\n";}

        
        
}


#sort files_to_run;
foreach $key (sort {$file_sizes{$b} <=> $file_sizes{$a}} (keys(%file_sizes))){
        push @sorted_chrom_keys, $key;
        push @sorted_chrom_values, $file_sizes{$key};
}

#schedule processors
for(my $i=1; $i<=$num_processors; $i++){
        $sizes{$i}=$sorted_chrom_values[$i-1];
        $index=$i-1;
      push @{  $processors{$i}}, $sorted_chrom_keys[$index];      
}
for(my $j=1; $j<=$num_processors; $j++){
         shift @sorted_chrom_keys;
      shift @sorted_chrom_values;
        
}

$max_size=$chr_count-$num_processors;

my $counter=1;
while($counter<=$max_size){  
    my @keys=();
      
        foreach $key (sort {$sizes{$a} <=> $sizes{$b} }(keys(%sizes))){
                push @keys, $key;
                
                }
      
      push @{  $processors{$keys[0]}}, $sorted_chrom_keys[0];
    
      
      
      $sizes{$keys[0]}=$sizes{$keys[0]}+$sorted_chrom_values[0];
   
     
    
       shift @sorted_chrom_keys;
      shift @sorted_chrom_values;
  
 $counter++;
}

my @chrom_sets=();


for(my $chrom_i=0; $chrom_i<=$max_size; $chrom_i++){

for ( my $count = 1; $count <= $num_processors; $count++) {
     
     
     if(defined($processors{$count}[$chrom_i])){
          
         
         $chrom_sets[$chrom_i]=$chrom_sets[$chrom_i]." ".$processors{$count}[$chrom_i];
     }
    
     
     
}
}
#start program


foreach(@chrom_sets){
    my $chrom_set=$_;
   
    
    my  @chroms=split(/\s+/, $chrom_set);
     shift @chroms;
   
    
#start program
for ( my $count = 1; $count <= min($#chroms+1,$num_processors); $count++) {
      
       my $pid = fork();
       
        if ($pid) {
        # parent
        
        push(@children, $pid);
        } elsif ($pid == 0) {
                # child
                if(defined($chroms[$count-1])){
               
               @args= ($count, $cur_wig_subfolder, $path_to_wig_files, $output_dir, $chroms[$count-1]); #array of chroms assigned to pid, pid, cur-wig-subfolder, path_to_wig_files, output directory
               runMinMedianMaxNonOverlapping(@args) if $has_overlapping eq "no" ;
               runMinMedianMaxOverlapping(@args) if $has_overlapping eq "yes" ;
                
                exit 0;
                }else {exit 0}
        } else {
                print "couldnt fork: $!\n";
        }



}

foreach (@children) {
        my $tmp = waitpid($_, 0);

}       


}
sub runMinMedianMaxOverlapping {
my $p_num=$_[0];
       @ARGV=(@_);
      eval { require "$SCRIPTS_PATH/MedianMaxMin_overlapping.pl" };
      
     # return $p_num;

}

sub runMinMedianMaxNonOverlapping {
my $p_num=$_[0];
       @ARGV=(@_);
      eval { require "$SCRIPTS_PATH/MedianMaxMin_non-overlapping.pl" };
      
     # return $p_num;

}
