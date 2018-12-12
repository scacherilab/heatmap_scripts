#!/usr/bin/perl -w
BEGIN { $^W = 0 }
#Author: Alina Saiakhova, Scacheri Lab, Case Western Reserve University
use Data::Dumper;
use List::Util qw(sum reduce);

my %peaks;
my %medians;
my %max;
my %min;



my $num=shift @ARGV;
my $cur_sgr_subfolder=shift @ARGV;
my $path_to_sgr_files=shift @ARGV;
my $output_dir=shift @ARGV;
my @chroms=@ARGV;
my $full_path_to_sample_sgr="$path_to_sgr_files$cur_sgr_subfolder";
my $sgr_open=1;
my $peakfile_open=1;




 my $flag=0;
  
 

my $start=();
my $end=();

foreach(@chroms){

   my ($chr)=$_;
   print "Process $num: Starting chromosome $chr...\n";
my $sgr_file=`ls $full_path_to_sample_sgr/| grep 'chr$chr\[.\_]' | head -1 `;


die "Unrecognized SGR file format" if $sgr_file == -1 ;
if ($sgr_file eq ""){
print "SGR file for chromosome$chr doesn't exist. Skipping chromosome $chr \n";
$sgr_open=0;

 }
chomp $sgr_file;



open(INPUTFILE,"<$full_path_to_sample_sgr/$sgr_file") || handleFileOpenError("sgr",$cur_sgr_subfolder, $chr);
open(PEAKFILE, "<$output_dir/chr$chr.txt") ||  handleFileOpenError("peak",$cur_sgr_subfolder, $chr);
open(SKIPPED, ">$output_dir/skipped_duplicate_peaks.txt");
if($sgr_open &&  $peakfile_open){
  
if ($num==1){open(OUTFILE1, ">$output_dir/chr".$chr."_OUT") || print ("Can't open peaks file: $output_dir/chr$chr.txt\n")}
elsif($num==2) {open(OUTFILE2, ">$output_dir/chr".$chr."_OUT")|| print ("Can't open peaks file: $output_dir/chr$chr.txt\n")}
elsif($num==3) {open(OUTFILE3, ">$output_dir/chr".$chr."_OUT")|| print ("Can't open peaks file: $output_dir/chr$chr.txt\n")}
elsif($num==4) {open(OUTFILE4, ">$output_dir/chr".$chr."_OUT")|| print ("Can't open peaks file: $output_dir/chr$chr.txt\n")}
my %existing_peaks=();
     while(<PEAKFILE>){
        my($peakline)=$_;
        my @columns=split (/\t|\s+/, $peakline);
        if(!defined($existing_peaks{"$columns[1]"."_"."$columns[2]"})){
         push @start, $columns[1];
push @end, $columns[2];
$existing_peaks{"$columns[2]"."_"."$columns[2]"}++;
         
        }
        else {
         
         print SKIPPED "$columns[0]\t$columns[1]\t$columns[2]\n";
         
         
        }
        
        



        
    }

while(my $line=<INPUTFILE>){
  
    chomp $line;
   
    @location=split (/\t|\s+/, $line);
shift(@location);
$sgr_coord=$location[0];
$peak_height=$location[1];
if($peak_height=~ m/E/){
   @argz=($peak_height);
   $peak_height=convertToInt(@argz);
}

 if($flag>0){
   
   splice @start, 0,$flag;
   splice @end, 0, $flag;
   $flag=0;
   
   }



for(my $i=0; $i<=$#start;$i++){
   $peak="$start[$i],$end[$i]";

 $start=$start[$i];
 $end=$end[$i];
 if(eof){
   
   
      @args= ($medians{$peak}, $peak, $num, $chr);
        calcMeanMedianMaxMin(@args);
     if($start<=$sgr_coord and $end >=$sgr_coord){
        push @{$medians{$peak}}, $peak_height;
       
        
     }
        
   
 }
     elsif($start<=$sgr_coord and $end >=$sgr_coord){
        push @{$medians{$peak}}, $peak_height;
        
     }
     elsif( $sgr_coord>$end  ){
     
        @args= ($medians{$peak}, $peak, $num, $chr);
        calcMeanMedianMaxMin(@args);
 $flag++;
       delete $medians{$peak};
       
      
     }
     elsif($sgr_coord<$start){
   last;
   

      }
     
   
   
     
    
  }



}
print "Process $num: Done with chromosome $chr...\n";

    

}
}


sub handleFileOpenError{
   
   my $type=$_[0];
   my $sample=$_[1];
   my $chr=$_[2];
   
   if($type eq "sgr"){
      $sgr_open=0;
      
   }
   elsif($type eq "peak"){
      
      $file_open=0;
   }
   print ("Skipping  $sample"."'s chromosome $chr\n");
}
 sub findNewPeakForLoc{
    my $sgr_coord=$_[0];
  for(my $i=0; $i<=$#start;$i++){
 my  $peak="$start[$i],$end[$i]";

 $start=$start[$i];
 $end=$end[$i];
     if($start<=$sgr_coord and $end >=$sgr_coord){
   
        push @{$medians{$peak}}, $peak_height;
        last;
       
    }
    
     }
    
 }
sub convertToInt{
   
   my $int=$_[0];
 
   $int =~ tr/.//d;
   $int =~tr/-//d;
   ($before, $after)=split(/E/,$int);
   $new_int=$before;
  for( $i=1; $i<$after;$i++){
 
  $new_int="0$new_int";
  }
  $new_int="0.$new_int";
  return $new_int
}

sub calcMeanMedianMaxMin {
   
  my  $peak2=$_[1] ;
    $num=$_[2];
    ($begin, $finish)=split(/,/, $peak2);

    my @pole =@{$_[0]};
    
    

    my $median;

  @pole=  sort { $a <=> $b }  @pole;
     if($#pole>=0 ){

    if( (@pole % 2) == 1 ) {
        $median = $pole[((@pole+1) / 2)-1];
        
    } else {
        $median = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
        
        
    }
    my @pole2=@pole;
    $mean= (reduce{ $a + $b }  @pole)/scalar @pole;
    $min= shift (@pole);
    $max= pop(@pole2);
    
    
    
        

    }
     else{
        $min="NA";
        $max="NA";
        $median="NA";
        $mean="NA";
     }
     
     
     
      if($num == 1){print OUTFILE1 "chr$_[3] \t$begin\t$finish\t$min\t$median\t$max\t$mean\n";}
    elsif($num == 2){print OUTFILE2 "chr$_[3] \t$begin\t$finish\t$min\t$median\t$max\t$mean\n";}
    elsif($num == 3){print OUTFILE3 "chr$_[3] \t$begin\t$finish\t$min\t$median\t$max\t$mean\n";}
    elsif($num==4){print OUTFILE4 "chr$_[3] \t$begin\t$finish\t$min\t$median\t$max\t$mean\n";}
    
    

}

