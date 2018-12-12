#!/usr/bin/perl -w
BEGIN { $^W = 0 }
#Author: Alina Saiakhova, Scacheri Lab, Case Western Reserve University
use Data::Dumper;

my %peaks;
my %wig_signal_vals;
my %max;
my %min;



my $pnum=shift @ARGV;
my $cur_wig_dir=shift @ARGV;
my $path_to_wigs=shift @ARGV;
my $output_dir=shift @ARGV;
my @chroms=@ARGV;
my $full_path_to_sample_wigs="$path_to_wigs$cur_wig_dir/";
my $wig_open=1;
my $peakfile_open=1;


foreach(@chroms){

my ($chr)=$_;
print "Process $pnum: Starting chromosome $chr...\n";
my $wig_file=validateChromWig($full_path_to_sample_wigs, $chr); 

open(CHR_WIG_FILE,"<$full_path_to_sample_wigs/$wig_file") || ($wig_open=0 && print("Can't open $full_path_to_sample_wigs/$wig_file"));
open(CHR_PEAK_FILE, "<$output_dir/chr$chr.txt") ||  ($peakfile_open=0 && print("Can't open $output_dir/chr$chr.txt"));

my $flag=0;
$num_wig_steps=0;

my ($success, $type_wig, $wig_step, $wig_start_coord)=getWigFileInfo($wig_file, $chr); #checks to make sure wig file is properly formatted

if($wig_open &&  $peakfile_open && $success=="true"){ #proceed only if both the wig and the peak files are open
if ($pnum==1){open(OUTFILE1, ">$output_dir/chr".$chr."_OUT") || print ("Can't open peaks file: $output_dir/chr$chr.txt\n")}
elsif($pnum==2) {open(OUTFILE2, ">$output_dir/chr".$chr."_OUT")|| print ("Can't open peaks file: $output_dir/chr$chr.txt\n")}
elsif($pnum==3) {open(OUTFILE3, ">$output_dir/chr".$chr."_OUT")|| print ("Can't open peaks file: $output_dir/chr$chr.txt\n")}
elsif($pnum==4) {open(OUTFILE4, ">$output_dir/chr".$chr."_OUT")|| print ("Can't open peaks file: $output_dir/chr$chr.txt\n")}

my ($peak_starts_str, $peak_ends_str)=populatePeakList(CHR_PEAK_FILE); #store peaks
my @peak_starts=@$peak_starts_str;
my @peak_ends=@$peak_ends_str;

while(my $wig_line=<CHR_WIG_FILE>){
chomp $wig_line;

if($wig_line =~m/^\d+/){ #if data line, then proceed

if($type_wig eq "fixedStep"){
     $wig_coord=$wig_start_coord+$num_wig_steps*$wig_step-1;
     $wig_signal_val=$wig_line;
     
}
else{ #variable step
     my @wig_liner=split (/\t|\s+/, $wig_line);
     $wig_coord=$wig_liner[0]-1;
     $wig_signal_val=$wig_liner[1];
     $wig_signal_val=convertToInt($wig_signal_val) if($wig_signal_val=~m/E/);
     
}
     
my $k=0;
      while(!defined($peak_starts[$k]) and $k<=$#peak_starts){ $k++;}
if($k>0){splice @peak_starts,  0, $k-1;splice @peak_ends, 0, $k-1;}


for(my $i=0; $i<=$#peak_starts;$i++){ #start looping through peaks at this wig coord
$peak="$peak_starts[$i],$peak_ends[$i]";

if (defined($peak_starts[$i])){
$peak_start=$peak_starts[$i];
$peak_end=$peak_ends[$i];


 if(eof){ #if end of file is reached
   calcMedianMaxMin($wig_signal_vals{$peak}, $peak, $pnum, $chr);
   push @{$wig_signal_vals{$peak}}, $wig_signal_val if $peak_start<=$wig_coord and $peak_end>=$wig_coord;
       
}
 elsif($peak_start<=$wig_coord and $peak_end >=$wig_coord){ #if wig coordinate falls within peak
        push @{$wig_signal_vals{$peak}}, $wig_signal_val;
      
        
}
elsif( $wig_coord>$peak_end  ){ #if wig coordinate is greater than the peak end coordinate
    
        calcMedianMaxMin($wig_signal_vals{$peak},$peak,$pnum,$chr);
        undef $peak_starts[$i];
        undef $peak_ends[$i];
        undef $wig_signal_vals{$peak}; #and free memory occupied by its wig_signal_vals
       
      
     }
elsif($wig_coord<$peak_start){ #if wig coordinate occurs before the current peak, no need to proceed with the rest of the peak list
          
   last;
   
}
else{
     
     #Not sure what the other cases would be right now
}
     
} #finish checking if peak start is undefined 
} #finish looping through peaks at this wig coord
    
  


$num_wig_steps++; #this only really applies to a fixedStep wig file
}
elsif ($wig_line=~/fixedStep|variableStep/){
     
     ($success, $type_wig, $wig_step, $wig_start_coord)=getWigLineInfo($wig_line, $chr);
     $num_wig_steps=0;
     
     if($success eq "false"){
          
          "Skipping the rest of chr$chr\n";
          last;
          seek(CHR_WIG_FILE, 0, 2); 
          
     }
}   else{}#if it's not a data line, then don't do anything

} #finish looping through wig file

print "Process $pnum: Done with chromosome $chr...\n"; #We are done with this chromosome
close(CHR_WIG_FILE);

}
else{ #Will end up in here if: wig file is malformed, can't open wig file, or can't open peak file
    
     }
     
} #finish looping through all the chromosomes in assigned to this process

######HELPER FUNCTIONS######
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

sub calcMedianMaxMin {
   
   my $peak2=$_[1];
  my  $num=$_[2];
 
    my ($begin, $finish)=split(/,/, $peak2);

    my @pole =@{$_[0]};
    
    

    my $ret;

  @pole=  sort { $a <=> $b }  @pole;
     if($#pole>=0 ){

    if( (@pole % 2) == 1 ) {
        $ret = $pole[((@pole+1) / 2)-1];
        
    } else {
        $ret = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
        
        
    }
    my @pole2=@pole;
    $min= shift (@pole);
    $max= pop(@pole2);
    
    
        

    }
     else{
        $min="NA";
        $max="NA";
        $ret="NA";
     }
     
     
     
      if($num == 1){print OUTFILE1 "chr$_[3] \t$begin\t$finish\t$min\t$ret\t$max\n";}
    elsif($num == 2){print OUTFILE2 "chr$_[3] \t$begin\t$finish\t$min\t$ret\t$max\n";}
    elsif($num == 3){print OUTFILE3 "chr$_[3] \t$begin\t$finish\t$min\t$ret\t$max\n";}
    elsif($num==4){print OUTFILE4 "chr$_[3] \t$begin\t$finish\t$min\t$ret\t$max\n";}
    elsif($num==5){print OUTFILE5 "chr$_[3] \t$begin\t$finish\t$min\t$ret\t$max\n";}
    elsif($num==6){print OUTFILE6 "chr$_[3] \t$begin\t$finish\t$min\t$ret\t$max\n";}
    elsif($num==7){print OUTFILE7 "chr$_[3] \t$begin\t$finish\t$min\t$ret\t$max\n";}
    

}
sub populatePeakList{
     my $CHR_PEAK_FILE=$_[0];
     my @peak_starts=@peak_ends=();
     my %existing_peaks=();
          while(<$CHR_PEAK_FILE>){
        my($peakline)=$_;
        my @columns=split (/\t|\s+/, $peakline);
        if(!$existing_peaks{"$columns[1]"."_"."$columns[2]"}++){
         push @peak_starts, $columns[1];
         push @peak_ends, $columns[2];
         
        }
        
    }
     return(\@peak_starts, \@peak_ends);
}
sub validateChromWig{
     my $full_path_to_sample_wigs=$_[0];
     my $chr=$_[1];
     
  my $wig_file=   `ls $full_path_to_sample_wigs| grep -E 'chr$chr\[.\_].*wig' | head -1`;
  if($wig_file eq ""){
  print "Could not find any chr$chr wig files in $full_path_to_sample_wigs in this format: *chr$chr.wig* or *chr$chr"."_*.wig. Skipping chr$chr\n";
  return 0;

  }
  else{
     chomp $wig_file;
    # if(`echo \$(wc -l $full_path_to_sample_wigs$wig_file |awk '{print \$1}')`==0){
         # print "WARNING:$wig_file appears to be empty\n"
    # }
     
     return $wig_file;
  }
     
}

sub getWigFileInfo{
     my $wig_file=$_[0];
     my $chr=$_[1];
     my $wig_type="";
     my $step="";
     my $start="";
     my $chrom="";
     
     my $wig_header_line=`grep -E 'fixedStep|variableStep' $full_path_to_sample_wigs$wig_file`;
     if($wig_header_line=~m/chrom=(\S+)\s/){ #get chromosome info
          if($1=="chr$chr"){
               $chrom=$1;
          }
         }
     else{
          $success="false";
          print "$wig_file does not contain chr$chr information. Chromosome $chr will be skipped\n";
          
     }
     
     if($wig_header_line=~m/fixedStep/){
          
         $wig_type="fixedStep";
         if($wig_header_line=~m/step=(\d+)\s/){
          
          $step=$1;
          
         }
         if($wig_header_line=~m/start=(\d+)\s/){
          
          $start=$1;
         }
         
         
         
     }
     elsif($wig_header_line=~m/variableStep/){
          
          $wig_type="variableStep";
          $step="NA";
          $start="NA";
     }
     else{
          
      print    "$wig_file is  malformed. Chromosome $chr will be skipped\n";
          $success="false";
     }
     
     if($wig_type ne "" and $step ne "" and $start ne "" and $chrom ne ""){
          $success="true";
     }
     else{
       print   "$wig_file is  malformed. Chromosome $chr will be skipped\n";
          $success="false";
     }
     
     
     return ($success, $wig_type,$step,$start);
}

sub getWigLineInfo{
     my $wig_header_line=$_[0];
     my $chr=$_[1];
     my $wig_type="";
     my $step="";
     my $start="";
     my $chrom="";
     
     if($wig_header_line=~m/chrom=(\S+)\s/){ #get chromosome info
          if($1=="chr$chr"){
               $chrom=$1;
          }
         }
     else{
          $success="false";
          print "$wig_header_line does not contain chr$chr information. Chromosome $chr will be skipped\n";
          
     }
     
     if($wig_header_line =~ /fixedStep/){
         $wig_type="fixedStep";
         if($wig_header_line=~/step=(\d+)/){
          
          $step=$1;
          
         }
         if($wig_header_line=~/start=(\d+)\s/){
          
          $start=$1;
         }
         
         
         
     }
     elsif($wig_header_line =~ /variableStep/){
          
          $wig_type="variableStep";
          $step="NA";
          $start="NA";
     }
     else{
          
      print   "$wig_header_line is  malformed. Chromosome $chr will be skipped\n";
          $success="false";
     }
     
     if($wig_type ne "" and $step ne "" and $start ne "" and $chrom ne ""){
          $success="true";
     }
     else{
       print   "$wig_header_line is  malformed. Chromosome $chr will be skipped\n";
      
          $success="false";
     }
     
     
     return ($success, $wig_type,$step,$start);
}
