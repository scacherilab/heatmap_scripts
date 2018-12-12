#!/usr/bin/perl -w
BEGIN { $^W = 0 }
#Author: Alina Saiakhova, Scacheri Lab, Case Western Reserve University

use Data::Dumper;


my $pnum=shift @ARGV;
my $cur_wig_subfolder=shift @ARGV;
my $path_to_wig_files=shift @ARGV;
my $output_dir=shift @ARGV;
my $num_windows=shift @ARGV;
my $dist_each_dir=shift @ARGV;
my @chroms=@ARGV;
my $full_path_to_sample_wigs="$path_to_wig_files$cur_wig_subfolder/";


foreach(@chroms){
my $wig_open=1;
my $peakfile_open=1;  
my ($chr)=$_;
print "Process $pnum: Starting chromosome $chr...\n";
my $wig_file=validateChromWig($full_path_to_sample_wigs, $chr);

open(CHR_WIG_FILE,"<$full_path_to_sample_wigs/$wig_file") ||  (print "Can't open $full_path_to_sample_wigs/$wig_file"  and $wig_open=0);
open(CHR_PEAK_FILE, "<$output_dir/chr$chr.txt") ||  (print "Can't open $output_dir/chr$chr.txt" and $peakfile_open=0);



my $flag=0;
$num_wig_steps=0;
my %wig_signal_vals=%peak_windows=();
my $peak=0;

my ($success, $type_wig, $wig_step, $wig_start_coord)=getWigFileInfo($wig_file, $chr); #checks to make sure wig file is properly formatted

if($wig_open &&  $peakfile_open && $success=="true"){ #proceed only if both the wig and the peak files are open
  
if ($pnum==1){open(OUTFILE1, ">$output_dir/chr".$chr."_OUT") || print ("Can't open peaks file: $output_dir/chr$chr.txt\n")}
elsif($pnum==2) {open(OUTFILE2, ">$output_dir/chr".$chr."_OUT")|| print ("Can't open peaks file: $output_dir/chr$chr.txt\n")}
elsif($pnum==3) {open(OUTFILE3, ">$output_dir/chr".$chr."_OUT")|| print ("Can't open peaks file: $output_dir/chr$chr.txt\n")}
elsif($pnum==4) {open(OUTFILE4, ">$output_dir/chr".$chr."_OUT")|| print ("Can't open peaks file: $output_dir/chr$chr.txt\n")}

my ($peaks_str)=populatePeakList(CHR_PEAK_FILE); #store peaks
my @peaks=@$peaks_str;
my $window_size=$dist_each_dir*2/$num_windows;

while(my $wig_line=<CHR_WIG_FILE>){
 chomp $wig_line;
 
 if($wig_line =~m/^\d+/){ #if data line, then proceed (must start with integer)
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

$num_p_window=1;

if($flag>0){splice @peaks, 0, $flag;$flag=0} #removes peaks that occur before current coordinate analyzed


foreach(@peaks){
   
my $peak=$_;

$peak_windows{$peak}='1' if !defined($peak_windows{$peak}) or !exists($peak_windows{$peak});
$num_p_window=$peak_windows{$peak};
$peak_start=$peak-$dist_each_dir+(($num_p_window-1)*$window_size);
$peak_end=$peak-$dist_each_dir+(($num_p_window)*$window_size);


if(eof){ #if end of wig file is reached
  
  
  push @{$wig_signal_vals{$peak}}, $wig_signal_val if $wig_coord < $peak_end and $wig_coord>=$peak_start;
  calcMediansInWindow($wig_signal_vals{$peak},$peak,$peak_start, $peak_end) if exists($wig_signal_vals{$peak});
  shift @peaks if scalar @{$wig_signal_vals{$peak}}>0; ##new

  undef $wig_signal_vals{$peak};
  $peak_end2=$peak_end;
  $peak_start2=$peak_start;

while($peak_windows{$peak}<$num_windows){
  
    calcMediansInWindow("NA", $peak, $peak_start2, $peak_end2);
    $peak_windows{$peak}++;
    $peak_end2=$peak-$dist_each_dir+(($peak_windows{$peak})*$window_size);
   $peak_start2=$peak-$dist_each_dir+(($peak_windows{$peak}-1)*$window_size);
   
}

$num_p_window2=$peak_windows{$peak};

if($num_p_window2>=$num_windows)
{
      if($wig_coord>=$peak_end2 and $num_p_window2==$num_windows){$final_end=$peak-$dist_each_dir+(($num_p_window2)*$window_size); $final_start=$peak-$dist_each_dir+(($num_p_window2-1)*$window_size);calcMediansInWindow("NA", $peak, $final_start, $final_end);}
$flag++;
delete $wig_signal_vals{$peak};
delete $peak_windows{$peak};
}
   
} #end if eof

elsif($wig_coord< $peak_end and $wig_coord>= $peak_start){ #if wig coord falls within given peak window, store wig signal
      
   
    push @{ $wig_signal_vals{$peak} }, $wig_signal_val;
   
  
}elsif($wig_coord>=$peak_end ){    ########## coordinate exceeds allowed range
   
   calcMediansInWindow($wig_signal_vals{$peak}, $peak, $peak_start, $peak_end) if exists($wig_signal_vals{$peak});
   $num_p_window=$peak_windows{$peak};
   undef $wig_signal_vals{$peak};
   $peak_end2=$peak_end;
   $peak_start2=$peak_start;

while($wig_coord>=$peak_end2 and $peak_windows{$peak}<$num_windows){

    calcMediansInWindow("NA", $peak, $peak_start2, $peak_end2);
    $peak_windows{$peak}++;
    $peak_end2=$peak-$dist_each_dir+(($peak_windows{$peak})*$window_size);
   $peak_start2=$peak-$dist_each_dir+(($peak_windows{$peak}-1)*$window_size);
    
   
}

if($peak_windows{$peak}<=$num_windows and $wig_coord<$peak_end2){

        push @{ $wig_signal_vals{$peak} }, $wig_signal_val;
       
        
    }else{ #begin else
   
   if($wig_coord>=$peak_end2 and $peak_windows{$peak}==$num_windows){$final_end=$peak-$dist_each_dir+(($peak_windows{$peak})*$window_size);$final_start=$peak-$dist_each_dir+(($peak_windows{$peak}-1)*$window_size);calcMediansInWindow("NA",$peak, $final_start, $final_end);}
   $flag++;
   delete $wig_signal_vals{$peak};
   delete $peak_windows{$peak};
    } #end else
    
    
    } else{ #wig coordinate is smaller than any of the coordinates in the current window of hte current peak. We can skip all subsequent peaks because that will also be true for those peaks. 
   
   
   last;}

} ##finish looping through peaks

$num_wig_steps++; 
 } else{ #if not data line, don't do anything
   
   
 }
 
} #finish looping through wig file


print "Process $pnum: Done with chromosome $chr...\n";

  foreach(@peaks){
  
   for($k=1; $k<=$num_windows; $k++){
       if($pnum == 1){print OUTFILE1 "$_\tNA\n";}
    elsif($pnum == 2){print OUTFILE2 "$_\tNA\n";}
    elsif($pnum == 3){print OUTFILE3 "$_\tNA\n";}
    elsif($pnum==4){print OUTFILE4 "$_\tNA\n";}
   }
  }
close(CHR_WIG_FILE);
close(CHR_PEAK_FILE);
    if($pnum == 1){close OUTFILE1;}
    elsif($pnum == 2){ close OUTFILE2}
    elsif($pnum == 3){close OUTFILE3;}
    elsif($pnum==4){close OUTFILE4;}
    


undef %wig_signal_vals;
undef %peak_windows;

} #finish checking if both wig and peak file are open
} #finish looping through chromosomes


######HELPER FUNCTIONS######
sub calcMediansInWindow {
   

    my @pole =@{$_[0]};
  
     my $ret;
    
    if($#pole>=0 and $pole[0] ne "NA"){
      
      
    

   

  @pole=  sort { $a <=> $b } @pole;
 

    if( (@pole % 2) == 1 ) {
        $ret = $pole[((@pole+1) / 2)-1];
       
    } else {
        $ret = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
       
        
    }}
    else{
      $ret="NA";
    }
    
   
    if($pnum == 1){print OUTFILE1 "$_[1]  $ret range: $_[2] and $_[3] \n";}
    elsif($pnum == 2){print OUTFILE2 "$_[1]  $ret range: $_[2] and $_[3] \n";}
    elsif($pnum == 3){print OUTFILE3 "$_[1]  $ret range: $_[2] and $_[3] \n";}
    elsif($pnum==4){print OUTFILE4 "$_[1]  $ret range: $_[2] and $_[3] \n";}
    
    


    
    
    
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
     if(`echo \$(wc -l $full_path_to_sample_wigs$wig_file |awk '{print \$1}')`==0){
          print "WARNING:$wig_file appears to be empty\n"
     }
     
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
     my $success="";
     
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
sub populatePeakList{
     my $CHR_PEAK_FILE=$_[0];
     my @peak_starts=@peak_ends=();
     my %existing_peaks=();
          while(<$CHR_PEAK_FILE>){
        my($peakline)=$_;
        my @columns=split (/\t|\s+/, $peakline);
        if(!$existing_peaks{"$columns[1]"}++){
         push @peaks, $columns[1];
         
        }
        
    }
     return(\@peaks);
}
