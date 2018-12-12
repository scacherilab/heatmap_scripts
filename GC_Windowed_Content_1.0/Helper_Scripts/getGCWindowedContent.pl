#!/usr/bin/perl -w
BEGIN { $^W = 0 }
use Data::Dumper;
#Author: Alina Saiakhova, Scacheri Lab, Case Western Reserve University



my $num=shift @ARGV;
my $path_to_fasta_files=shift @ARGV;
my $output_dir=shift @ARGV;
my $num_windows=shift @ARGV;
my $dist_each_dir=shift @ARGV;
my @chroms=@ARGV;
my $fasta_open=1;
my $regions_open=1;



foreach(@chroms){
   my ($chr)=$_;
   if($chr ne ""){
   
print "Process $num: Starting chromosome $chr...\n";
my $fasta_file=`ls $path_to_fasta_files/| grep 'chr$chr.fa' | head -1 `;
chomp $fasta_file;
if ($fasta_file eq "" and $chr ne ""){
print "FASTA file for chromosome$chr (*chr$chr.fa) doesn't exist. Skipping chromosome $chr. \n";
$fasta_open=0;

 }
else{
open(INPUTFILE,"<$path_to_fasta_files$fasta_file") || handleFileOpenError("fasta", $chr);
}

open(REGIONSFILE, "<$output_dir/chr$chr.txt") ||  handleFileOpenError("regions", $chr);

my @regions;
my %num_g_or_c;
my %num_bases;
my %num_N;
my %windows;


while(<REGIONSFILE>){
my $region_line = $_;
my @liner=split (/\t|\s+/, $region_line);
 push @regions, $liner[1];

 $num_g_c{$liner[1]}=0;
 $num_bases{$liner[1]}=0;
 $windows{$liner[1]}=1;

}

if($fasta_open &&  $regions_open){
if ($num==1){open(OUTFILE1, ">$output_dir/chr".$chr."_OUT") || print ("Can't open output file: $output_dir/chr$chr.txt\n")}
elsif($num==2) {open(OUTFILE2, ">$output_dir/chr".$chr."_OUT")|| print ("Can't open output file: $output_dir/chr$chr.txt\n")}
elsif($num==3) {open(OUTFILE3, ">$output_dir/chr".$chr."_OUT")|| print ("Can't open output file: $output_dir/chr$chr.txt\n")}
elsif($num==4) {open(OUTFILE4, ">$output_dir/chr".$chr."_OUT")|| print ("Can't open output file: $output_dir/chr$chr.txt\n")}




$window_size=$dist_each_dir*2/$num_windows; 
my $base_coord=$flag=0;
my $line_num=0;
while(my $line=<INPUTFILE>){
   if($line_num>0){
      chomp $line;
my @fasta_liner=split(//, $line);
my $fasta_line_counter=0;
foreach( $fasta_line_counter=0; $fasta_line_counter<=$#fasta_liner;$fasta_line_counter++){
   
   
my $base=$fasta_liner[$fasta_line_counter];
$base_coord++;


if($flag>0){
   
splice @regions, 0,$flag;$flag=0;
}
 

foreach(@regions){
   
my $region=$_;

my $num_window=$windows{$region};
$start=$region-$dist_each_dir+(($num_window-1)*$window_size);
$end=$region-$dist_each_dir+(($num_window)*$window_size);

if(eof and $fasta_line_counter eq $#fasta_liner){ #if end of FASTA file is reached
  
if($base_coord< $end and $base_coord>= $start){
$num_bases{$region}++;
if($base=~/[GgCc]/){
   $num_g_c{$region}++;
}
 elsif($base =~/[Nn]/){
   $num_N{$region}++;
 }
 }

my $gc_percent;
   if($num_N{$region}<$num_bases{$region}){
     $gc_percent=$num_g_c{$region}/$num_bases{$region}; 
   }
   else{
      
      $gc_percent="NA";
   }
calcGCContentInWindow($gc_percent, $region, $start, $end) if defined($num_bases{$region});
undef $num_bases{$region};
undef $num_g_c{$region};
undef $num_N{$region};


$end2=$end;
$start2=$start;

while($windows{$region}<$num_windows){
  
    calcGCContentInWindow("NA", $region, $start2, $end2);
    $windows{$region}++;
   $end2=$region-$dist_each_dir+(($windows{$region})*$window_size);
   $start2=$region-$dist_each_dir+(($windows{$region}-1)*$window_size);
   
}

if($windows{$region}==$num_windows)
{
   if($base_coord>=$end2){
    $final_end=$region-$dist_each_dir+(($windows{$region})*$window_size);
    $final_start=$region-$dist_each_dir+(($windows{$region}-1)*$window_size);
     calcGCContentInWindow("NA", $region, $final_start, $final_end);
           

   }
   $flag++;
undef $num_bases{$region};
undef $num_g_c{$region};
undef $windows{$region};
undef $num_N{$region};
    }
   
}

   elsif($base_coord< $end and $base_coord>= $start){
      
      $num_bases{$region}++;
if($base=~/[GgCc]/){
   $num_g_c{$region}++;
}
 elsif($base =~/[Nn]/){
   $num_N{$region}++;
 }
    
  
    }
    ########## coordinate exceeds allowed range
elsif($base_coord>=$end ){
   
if(defined($num_bases{$region})){
   my $gc_percent;
   if($num_N{$region}<$num_bases{$region}){
     $gc_percent=$num_g_c{$region}/$num_bases{$region}; 
   }
   else{
      
      $gc_percent="NA";
   }
   
 calcGCContentInWindow($gc_percent, $region, $start, $end);}
$num_window=$windows{$region};
undef $num_bases{$region};
undef $num_g_c{$region};
undef $num_N{$region};


$end2=$end;
$start2=$start;

while($base_coord>=$end2 and $windows{$region}<$num_windows){
 
    
     calcGCContentInWindow("NA", $region, $start2, $end2);
    $windows{$region}++;
    $end2=$region-$dist_each_dir+(($windows{$region})*$window_size);
   $start2=$region-$dist_each_dir+(($windows{$region}-1)*$window_size);
    
   
}

if($windows{$region}<=$num_windows and $base_coord<$end2){
   
   $num_bases{$region}++;
      if($base=~/[GgCc]/){
   $num_g_c{$region}++;
}
 elsif($base =~/[Nn]/){
   $num_N{$region}++;
 }
       
        
    }

else
{
   
  
   if($base_coord>=$end2 and $windows{$region}==$num_windows){
      
    $final_end=$region-$dist_each_dir+(($windows{$region})*$window_size);
    $final_start=$region-$dist_each_dir+(($windows{$region}-1)*$window_size);
    
      calcGCContentInWindow("NA", $region, $final_start, $final_end);
           

   }

         $flag++;
undef $num_bases{$region};
undef $num_g_c{$region};
undef $windows{$region};
undef $num_N{$region};

    }
    }


else{ #sgr coordinate is smaller than any of the coordinates in the current window of hte current peak. We can skip all subsequent peaks because that will also be true for those peaks. 
   
   
   last;}


#####end for and while loops
}
}
}
   else{
      $line_num++;
   }
   
  last if(scalar @regions == 0);
}

print "Process $num: Done with chromosome $chr...\n";

  foreach(@regions){
 
   for(my $k=1; $k<=$num_windows; $k++){
      
    
       if($num == 1){print OUTFILE1 "$_\tNA\n";}
    elsif($num == 2){print OUTFILE2 "$_\tNA\n";}
    elsif($num == 3){print OUTFILE3 "$_\tNA\n";}
    elsif($num==4){print OUTFILE4 "$_\tNA\n";}
    
      }
  }


close(INPUTFILE);
close(REGIONSFILE);
    if($num == 1){close OUTFILE1;}
    elsif($num == 2){ close OUTFILE2}
    elsif($num == 3){close OUTFILE3;}
    elsif($num==4){close OUTFILE4;}
  
    
}

undef %num_bases;
undef %num_g_c;
undef %windows;
undef %num_N;

}
}

sub  calcGCContentInWindow {

my $gc_perc=$_[0];
if($gc_perc ne ""){
    if($num == 1){print OUTFILE1 "$_[1]  $gc_perc range: $_[2] and $_[3] \n";}
    elsif($num == 2){print OUTFILE2 "$_[1]  $gc_perc range: $_[2] and $_[3] \n";}
    elsif($num == 3){print OUTFILE3 "$_[1]  $gc_perc range: $_[2] and $_[3] \n";}
    elsif($num==4){print OUTFILE4 "$_[1]  $gc_perc range: $_[2] and $_[3] \n";}
    elsif($num==5){print OUTFILE5 "$_[1]  $gc_perc range: $_[2] and $_[3] \n";}
    
    

}
    
    
    
}
sub handleFileOpenError{
   
   my $type=$_[0];
   my $chr=$_[1];
   
   if($type eq "fasta"){
      $sgr_open=0;
      
   }
   elsif($type eq "region"){
      
      $file_open=0;
   }
   print ("Skipping chromosome $chr\n");
}


