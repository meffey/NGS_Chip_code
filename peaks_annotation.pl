#!/usr/bin/perl -w  
my $in=$ARGV[0];
my $in_1=$ARGV[1];
if (! defined $in){die "Usage:$0 filename";}
 my $out=$in_1;
 my $out_1=$in_1;
 my $out_2=$in_1;
 $out=~s/(\.\w+)?$/.gene/;
 $out_1=~s/(\.\w+)?$/.peaks/;
 $out_2=~s/(\.\w+)?$/.genes/;
if(! open $in_fh, '<', $in){die "Can't open '$in':$!";}
if(! open $in_1_fh, '<', $in_1){die "Can't open '$in_1':$!";}
if(! open $out_fh, '>', $out){die "Can't write '$out':$!";}
if(! open $out_1_fh, '>', $out_1){die "Can't write '$out_1':$!";}
if(! open $out_2_fh, '>', $out_2){die "Can't write '$out_1':$!";}
@gff=<$in_fh>;
close $in_fh;

@summit=<$in_1_fh>;

foreach$s(@summit){
  chomp$s;
  $peaks++;
  @r=split/\t/,$s;
  $h=int(($r[1]+$r[2])/2);
  $top{$r[0]}{$h}=$peaks;
 }
 
 foreach$p(@gff){
  chomp$p;
  @t=split/\t/,$p;
  $a=0;
if(($p!~/transpo/)&&($p=~/Chr/)){
   if(($t[2]=~/gene/)&&($t[6]eq"+")){
   $kb=$t[3]-1000;
   $t[8]=~/ID=(.*?);/;
   for($a=$kb;$a<$t[3];$a++){
     if(defined$top{$t[0]}{$a}){
	 $ann{$1}{$top{$t[0]}{$a}}.="promoter";
	 $two{$1}.="promoter";
	 $one{$top{$t[0]}{$a}}="promoter";
	 }
   }
   for($a=$t[3];$a<$t[4];$a++){
        if(defined$top{$t[0]}{$a}){
	 $ann{$1}{$top{$t[0]}{$a}}.="gene";
	 $two{$1}.="gene";
	 $one{$top{$t[0]}{$a}}="gene";
   }
   }
  }elsif(($t[2]=~/gene/)&&($t[6]eq"-")){
   $kb=$t[4]+1000;
   $t[8]=~/ID=(.*?);/;
      for($a=$t[4];$a<$kb;$a++){
     if(defined$top{$t[0]}{$a}){
	 $ann{$1}{$top{$t[0]}{$a}}.="promoter";
	 $two{$1}.="promoter";
	 $one{$top{$t[0]}{$a}}="promoter";
	 }
   }
   for($a=$t[3];$a<$t[4];$a++){
        if(defined$top{$t[0]}{$a}){
	 $ann{$1}{$top{$t[0]}{$a}}.="gene";
	 $two{$1}.="gene";
	 $one{$top{$t[0]}{$a}}="gene";
   }
   }
  }elsif(($t[2]eq"exon")||($t[2]eq"five_prime_UTR")||($t[2]eq"three_prime_UTR")){
  $t[8]=~/Parent=(.*?)\.[0-9]/; 
    $m=$1;
   if((defined $two{$m})&&($two{$m}=~/gene/)){
    for($a=$t[3];$a<$t[4];$a++){
        if(defined$top{$t[0]}{$a}){
      $ann{$m}{$top{$t[0]}{$a}}.="$t[2]";
	  $one{$top{$t[0]}{$a}}="$t[2]";
	   }
	  }
     }
	}	
   }
  }

foreach$key(sort keys %ann){
   foreach $ppp(keys %{$ann{$key}}){
   print $out_fh "$key\t$ppp\t$ann{$key}{$ppp}\n";
   }
 }

foreach$key(sort keys %one){
 print $out_1_fh "$key\t$one{$key}\n";
 }
 
 foreach$key(sort keys %two){
 print $out_2_fh "$key\t$two{$key}\n";
 }
 
close $in_1_fh;
close $out_fh;
close $out_1_fh;
close $out_2_fh;
