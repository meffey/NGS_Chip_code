#!/usr/bin/perl -w 
my $in=$ARGV[0];
my $in_1=$ARGV[1];
if (! defined $in){die "Usage:$0 filename";}
 $out=$ARGV[2];
if(! open $in_fh, '<', $in){die "Can't open '$in':$!";}
if(! open $in_1_fh, '<', $in_1){die "Can't open '$in_1':$!";}
if(! open $out_fh, '>', $out){die "Can't write '$out':$!";}
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
   $kb=$t[3]+400;
   $t[8]=~/ID=(.*?);/;
   for($a=$t[3];$a<$kb;$a++){
     if(defined$top{$t[0]}{$a}){
	 $ann{$1}{$top{$t[0]}{$a}}.="tss";
	 }
   }
  }elsif(($t[2]=~/gene/)&&($t[6]eq"-")){
   $kb=$t[4]-400;
   $t[8]=~/ID=(.*?);/;
      for($a=$kb;$a<$t[4];$a++){
     if(defined$top{$t[0]}{$a}){
	 $ann{$1}{$top{$t[0]}{$a}}.="tss";
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

 
close $in_1_fh;
close $out_fh;

