#!/usr/bin/perl -w
open IN,"$ARGV[0]" or die;
open IN_1,"$ARGV[1]" or die;
$length=$ARGV[2]
$out=$ARGV[3]
if(! open $out_fh, '>', $out){die "Can't write '$out':$!";}
while(<IN>){
 if (/\>(Chr[\d]+)/){ 
		$name=$1;
	}else{
		s/[^ATCGNatcgn]//g;
		$seq{$name}.=$_;
	}
}
close IN;

while(<IN_1>){
s/\n//g;
@t=split /\t/;
if(defined $seq{$t[0]}){
$start=int(($t[1]+$t[2])/2)-int($length/2);
$s= substr($seq{$t[0]},$start,$length);
print $out_fh ">$t[0]:$t[1]:$t[2]\n$s\n";
}
}
close IN_1;
