#!/usr/bin/perl -w   
open IN,"$ARGV[0]" or die;
open IN_1,"$ARGV[1]" or die;
$out=$ARGV[2];
if(! open $out_fh, '>', $out){die "Can't write '$out':$!";}
 
 while(<IN>){
 s/\n//g; 
 @a=split /\t/;
   for($s=$a[1];$s<=$a[2];$s++){
  $one{$a[0]}{$s}=1;
    }
  }
  
   while(<IN_1>){
 s/\n//g; 
 @a=split /\t/;
   for($s=$a[1];$s<=$a[2];$s++){
    if(defined $one{$a[0]}{$s}){
     print $out_fh "$_\n";
	   last;
      }
    }
}

close IN;
close IN_1;