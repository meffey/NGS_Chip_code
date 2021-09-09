#!/usr/bin/perl -w  
open IN,"$ARGV[0]" or die;
$out=$ARGV[1];
if(! open $out_fh, '>', $out){die "Can't write '$out':$!";}
while(<IN>){
  s/\n//g;
  @t=split /\t/; 
  @a=split /:/,$t[0];
  @b=split /#/,$t[1];
   @n=split /:/,$b[0];
   @n1=split /:/,$b[1];
   @n2=split /:/,$b[2];
     if($a[3] eq "-"){
    $x=$n[2]+1;$x1=$n1[1]-1;$x2=$n1[2]+1;$x3=$n2[1]-1;
	    	$first="$a[0]:$x:$x1:$a[3]";
	$secend="$a[0]:$x2:$x3:$a[3]";
      }else{
	$x=$n[1]-1;$x1=$n1[2]+1;$x2=$n1[1]-1;$x3=$n2[2]+1; 	
	$first="$a[0]:$x1:$x:$a[3]";
	$secend="$a[0]:$x3:$x2:$a[3]";  
	  }

	print $out_fh "$t[0]\t$b[0]#$first#$b[1]#$secend#$b[2]\n";
  }
close IN;  
