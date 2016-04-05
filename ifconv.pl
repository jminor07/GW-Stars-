#!/usr/bin/perl -w

#IFCONV: determines if the model is converging
#LS 2005 v1.0

#NOTES:   it works with a 72 layer model, no matter how much iterations you did since
#         it uses the very last error report. It assumes you have an extra
#         line after the 72th layer (e.g. the output from a DATE shell command),
#         so yu do not have to worry about it; should you have more than one,
#         edit the line:
#         `tail -n 73 $ARGV[0] > ifconv.prov`;
#         putting the right number instead of 73.
#         Should you like to vary the convergence criteria, simply change the 
#         two variables $max_flux and $max_deri.
#         You have to pass the filename as an argument: ifconv model.dat

#setting the convergence criteria
$max_flux=1.0;      #max flux error 1%
$max_deri=10.;      #max flux derivative error 10%

#Some subroutines...

sub vecminmax{
#Subroutine vecminmax: looks for ABSOLUTE max and min of a vector and puts them in
#$v_min  and $v_max variables. Also returns $n_min and $n_max, positions in the vector
my @vec=@_;
$v_max=abs($vec[0]);$v_min=$v_max;
$n_min=0;$n_max=0;my $j=0;
foreach $vec (@vec){
	$vec=abs($vec);
	if($v_min>$vec){$v_min=$vec;
	$n_min=$j;}
	if($v_max<$vec){$v_max=$vec;
	$n_max=$j;}
	$j++;
}
return $v_min,$v_max,$n_min,$n_max;
}

sub whichlay{
#This program speciifc subroutine searches and prints out the 
#layers not respecting the convergence criteria
print "LAYER TEMP    K_ROSS          ERR      DERIV \n \n";
foreach $lay (@lay){
	$lay--;
	if ((abs($flux[$lay]) >= $max_flux) or (abs($deri[$lay]) >= $max_deri)){
	my $layp=$lay+1;
		print "$layp $teff[$lay] $ross[$lay] $flux[$lay] $deri[$lay] \n";
	}
}	
}

#now we read in the filename from the command line
if(-e $ARGV[0]){
	print qq(Checking convergence from file $ARGV[0]\n);}
else
	{die qq(You gave me the wrong filename);}

#Now we look at what we have in the last line
#To be done


#Since we have one extra line...
`rm -f ifconv.prov ifconv.in`;
`tail -n 72 $ARGV[0] > ifconv.prov`;
`head -n 72 ifconv.prov > ifconv.in`;

#readin. We read the significant valuess
open(INFI, "ifconv.in");
@lay=1..72;@teff=@ross=@flux=@deri=@lay;
$i=0;
while (<INFI>) {
	chomp;
	$lay[$i]=substr $_,0,4; #Layer number
	$teff[$i]=substr $_,15,9; #Temperature in the layer
	$ross[$i]=substr $_,79,11; #Rosseland optical depth
	$flux[$i]=substr $_,112,12; #Flux error
	$deri[$i]=substr $_,124,9; #Flux derivative error.
	#print $lay[$i],$teff[$i],$ross[$i],$flux[$i],$deri[$i],"\n";
	$i++;
}
close (INFI); `rm -f ifconv.in ifconv.prov`;
#print "flux 0 is $flux[0] \n";

#Now we check the convergence and print out.
vecminmax @flux;$fmax=$v_max;$nfmax=$n_max;
$rsf=$ross[$nfmax];$tef=$teff[$nfmax];
vecminmax @deri;$dmax=$v_max;$ndmax=$n_max;
$rsd=$ross[$ndmax];$ted=$teff[$ndmax];
$nfmax++;$ndmax++;
print "\n wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww \n \n";
if(($fmax<$max_flux) and ($dmax<$max_deri)){
	print "The model is CONVERGING \n \n";
	print "Max flux error $fmax at layer $nfmax \n";
	print "Where T = $tef and Rosseland Depth is $rsf \n";
	print "Max flux derivative error $dmax at layer $ndmax \n";
	print "Where T = $ted and Rosseland Depth is $rsd \n \n";
}else{
	print "The model is NOT CONVERGING \n \n";
	print "Max flux error $fmax at layer $nfmax \n";
	print "Where T = $tef and Rosseland Depth is $rsf \n";
	print "Max flux derivative error $dmax at layer $ndmax \n";
	print "Where T = $ted and Rosseland Depth is $rsd \n \n";
	print "Here follow the nonconverging layers... \n \n";
	whichlay;
}
#print $v_min,$v_max,$n_min,$n_max,"\n";






