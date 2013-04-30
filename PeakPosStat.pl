# Collect Statistics for Peak Positions
#!/usr/bin/perl
my $VERSION = "1.02"; my $AUTHOR = "Prilusky"; my $YEAR = "2007";
# from Amnon Amir amnon.amir@weizmann.ac.il
$|++;
use File::Spec;
use Bio::Trace::ABIF;

# The input directory to scan
# $indir="E:\\Amnon\\MultiSequencing\\Seqs\\";
my $indir = $ARGV[0];

# The output directory where to save the statistics file
# $outdir=">E:\\Amnon\\MultiSequencing\\SeqsOut\\";
  my $outdir = $ARGV[1] || "."; 

my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $date = sprintf("%04d%02d%02d",$year + 1900,$mon + 1,$mday);
my $outF = File::Spec->catfile($outdir,"${date}_stat.txt");

# Init the ABI library    
my $abif = Bio::Trace::ABIF->new();
   
# The number of nucleotides before
$nbef=5;
  
# The number of nucleotides after (& including current)
$naf=1;

# The total number of nucleotides to consider
$ntot=$nbef+$naf;
  
# initialize the result array lsb=smallest position 
for ($i=0;$i<4**$ntot;$i++) {
	# The relative peak pos sum
	$farray[$i]=0;
	# The relative peak pos square sum
  	$fsqarray[$i]=0;
	# The number of occurences of this sequence
  	$fnum[$i]=0;
  	# The normalized peak pos sum
  	$nfarray[$i]=0;
	# The normalized peak pos square sum
  	$nfsqarray[$i]=0;
}

# expected two level dirs completed_runs/0209/xxx.ab1
  my $numIncludedFiles;
  foreach my $f (glob("$indir/*/*ab1")) {
       $abif->open_abif($f);
       print "Testing Sequence ",$abif->sample_name(), "\n";
       $numIncludedFiles++;

	my @quality_values = $abif->quality_values();
	my $sequence = $abif->sequence();
	my @bl = $abif->base_locations();

	$len=@bl;

	print "sequence length - $len \n";


	# number of subsequences which contain a "N"
	$totnotok=0;
	# The position from which to start scanning (to eliminate initial intensity increase)
	$cnucstart=250;
	# The position from which to end scanning (to eliminate end decay)
	$cnucend=100;
	
	# Get the scaling factor
	$sf = $bl[$len-$cnucend]-$bl[$cnucstart];
	$sf = $sf / ( ($len-$cnucend) - $cnucstart);
	
	# scan all the positions in the sequence
	for ($cnuc=$cnucstart+1;$cnuc<($len-$ntot-$cnucend);$cnuc++) {
		$cval=0;
		$ok=1;
		# scan the bases of the local nucleotide
		for ($cpos=0;$cpos<$ntot;$cpos++) {
			if ($quality_values[$cpos+$cnuc]<0) {
				$ok=0;
			} else {
				$tmpval = (index("ACGT",substr($sequence,$cnuc+$cpos,1)) * (4 ** $cpos));
				# if $tmpval<0 means it is not ACGT
				if ($tmpval<0) {
					$ok=0;
				} else {
					$cval +=$tmpval; 
				}
			}
		}
		# if all subsequence was ok (no N)
		if ($ok==1) {
			$fnum[$cval]++;
			$newval = $bl[$cnuc+$nbef]-$bl[$cnuc+$nbef-1];
			# And add the new normalized intensity to the local sequence intensity
	    	if ($newval>0) {
			  $farray[$cval] += $newval;
			  $fsqarray[$cval] = $fsqarray[$cval] + (($newval ** 2));
			  $nfarray[$cval] += $newval/$sf;
			  $nfsqarray[$cval] = $nfsqarray[$cval] + ((($newval/$sf) ** 2));
	    	}
		} else {
			$totnotok++;
		}
	}

	print "total bad bases - ",$totnotok,"\n";
	$abif->close_abif();
}


# Save the results
open(ofile, ">$outF");
print ofile "# Total $numIncludedFiles files\n";
for ($i=0;$i<4**$ntot;$i++) {
	if ($fnum[$i]>0) {
		print ofile $i," ",$fnum[$i]," ",($farray[$i]/$fnum[$i]);
		print ofile " ",(($fsqarray[$i] - ($farray[$i] ** 2)/$fnum[$i]) / $fnum[$i]);
		print ofile " ",($nfarray[$i]/$fnum[$i]);
		print ofile " ",(($nfsqarray[$i] - ($nfarray[$i] ** 2)/$fnum[$i]) / $fnum[$i]),"\n";
	} else {
		print ofile $i," ",$fnum[$i]," 0 0 0 0\n";
	}
}
close(ofile);

  
  
