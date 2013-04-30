# Collect Statistics for Peak Heights
use Bio::Trace::ABIF;

# The input directory to scan
$indir="E:\\Amnon\\MultiSequencing\\Seqs\\";

# The output directory where to save the statistics file
$outdir=">E:\\Amnon\\MultiSequencing\\SeqsOut\\";


# Init the ABI library    
my $abif = Bio::Trace::ABIF->new();
   
# The number of nucleotides before
$nbef=5;
  
# The number of nucleotides after (& including current)
$naf=1;

# The total number of nucleotides to consider
$ntot=$nbef+$naf;
  
# The window size for normalizing the peak hight
$avgwinsize=10;

   
# initialize the result array lsb=smallest position 
for ($i=0;$i<4**$ntot;$i++) {
	# The relative peak height sum
	$farray[$i]=0;
	# The relative peak height square sum
  	$fsqarray[$i]=0;
	# The number of occurences of this sequence
  	$fnum[$i]=0;
}


# Open the input directory and go over all files
opendir(DIR,$indir) 
   || die "NO SUCH Directory: $indir";

while ($filename = readdir(DIR) )
{
	$abif->open_abif($indir.$filename);
	print "Testing Sequence ",$abif->sample_name(), "\n";

	my @quality_values = $abif->quality_values();
	my $sequence = $abif->sequence();
	my @bl = $abif->base_locations();

	$len=@bl;

  	# get the traces (my annotation a=1,c=2,g=3,t=4)
	@trace1 = $abif->analyzed_data_for_channel(2);
	@trace2 = $abif->analyzed_data_for_channel(4);
	@trace3 = $abif->analyzed_data_for_channel(1);
	@trace4 = $abif->analyzed_data_for_channel(3);

	# calculate the total amplitude trace (for relative amplitude)
	for ($i=0;$i<@trace1;$i++) {
  		$ttrace[$i] = $trace1[$i]+$trace2[$i]+$trace3[$i]+$trace4[$i];
	}  
	print "sequence length - $len \n";


	# number of subsequences which contain a "N"
	$totnotok=0;
	# The position from which to start scanning (to eliminate initial intensity increase)
	$cnucstart=250;
	# The position from which to end scanning (to eliminate end decay)
	$cnucend=100;
	# scan all the positions in the sequence
	for ($cnuc=$cnucstart+$avgwinsize;$cnuc<($len-$ntot-$cnucend)-$avgwinsize;$cnuc++) {
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
			# Calculate the local intensity avg
			$avginten=0;
			$avgpoints=0;
			for ($cp=($cnuc+$nbef-$avgwinsize);$cp<=($cnuc+$nbef+$avgwinsize);$cp++) {
				if ($cp!=$cnuc+$nbef) {
					$avginten += $ttrace[$bl[$cp]+1];
					$avgpoints++;
				}
			}
			$avginten = $avginten / $avgpoints;
			$fnum[$cval]++;
    		$newval=0;
			if (index("ACGT",substr($sequence,$cnuc+$nbef,1))==0) {
				$newval = $trace1[$bl[$cnuc+$nbef]+1] / $avginten;
			}
			if (index("ACGT",substr($sequence,$cnuc+$nbef,1))==1) {
				$newval = $trace2[$bl[$cnuc+$nbef]+1] / $avginten;
			}
			if (index("ACGT",substr($sequence,$cnuc+$nbef,1))==2) {
				$newval = $trace3[$bl[$cnuc+$nbef]+1] / $avginten;
			}
			if (index("ACGT",substr($sequence,$cnuc+$nbef,1))==3) {
				$newval = $trace4[$bl[$cnuc+$nbef]+1] / $avginten;
			}
			# And add the new normalized intensity to the local sequence intensity
	    	if ($newval>0) {
			  $farray[$cval] += $newval;
			  $fsqarray[$cval] = $fsqarray[$cval] + (($newval ** 2));
	    	}
		} else {
			$totnotok++;
		}
	}

	print "total bad bases - ",$totnotok,"\n";
	$abif->close_abif();
}
closedir(DIR);


# Save the results
open(ofile, $outdir.'stat.txt');
for ($i=0;$i<4**$ntot;$i++) {
	if ($fnum[$i]>0) {
		print ofile $i," ",$fnum[$i]," ",($farray[$i]/$fnum[$i]);
		print ofile " ",(($fsqarray[$i] - ($farray[$i] ** 2)/$fnum[$i]) / $fnum[$i]),"\n";
	} else {
		print ofile $i," ",$fnum[$i]," 0 0\n";
	}
}
close(ofile);

  
  
