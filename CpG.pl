#!/usr/bin/perl -w


use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;


%opt = (
	outfile => "CpG.out", 
	pattern => "",
	);
            
GetOptions(\%opt, "infile=s", "outfile=s", "nthresh=f", "pattern=s","rc");

die "ERROR: must use '--infile' to set name of your sequence file"
  if !$opt{infile};

open (OUT, ">$opt{outfile}") || die "cannot write to  $opt{outfile}";

$in  = Bio::SeqIO->new('-file' => $opt{infile}, 
                        '-format' => 'Fasta'); 

print OUT join("\t", 
	       "#Chr",
	       "CpGo",	       
	       "CpGoe",
	       "CpGdensity",
	       "GC",
	       "GCskew",
	       "Ns",
	       "Pattern",
	       "length"
	       ),"\n";


##print search pattern and convert to upper case
if ($opt{pattern} ne ""){ 
	print "pattern = $opt{pattern}";
	$opt{pattern} =~ tr/a-z/A-Z/;
}


#read in sequence
while (my $seqobj = $in->next_seq()){
    my $seq = $seqobj->seq; #get sequence
    $seq =~ tr/a-z/A-Z/;

    my $n = $seq =~ tr/N/N/; #count Ns

    #check sequence length
    my $seqlen = length($seq);
    unless ($seqlen){
	print  $seqobj->display_id, " has zero length\n" ;
	next;
    }

    #get stats
    my $gc;
    my $skew;
    my $e;

    if ($n < $seqlen){
	$e = expected_CG($seq); # calculate expected and gc content with Ns removed
	$gc = ($seq =~ tr/GC/GC/);
	my $g = ($seq =~ tr/G/G/);
	if($gc==0){$skew=0.5;}
	elsif($g==0){$skew=0}
	else{
	    $skew = $g/$gc;
	}
	$gc = sprintf("%.3f",($gc/ ($seqlen - $n)));
        $skew = sprintf("%.3f",$skew);
    }
    else{   #sequence is all Ns
	$gc = 0; 
	$skew = 0.5;
	$e = 0;
    }



    my $o = observed_CG($seq);

    if ($e==0){
	$oe = "NA";
    }
    else{
    	$oe = sprintf("%.3f",$o/$e);
    }
    my $d = ($o/$seqlen)*100;

    #search for pattern
    my $patt = 0;
    if ($opt{pattern} ne ""){
	my $revPattern = reverseComp($opt{pattern});
	$patt = observed_pattern($seq,$revPattern); #calculate observed pattern on window seq including Ns
    }

    #print table
    my $nPer = sprintf("%.3f",($n/$seqlen));
    if (($nPer) <= $opt{nthresh}){ #not too many Ns
    	print OUT join("\t", 
	       $seqobj->display_id,
	       $o,
	       $oe,
	       $d,
	       $gc,
	       $skew,
	       $nPer,
	       $patt,
	       $seqlen,
	       ),"\n";
    }

}

close (OUT);


sub observed_pattern{
	my $seq = shift;
	my $rc = shift;
	my $pattern = 0;
	$pattern += $seq =~ s/$opt{pattern}/$opt{pattern}/g;
        if(exists $opt{rc}){
	    $pattern += $seq =~ s/$rc/$rc}/g;
	}
	return $pattern;
}

sub observed_CG{
    my $seq = shift;
    my $CG = 0;
    $CG += $seq =~ s/CG/CG/g;
    return $CG;
}

sub expected_CG{
    $e = 0;
    my $seq = shift;
    my $G = $seq =~ tr/G/G/;
    my $C = $seq =~ tr/C/C/;
    my $n = $seq =~ tr/N/N/;
    my $len = length($seq) - $n;
    if ($len>0){
    	$e = ( ($G/$len) * ($C/$len) * ($len-1));
    }
    else{
	$e = 0;
    }	
    return sprintf("%.3f",$e);
}

sub reverseComp{
    my $string = $_[0];
    my $revcom = reverse $string;
    $revcom =~ tr/ACGT/TGCA/;
    return $revcom;
}
