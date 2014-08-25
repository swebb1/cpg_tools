#!/usr/bin/perl -w
use Getopt::Long;
use strict;

#
# CpG Calculator that deals with Very large files 
#
#

my %opt = (
	outfile => "CpGout.tsv", 
	window  =>  500,
	step    =>  50,
	pattern => "",
	nthresh => 0.5,
	all => 0,
	);
            
GetOptions(\%opt, "infile=s", "outfile=s", "window=i", "step=i", "pattern=s", "nthresh=f","rc");

die "ERROR: must use '--infile' to set name of your sequence file"
  if !$opt{infile};

my $rc;
##print search pattern and convert to upper case
if ($opt{pattern} ne ""){ 
	print "pattern = $opt{pattern}";
	$opt{pattern} =~ tr/a-z/A-Z/;
	$rc = reverseComp($opt{pattern});
}


open  (IN, $opt{infile})|| die $!;

open (OUT, ">$opt{outfile}")|| die "cannot write to  $opt{outfile}";
#Print headings
print OUT join("\t", 
	       "#Chr",
	       "Start",
	       "Stop",
	       "CpGo",	       
	       "CpGoe",
	       "CpGdensity",
	       "GC",
	       "GCskew",
	       "Ns",
	       "Pattern",
	       "Length"."\n",
	       );

my $start=1;
my $seq = "";
my $chr;
my ($e,$g_c,$o,$skew,$patt,$CG,$oe);


MAIN: while (<IN>){

    if (/>/){
	$start =0;
	$seq = "";	
	chomp;
	$chr = $_;
	$chr =~ tr/>//d;
	$chr =~ s/^\s+//;
	$chr =~ s/\s.*//;
	
	next MAIN;
    }

    #get seq and convert to UC
    chomp;
    my $s = $_;
    $s =~ tr/a-z/A-Z/;
    $seq .= $s;

    #get sequence length and move on if length is smaller than window
    my $len = length($seq);
    
    next MAIN if($len < $opt{window});

    while($len >= $opt{window} ){
	
	my $winSeq;
	if ($seq =~ /^(.{$opt{window}})/){ 
	    $winSeq = $1;
	} 

	$o = observed_CG($winSeq); #calculate observed CpG on window seq including Ns
	my $lw = length($winSeq);
	my $d = ($o/$lw)*100;

	$patt = 0;
	if ($opt{pattern} ne ""){
		$patt = observed_pattern($winSeq,$rc); #calculate observed pattern on window seq including Ns
	}

	#remove Ns
	$winSeq =~ tr/N//d;
	my $n = length($winSeq);
	if ($n > 0){
		($e, $g_c,$skew) = expected_CG($winSeq); # calculate expected and gc content with Ns removed
	}
	else{   #sequence is all Ns
		$g_c = 0; 
		$e = 0;
	}



	if ($e == 0){ #no Gs or no Cs so expected is zero, cannot divide by 0
		$oe = "NA";
	}
	else{
		$oe = sprintf("%.3f", $o/$e);
	}

	my $nPer = sprintf("%.3f",1-($n/$opt{window}));
	if ( ($nPer) > $opt{nthresh}){ #too many Ns

	    $start += $opt{step};
	    $seq =~ s/^.{$opt{step}}//;
	    $len = length($seq);

	    next;
	}


	#print table
	print OUT join("\t", 
		       $chr,	
		       $start,
		       $start + ($opt{window}),	       
		       $o,
		       $oe,
		       $d,
		       $g_c,
		       $skew,
		       $nPer,	
		       $patt,
		       $lw."\n",
		       );
	$start += $opt{step};
	$seq =~ s/^.{$opt{step}}//;
	$len = length($seq);
    }
}


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
    my $e = 0;
    my $seq = shift;
    my $G = $seq =~ tr/G/G/;
    my $C = $seq =~ tr/C/C/;
    my $len = length($seq);
    $e = ( ($G/$len) * ($C/$len) * ($len-1));

    my $gcskew;
    if($G+$C==0){$gcskew=0.5;}
    elsif($G==0){$gcskew=0}
    else{
       $gcskew = $G/($G+$C);
    }
    $gcskew = sprintf("%.3f",$gcskew);


    return($e, sprintf("%.3f", (($G + $C)/$len)),$gcskew);
}

sub reverseComp{
    my $string = $_[0];
    my $revcom = reverse $string;
    $revcom =~ tr/ACGT/TGCA/;
    return $revcom;
}
