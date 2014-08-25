#!/usr/bin/perl -w
use strict;
use Getopt::Long;

#
# Set Defaults
#

my %opt = (
    length => 1000, 
    GC => 70, 
    CpG => 0.8,
    out => "randomseq.fa",	
    );

my %cutsites;


GetOptions(\%opt, "out=s", "length=i", "GC=i", "CpG=f", "cutsites=s");

die "ERROR: must use '--out' to set name of your sequence files"
    unless ($opt{out}) ;


if($opt{CpG} > 1){
    $opt{CpG} =1;
    print STDERR "Max CpG value is 1,  set to 1"
}

#Set the ranges of allowed CpG values
$opt{dmaxCpG} = $opt{CpG};
$opt{dminCpG} = $opt{CpG} - 0.05;
$opt{dminCpG} = 0 if $opt{dminCpG} < 0;


if ($opt{cutsites}){ # use file if present
    $opt{cutsites} = "/storage/home/galaxy/galaxy_data/tool-data/cutsites" if $opt{cutsites} eq "default"; #set to default if 'default' supplied as input
    open (CUT, $opt{cutsites}) or die "cannot open cutsite "; 
    
    while(<CUT>){ # in fatab format
	chomp;
	tr/a-z/A-Z/;
	my ($site, $cseq) = split;
	if(length($cseq) < 4){
	    print STDERR "ignoring $site as $cseq < 4 characters";
	    next;
	}
	$cutsites{$site} = $cseq; #allow for future changes to add what was found
	my $rsite = "r" . $site; # allow forward and reverse sequences
    my $rseq =  $cseq;
	$rseq =~ tr/ATGC\]\[/TACG\[\]/; # should word on basic regexps
	$rseq = reverse($rseq);
	$cutsites{$rsite} = $rseq;

    }
}



my $out = open(OUT, ">" . $opt{out}) or die "cannot write to $opt{out}";

my $new_seq; 
my $ok = 0;

$new_seq = &rand_seq(&create_seq); # create a seq and randomise it

 SEQ: until( $ok > 0){
     my $eCpG = &expected_CG($new_seq); # observed CpG 
     my $oCpG = &observed_CG($new_seq); # expected  CpG 
#Calculate required observed 
#
# a bit too much kruft, I could simplify
     my $remove;
     if($opt{dminCpG} > 0){
	 my $eoCpG = int($eCpG * ($opt{dminCpG} + 
				  ( ($opt{dmaxCpG} -$opt{dminCpG})/2 )));
	 $remove = $oCpG-$eoCpG;
     }else{
	 $remove = $oCpG;
     }

     die "Sorry: cannot add CpGs yet" if($remove < 0);
     $new_seq = &edit_dna($new_seq, $remove);  # find  and change all CpG's 
     $eCpG = &expected_CG($new_seq); #cal CpG stats again to check if all ok
     $oCpG = &observed_CG($new_seq);

     redo SEQ unless # if not ok then try again
	 (
	  ($oCpG/$eCpG <= $opt{dmaxCpG})
	  and ($oCpG/$eCpG >= $opt{dminCpG})
	 );
     
     if($opt{cutsites}){
	 my $changed = 1; # nos same for cutsites
	 while($changed == 1){ 
	 # keep looking as randomly changing one might create another 
	     ($new_seq, $changed) = &remove_cutsites($new_seq);
	 }
     }
     $ok=1;
} 



my $title = ">RAND\n"; # print
print OUT "$title", &tidy_seq($new_seq), "\n";    




sub remove_cutsites{
    #wrappe for the mutate sub, allows initiation so then this sub cn be looped
    my $seq = shift;
    my $diff =0;
 
  MUT: foreach my $site (values %cutsites){ # check for sites
      if (my @sites= ($seq =~ /($site)/g)){ #redo if found
	  @sites = &uniq(@sites);
	  foreach $site (@sites){
	      $seq = &mutate($seq, $site); # mutatetakes the seq and the string to change
	      $diff =1; # if found then record that it was changed
	  }
      }	 
  }
    return($seq, $diff);
}


sub mutate{
    my $seq = shift; # seq to change
    my $site = shift; # sub-sequence to change 
    my @locs; 
    my $p = 0;
    
    $p = index($seq, $site);
    if($p > -1){
	push(@locs, $p);
    }else{
	die "no site found for $site in\n $seq" ;
    }	
    while ($p = index($seq, $site, $p+1)){ # index all locations
	last if $p < 0;
	push(@locs, $p);
    }

    my @seq = split //, $seq; # convert scalar to list, 1 base per index

    #change 2 bases at random;    
   while(@locs){
	my $change = shift @locs;
	my $len = length($site) + $change;
        my @pos = ($change..$len);
	my $mutate1 = $pos[int(rand($#pos))];
	my $mutate2 = $pos[int(rand($#pos))];
	while($mutate1 == $mutate2){ # make sure the 2 positions are different
	    $mutate2 = $pos[int(rand($#pos))];
	}
#	print "$mutate1 $mutate2\n";
	foreach my $i ($mutate1, $mutate2){		   
	    $seq[$i] =~ tr/GCTA/CGAT/; # convert keeping same G+C content
	}
    }
    return(join("", @seq));
}


sub create_seq{
    my $T;
    my $A;
    my $C;
    my $G;
    my $nseq = "";
    my $gc = $opt{length} * ($opt{GC}/100);
    my $at = $opt{length} *  ((100- $opt{GC}) /100);
    
    foreach (0..$gc ){
	if(rand(1)> 0.5){
	    $nseq .= "G" ;
	}else{
	    $nseq .= "C" ;
	}
    }


    foreach (0..$at){
	if(rand(1)> 0.5){
	    $nseq .= "T" ;
	}else{
	    $nseq .= "A" ;
	}
    }
    
    $nseq;
}

sub  rand_seq{
    my $seq = shift;
    my @seq = split //, $seq;
    my $new;
    while(@seq){
	my $nuc
	    = splice(@seq, int (rand($#seq +0.5)), 1);
	$new .= $nuc;
    }
    return($new);
}


sub tidy_seq{                   #adds a new line character every 60 bases  
    my $seq = shift;
    my ($new_seq) = "";
    my $start;
    my $end;
    my $seq_len = length($seq);
    my $to_add = int($seq_len/60);

    $end = $start= 0;

    foreach (1..$to_add){
        $new_seq .= substr($seq,$start,60);
        $new_seq .= "\n";
        $start += 60;
    }
    $new_seq .= substr($seq,$start);
    return ($new_seq);
}

    

sub observed_CG{
    my $seq = shift;
    my $CG = 0;
    $CG += $seq =~ s/CG/CG/g;
    $CG;
}


sub expected_CG{
    my $CG = 0;
    my $seq = shift;
    die "zero length sequ deteted" unless length($seq);
    my $G = $seq =~ tr/G/G/;
    my $C = $seq =~ tr/C/C/;
    my $n = $seq =~ tr/N/N/;
    my $len = length($seq) - $n;
    $CG = ( ($G/$len) * ($C/$len) * ($len-1));
    $CG;
}

sub edit_dna{
    my @replace;
    my $seq = shift;
    my $remove = shift;
    my @CGpos;
    my $p = 0;

    while ($p = index($seq, "CG", $p+1)){
	last if $p < 0;
	push(@CGpos, $p);
    }
    my @seq = split //, $seq;
    foreach(1..$remove){
	my $change = splice(@CGpos, int (rand($#CGpos +0.5)), 1);
	$replace[1] = "T";
	$replace[2] = "A";
	if (rand(1) > 0.5){
	    $change++;
	    die $seq[$change] unless( $seq[$change] eq "G");
	    $replace[0] = "C";
#	    $seq[$change] = $replace[int(rand(3))];
	    $seq[$change] = $replace[0];
	}else{
	    die "not C $seq[$change]" unless( $seq[$change] eq "C");
	    $replace[0] = "G";
#	    $seq[$change] = $replace[int(rand(3))];
	    $seq[$change] = $replace[0];
	}
    }
    return(join("", @seq));
}


    
sub uniq{
    my %hash;
    $hash{$_} = 1 foreach (@_);
    return (keys %hash);
}
