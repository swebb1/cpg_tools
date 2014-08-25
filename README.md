cpg_tools
=========

CpG.pl
===
Calculate CpG statistics over entire fasta record.

--infile: Input fasta file
--outfile: Output file
--nthresh: Include fasta records with this fraction of Ns or less (use a fraction of 1).


CpGwindow.pl
===
Split fasta record in to a window of defined size and calculate CpG statistics per window.

--infile: Input fasta file
--outfile: Output file
--window: Window length
--step: Step size
--nthresh: Include fasta records with this fraction of Ns or less (use a fraction of 1)

CpGseq.pl
===
Generate random sequences using parameters of length, GC content and CpG observed/expected ratio.

--out: Output file
--length: Approxiamate length of output sequence
--GC: GC% content of output sequence
--CpG: CpG observed/expected ratio of output sequence
