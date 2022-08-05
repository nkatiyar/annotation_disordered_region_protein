#!/usr/bin/perl
use strict;
# Neerja Katiyar and Phil Goetz
# For each gene in a genome,
#   For each BLASTP hit,
#     Retrieve BLASTP homology string
#     Compute correlation, over all residues, between BLASTP homology score, & max(anchor, 1-IUPred)
#

use English;
use File::Basename;
use File::Which;
use FindBin qw($Bin);     # Find where this binary is
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use List::Util qw(max min);
use Math::Complex;
use lib "$Bin/lib";
#use lib "$ENV{AA}/src/lib";
use pgoetzUtils;
use Anno::autoAnnotate;
use Anno::Evidence; 

select STDOUT; $OUTPUT_AUTOFLUSH = 1;

# These constants start from 0
use constant HT_HMM_ACCESSION  => 0;
use constant HT_HMM_LEN        => 2;
use constant HT_MATCH_START    => 6;
use constant HT_MATCH_END      => 7;
use constant HT_PROTEIN_START  => 8;
use constant HT_PROTEIN_END    => 9;
use constant HT_DOMAIN_SCORE   => 11;
use constant HT_TOTAL_SCORE    => 12;
use constant HT_DOMAIN         => 13;
use constant HT_DOMAIN_COUNT   => 14;
use constant HT_NAME           => 15;
use constant HT_TRUSTED_CUTOFF => 17;
use constant HT_NOISE_CUTOFF   => 18;

## BER BTAB file columns
use constant BT_ORF              => 0;
use constant BT_QUERY_LEN        => 2;  # Length of query protein + 600 in nt
use constant BT_METHOD           => 3;
use constant BT_ACCESSION        => 5;
use constant BT_QUERY_START      => 6;  # Start on (query+300 each end) of alignment, in nt
use constant BT_QUERY_STOP       => 7;  # Stop on (query+300 each end) of alignment, in nt
use constant BT_HIT_START        => 8;  # Start on hit of alignment, in amino acids
use constant BT_HIT_STOP         => 9;  # Stop on hit of alignment, in amino acids
use constant BT_PCNT_IDENT       => 10;
use constant BT_PCNT_SIM         => 11; # similarity
use constant BT_PRAZE_SCORE      => 12; # or blast score
use constant BT_COMMENT          => 15; # name {species;}  (exp=-.; wgp=.; cg=; genomes=; )
use constant BT_HIT_LEN          => 18; # length of the target gene.  Really.
use constant BT_EVALUE           => 19; # this column is actually blank
use constant BT_PVALUE           => 20; # this column is probably actually EVALUE

# Locate anchor and iupred and their data
my $DisoDir = "/usr/local/devel/prok/function/disorder";
my $motifs_file = "/local/archive_home/pgoetz/disorder/ANCHOR/motifs_elm.txt"; #Motifs file for ANCHOR
$ENV{IUPred_PATH} = "$DisoDir/iupred";

my $DEBUG = 0;       # 1: print info, 2: process fewer items, 3: skip loading dictionary
my $Datadir = '';
my $GapInLen = 0;    # 1 => don't penalize gaps in %ID; penalize then in length
my @RightCount;      # $RightCount[0] = # of times method 0 top blast hit was right
my $RightTries = 0;  # number of tries for @RightCount
my %ResidueCounts;   # $ResidueCounts{A} = # of 'A' in all query proteins
my %ResidueCountsD;  # $ResidueCountsD{A} = # of disordered 'A' in all query proteins
my $Write = 0;       # write to residues files

GetOptions(
	'debug|d!'            => \$DEBUG,
	'data:s'              => \$Datadir,
	'w!'                  => \$Write,
);

# Read in previous values from residue counts
my $ResidueFile = "residues";
&tsvFileToHash("$ResidueFile", 0, 1) if (-e $ResidueFile);
&tsvFileToHash("${ResidueFile}D", 0, 1) if (-e "${ResidueFile}D");

my $blastfile_path = "$Datadir/ber";  #path to BLAST output files
my $iupred_path = "$Datadir/iupred";  #path to IUPRED output files
my $anchor_path = "$Datadir/anchor";  #path to ANCHOR output files
my $cor_path = "$Datadir/disorder_cor";
my $htab_path = "$Datadir/hmm";

my $HmmInfoRef;
print "Loading HMM data from egad\n";
my $egadH = &connectToDB("Sybase", 'SYBPROD', 'access', 'access', '', 'egad');
$HmmInfoRef = &loadHMMinfo($egadH); # hash stores info from egad..hmm2 and egad..hmm3
$egadH->disconnect();

# Dictionary from GO; words counted and compiled by Kevin Galinsky; 7579 entries
# File format: Each line is word\tcount
my $Datadir = "/usr/local/archive/home/pgoetz/anno/data";
my %Dictionary;
my %StrippedToStripped;   # $StrippedToStripped{$name} = $synonym
if ($DEBUG < 3) {
	my $DictFile = "$Datadir/words";
	print "Reading dictionary from $DictFile\n";
	%Dictionary = &tsvFileToHash($DictFile, 0, 1);

	my $SynonymsH = &connectToDB('SQLite', '', '', '', "$Datadir/db", 'thesaurus');
	print "Reading $Datadir/db/thesaurus\n";
	%StrippedToStripped = %{&readSynonyms($SynonymsH)};
	my $numSynos = scalar(keys(%StrippedToStripped));
	print "$numSynos entries in thesaurus\n";
	$SynonymsH->disconnect();
}

# Get SP names from uniprot, not from Panda
my $UniprotH = &connectToDB('SQLite', '', '', '', "$Datadir/db", 'uniprot');

# Get hash of known gene_sym; there are 40,924 of them as of 090626
print "Reading gene symbols from omnium\n";
my $OmniumH = &connectToDB("Sybase", 'SYBPROD', 'access', 'access', '', 'omnium');
my $GSHashRef:shared = &getOmniGenesymHashRef($OmniumH);
$OmniumH->disconnect();

# $PNUHashRef is for renaming
print "Reading PNU\n";
my $PNUfile = "$Datadir/pnu.sto";         # contains Storable hash from PNU
my $PNUHashRef = &getPNUHashRef(1, $PNUfile);    # ref to hash of PNU fullnames

#Reading BLAST files
print "Loading BLASTP output from $blastfile_path\n";
opendir(INFILE, $blastfile_path) or die "Can't open BLASTFILE PATH $blastfile_path\n";
my @files = grep(/\.blastp$/,readdir(INFILE));
closedir(INFILE);

print "Read " .scalar(@files) . " feature names from $blastfile_path\n\n";

my @orfs;
foreach my $blastfile(@files)
{
	my $orf_name = $blastfile;
	if($orf_name =~/^ORF.*/)
	{
		my @orf_nam = split('\.',$orf_name);
		push(@orfs,$orf_nam[0]);
	}	   
}

my $totalCorSum = 0;
my $totalResidues = 0;   # total residues in alignments in blast hits
my $totalCorSumC = 0;
my $totalResiduesC = 0;  # total residues in alignments in blast hits judged correct
my $queryResidues = 0 ;  # total residues in queries
my $queryResiduesD = 0;  # total disorded residues in queries
my $numWiEquiv = 0;      # number of ORFS with equivalogs
my $numDisorderedGenes = 0;
my $numOrfsTried;            # number of orfs not skipped
my @comparisons;         # @comparisons[3] = \{@result = &compareMethods on orf 3}
my @CompareStats;        # @CompareStats[0]->[1] = number of times method 0 beat method 1

if ($DEBUG > 1) { @orfs = @orfs[0 .. 50] }
foreach my $orf (@orfs)
{
	my $orf_cor = "$cor_path/$orf.cor";
	my $iupredPrefix = "$iupred_path/$orf.iupred";
	my $anchor_file = "$anchor_path/$orf.anchor";
	my $btab_file = "$blastfile_path/$orf.nr.btab";
	my $htabFile = "$htab_path/$orf.hmmpfam.htab";

	# Find equivalogs that hit $orf
	# Do this first, since we skip orfs with no equivalogs
	my @equivalogs = &parseHtab($orf, $htabFile, $HmmInfoRef);
	my @equivalogNames;
	if (@equivalogs) {
		@equivalogNames = &getEquivNames(\%Dictionary, $GSHashRef, $PNUHashRef, \%StrippedToStripped, @equivalogs);
		$numWiEquiv++;      # count of ORFs with equivalogs
	}

	#Read and parse IUPRED output file

	my @iupPosn;        # @{$iupPosn[1]} stores iupred result residue positions (long file)
	my @iupProb;        # Stores the disorder probability of iupred residues
	my @iupResidue;     # Just a list
	# Fill in @iupPosn, @iupProb, @iupResidue
	my @iupSuffix = ('short', 'long');  # type of iupred output file to use
	for (my $i=0; $i<@iupSuffix; $i++) {
		my $suffix = $iupSuffix[$i];
		my $iupredFile = "$iupredPrefix.$suffix";
		print "Parsing $iupredFile\n" if $DEBUG;
		my ($posRef, $resRef, $probRef) = &parseIupred($iupredFile);
		$iupPosn[$i] = $posRef;
		$iupProb[$i] = $probRef;
		@iupResidue = @$resRef;
	}

	my $numIuLines = scalar(@{$iupProb[1]});
	die "ERR: Unequal number of entries in the different $iupredPrefix files"
			if scalar(@{$iupProb[0]}) != $numIuLines;

	# Read and parse ANCHOR output file
	print "Parsing ANCHOR output file\n" if $DEBUG;
	my ($posRef, $resRef, $probRef) = &parseAnchor($anchor_file);
	my @anchor_pos = @$posRef;
	my @anchor_res = @$resRef;
	my @anchor_probt = @$probRef;

	# Verify that residues match from iupred and anchor
	for(my $k=0; $k<scalar(@anchor_res); $k++) {
		my $pos = $k+1;   # first position is 1
		my $ares = $anchor_res[$k];
		my $ires = $iupResidue[$k];
		die "Conflicting residue in $anchor_file at pos=$pos: " .
		    "IUPRED[$k]=$ires, ANCHOR[$k]=$ares"
				if $ires ne $ares;
	}

	# Check if number of lines used is same in IUPRED and ANCHOR
	my $protLen = scalar(@anchor_probt);
	die "ERR: $anchor_file has $protLen lines, $iupredPrefix files have $numIuLines lines"
			if $protLen != $numIuLines;

	# Combine iupred.short, iupred.long, and anchor into a disorder score
	# Test whether there is any disorder in this gene
	# Also count whether this is a "disordered" gene
	my $disorderSum = 0;
	my $disoResidueCount = 0;   # number of predicted disordered residues in the query
	my $isDisordered = 0;  # 1 => this gene has 11 disordered residues in a row
	my $runLength = 0;     # number of residues in a row that are disordered
	my @orderProb;
	for (my $j=0; $j<$protLen; $j++) {
		my $posn1 = $iupPosn[0]->[$j];
		my $posn2 = $iupPosn[1]->[$j];
		my $prob1 = $iupProb[0]->[$j];
		my $prob2 = $iupProb[1]->[$j];
		my $posnA = $anchor_pos[$j];
		if (($j+1 == $posn1) && ($posn1 == $posn2) && ($posn1 == $posnA)) {
			my $iup_prob_grt = &max($prob1,$prob2);
			my $orderScore = &getOrder($iup_prob_grt, $anchor_probt[$j]);
			$orderProb[$j] = $orderScore;

			$disorderSum += (1 - $orderScore);
			$ResidueCounts{$iupResidue[$j]}++;
			if ($orderScore < .5) {
				$disoResidueCount++;
				$runLength++;
				$isDisordered++ if ($runLength > 10);
				$ResidueCountsD{$iupResidue[$j]}++;
			}
			else {
				$runLength = 0;    # start over
			}
		}
		else {
			die "ERR: index=$j pos[0]=$posn1 pos[1]=$posn2 pos[ANCHOR]=$posnA";
		}
	}
	$numDisorderedGenes++ if ($isDisordered);
	my $disoFrac = $disoResidueCount/$protLen;
	$queryResiduesD += $disoResidueCount;
	$queryResidues += $protLen;
	print "Disorder in $orf is $disorderSum; disordered residues = $disoResidueCount/$protLen = $disoFrac\n";
	next if ($disoFrac < .05);   # next $orf
	next if !@equivalogs;

	print "Parsing BTAB file $btab_file\n" if $DEBUG;
	my @accToPcntId;
	$accToPcntId[0]->{hashname} = 'ID*len';      # %ID reported by blast, * $fLen
	$accToPcntId[1]->{hashname} = 'fixed';    # %ID*$fLen, with adjustments for XXXXXXX (low-complexity) and ----- (gap)
	$accToPcntId[2]->{hashname} = 'disorder'; # Like 'fixed', but also with disorder adjustment
	my %accToProt;
	my %accIsRight;   # $accIsRight{$acc} = 1 if &member($accToProt{$acc}, @equivalogNames)
	my %accToEvRef;   # $accToEvRef{$acc} = pointer to EvidenceBER

	# Read btab file, record %ID, get & record names, judge whether right or not
	&getNames($btab_file, \%accToProt, \%accIsRight, \%accToEvRef, \@accToPcntId, @equivalogNames);

	my @acc_arr = keys(%accToProt);
	# Sort in decreasing order of %ID, just for debugging, so generally good hits are first
	my @acc_arr = sort { $accToPcntId[0]->{$b} <=> $accToPcntId[0]->{$a} } @acc_arr;
	# Reduce size of list; no need to look at very bad hits
	@acc_arr = @acc_arr[0 .. 99];

	my $blastfile = "$blastfile_path/$orf.blastp";	
  	open(READFILE, $blastfile) or die "Can't open BLASTFILE $blastfile";
	print "Parsing alignments in $blastfile\n" if $DEBUG;
	my $alignmentsRef = &alignmentsHash($blastfile, 2000);
	#my @acc_arr = keys(%$alignmentsRef);
	#print "Acc @acc_arr\n";
	my $orfCorSum = 0;      # sum of all correlations for this orf
	my $orfResidues = 0;    # sum of residues in alignments for all hits for this orf
	my $orfCorSumC = 0;     # sum of all correlations for this orf for correct accessions
	my $orfResiduesC = 0;   # sum of residues in alignments for all hits for this orf for correct accessions
	# Compute probability of each residue
	my $residues = &listSum(values %ResidueCounts);
	my $residuesD = &listSum(values %ResidueCountsD);
	my (%residueP, %residueDP);
	foreach my $r (keys %ResidueCounts) {
		$residueP{$r} = $ResidueCounts{$r} / $residues;
		$residueDP{$r} = $ResidueCountsD{$r} / $residuesD;
	}

	$numOrfsTried++;
	open(CORFILE, ">$orf_cor") or die "Can't open $orf_cor for writing: $!";
	if ($DEBUG > 1) { @acc_arr = @acc_arr[0 .. 49] }
	print "Processing " . scalar(@acc_arr) . " accessions for $orf\n" if $DEBUG;
	my $hitsDir = "$blastfile_path/$orf";
	mkdir($hitsDir) if !-d $hitsDir;
	my $fastaDir = "$hitsDir/fasta";
	mkdir($fastaDir) if !-d $fastaDir;
  	foreach my $acc (@acc_arr) 
	{
		next if $acc !~ m/^SP/;
		my $hashRef = $alignmentsRef->{$acc};
		my $qstring = $hashRef->{query_string};
		if (length($qstring) < 10) {
			print "Skipping $acc; length of query string is only \n" . length($qstring);
			next;
		}
		my $hstring = $hashRef->{hit_string};
		if (length($hstring) < 10) {
			print "Skipping $acc; length of hit string is only \n" . length($hstring);
			next;
		}

		$acc =~ m/^SP\|(.*)$/; my $accid = $1;
		my $hitFasta = "$hitsDir/fasta/$accid.fasta";
		# Run IUPRED if it has not been run on blast hit $acc
		my $hitIuPrefix = "$hitsDir/iupred/$accid.iupred";
		my $fasta = 0;
		foreach my $suffix ('short', 'long') {
			my $iupFile = "$hitIuPrefix.$suffix";
			if (-e $iupFile) {
				$fasta = 1;
			}
			else {
				if (!-e $hitFasta) {
					# Create FASTA for hit
					$fasta = &writeFasta($acc, $hitFasta);  # returns 0 if failure
					if (!$fasta) {
						print "ERR: Could not create $hitFasta\n";
						unlink($hitFasta);
						next;
					}
				}

				# Call script to generate $iupFile
				my $cmd = "$DisoDir/iupred/iupred $hitFasta $suffix > $iupFile";
				my @output = `$cmd`;
			}
		}
		# Failure to find or create FASTA file
		if (!$fasta) {
			print "ERR: No $hitFasta\n";
			next;
		}

		# Parse $iupFile
		my @hitIupPosn;        # @{$iup_positions[1]} stores iupred result residue positions (long file)
		my @hitIupProb;        # Stores the disorder probability of iupred residues
		# Fill in @iupPosn, @iupProb, @iupResidue
		#my @hiupResidue = &parseIup($hitIuPrefix, \@hitIupPosn, \@hitIupProb);
		#my ($posRef, $resRef, $probRef) = &parseIup($anchor_file);

		#Run ANCHOR if it has not been run on blast hit $acc
		my $hitAnchorFile = "$hitsDir/anchor/$accid.anchor";
		if (!-e $hitAnchorFile) {
			my $cmnd = "perl $DisoDir/ANCHOR/anchor.pl $hitFasta -m $motifs_file > $hitAnchorFile";
			my @outnow = `$cmnd`;
		}

		#my $evRef = $accToEvRef{$acc};
		my $homolstring = $hashRef->{homology_string};
		my $qstart = $hashRef->{query_start};   # $qstart == 1 corresponds to $blast_val[0]
		my $qend = $hashRef->{query_end};
		my @queryChars = split('', $qstring);  
		my $hstart = $hashRef->{hit_start};
		my $hend = $hashRef->{hit_end};
		my $hgaps = $hashRef->{hit_gaps};
		#my $hlen = $evRef->len();
		my $hlen = $hashRef->{hit_length};  # Length of hit protein, NOT length of alignment on hit
		my $hitAlignLen = $hashRef->{hit_aligned};  # Length of hit protein in alignment
		my $fracId = $hashRef->{frac_identical};
		my $fracConserved = $hashRef->{frac_conserved};
		# Test on ann_test 1 showed that, for a gene with %ID=$id, of the unmatched positions,
		#   .876*$id were half-matches.
		my $expectedConserved = .876 * $fracId;
		my $fLen;
		if ($GapInLen) { $fLen = ($hitAlignLen - $hgaps) / $hlen }
		else { $fLen = ( $hitAlignLen ) / $hlen }
		die "ERR: fLen=$fLen" if $fLen < 0 || $fLen > 1;
		my @hitChars = split('', $hstring);
		my @homo_string = split('',$homolstring);
		my $homo_len = scalar(@homo_string);
		my $gapCount = 0;        # gap count
		my ($val, $valD);
		my $p;    # probability of match
		my $pD;   # probability of match in disordered region
		my @orderScores = ();
		my @blast_val = ();
		my $blastSum = 0;
		my $blastDSum = 0;
		my $orderSum = 0;
		my $conserved = 0;          # number of '+'
		my $blastTimesOrderSum =0;
		for(my $j=0; $j<scalar(@queryChars); $j++)
		{
			# Could have used frac_identical, frac_conserved, gaps instead
			# Compute distribution-based ID score
			my $queryChar = uc($queryChars[$j]);
			my $hitChar = uc($hitChars[$j]);

			if ($queryChar eq '-')
			{
				# Add gap penalty?
				next if $GapInLen;   # 'next' means "don't count gaps against the score'
				$gapCount++;         # used to match @blast_val to @orderProb
			}
			if ($hitChar eq '-')
			{
				next if $GapInLen;
			}
			next if ($queryChar eq 'X');
			next if ($hitChar eq 'X');

			#if($homo_string[$j] eq " ") { $val = 0; }
			if($homo_string[$j] eq "+") {
				my $p = $residueP{$queryChar};
				my $q = $residueP{$hitChar};
				$val =0.4 - $p*$q;
				my $pD = $residueDP{$queryChar};
				my $qD = $residueDP{$hitChar};
				$valD = 0.4 - $pD*$q;
				$conserved++;
			}
			elsif($homo_string[$j] =~ m/[A-Z]/) {
				my $p = $residueP{$queryChar};
				$val = 1 - $p*$p;
				my $pD = $residueDP{$queryChar};
				$valD = 1 - $pD*$p;
				#$blastSum += $val;
			}
			#else { die "ERR: Bad character in homology string: $homo_string[$j]" }
			else { $val = 0 }
			push(@blast_val, $val);
			my $index = $qstart+$j-($gapCount+1);
			#print "Query $queryChars[$j] \t Hit $hitChars[$j]\tHomo=$homo_string[$j]\tVal=$val\tIUPRED_Posn=$iupred_posn\n";
			my $order = $orderProb[$index];
			push(@orderScores, $order);
			die "ERR: Residue $queryChar reported in alignment[$j] with $acc does not match residue $anchor_res[$index] in anchor file at position $index+1"
					if ($queryChar ne '-') && ($queryChar ne $anchor_res[$index]);

			$orderSum += $order;
			$blastSum += $val;
			$blastDSum += ($order > .5) ? $val : $valD;
			$blastTimesOrderSum += ($order > .5) ? ($val * $order) : ($valD * (1-$order));
		} # for(my $j=0; $j<scalar(@queryChars); $j++)

		my $alignLen = scalar(@blast_val);
		if ($DEBUG) {
			#@blast_val=(1,3,3,5);
			#@orderScores=(12,12,11,7);
			my $orderLen = scalar(@orderScores);
			die "ERR: scalar(\@orderScores)=$orderLen <> scalar(\@blast_val)=$alignLen" if $alignLen != $orderLen;
		}

		# Test on ann_test 1 showed that, for a gene with %ID=$id, of the unmatched positions,
		#   .876*$id were half-matches.
		my $pcntId = $accToPcntId[0]->{$acc};
		my $expectedConserved = ($alignLen - $blastSum) * .876 * $pcntId;
		# Add half of excess over expected number of conserved sites to blastSum
		#$blastSum += .5 * ($conserved - $expectedConserved);
       
		$accToPcntId[0]->{$acc} *= $fLen;
		$accToPcntId[1]->{$acc} = $fLen * $blastSum/$alignLen;
		$accToPcntId[2]->{$acc} = $fLen * $blastDSum/$alignLen;
		$accToPcntId[3]->{$acc} = $fLen * $blastTimesOrderSum/$alignLen;
		if ($DEBUG) {
			my @ids;
			foreach my $ref (@accToPcntId) { push( @ids, &precision($ref->{$acc}, 4) ); }
			print "$acc: qs=$qstart qe=$qend bs=$blastSum al=$alignLen hlen=$hlen halen=$hitAlignLen hgaps=$hgaps flen=$fLen ids=@ids\n";
		}

		#$ psf{%len,%ID}:
		# %ID:    >90       >80       >70       >60       >50       >40       >30       >25   
		# >90  0.99959   0.99978   0.98493   0.91286   0.78210   0.73260   0.68310   0.63360   
		# >75  1.00000   0.99786   0.95033   0.81920   0.74260   0.69560   0.64860   0.60160   
		# >60  1.00000   0.99888   0.87532   0.58692   0.70310   0.65860   0.61410   0.56960   
		# >45  0.83160   0.78960   0.74760   0.70560   0.66360   0.62160   0.57960   0.53760   
		# >30  0.78210   0.74260   0.70310   0.66360   0.62410   0.58460   0.54510   0.50560   
		# >15  0.73260   0.69560   0.65860   0.62160   0.58460   0.54760   0.51060   0.47360   


		#Compute correlation between @blast_val[$qstart .. $]
		my $sd2;
		my $cov;
		my $mean1 = $blastSum / $alignLen;
		my $mean2 = $orderSum / $alignLen;
		#next if ($mean2 == 1);   # Not interested in completely-ordered proteins
		my $sd1 = &SD($mean1, @blast_val);
		my $sd2 = &SD($mean2, @orderScores);
		my $cov_sum=0;
		my $cor;
        
		for(my $d=0; $d<scalar(@blast_val); $d++)
		{
			$cov_sum += (($blast_val[$d]-$mean1)*($orderScores[$d]-$mean2));
		}
		$cov=$cov_sum/$alignLen;
		if(($sd1==0) || ($sd2 ==0)) {
			$cor =0;
			next;    # skip this accession
		}
		else { $cor=$cov/($sd1*$sd2); }
       
		$mean1 = &precision($mean1, 3);
		$mean2 = &precision($mean2, 3);
		$sd1 = &precision($sd1, 3);
		$sd2 = &precision($sd2, 3);
		$cor = &precision($cor, 5); 
		my $pcntId = &precision($accToPcntId[0]->{$acc}, 5);
		my $pcntIdFixed = &precision($accToPcntId[1]->{$acc}, 5);
		my $pcntIdDisorder = &precision($accToPcntId[2]->{$acc}, 5);
		print "$acc: alignLen=$alignLen blast(mean=$mean1 stdv=$sd1) order(mean=$mean2 stdv=$sd2) correct=$accIsRight{$acc} corr=$cor\n" if $DEBUG;
		print CORFILE "$orf\t$acc\t$cor\t$pcntId\t$pcntIdFixed\t$pcntIdDisorder\t$accIsRight{$acc}\t$accToProt{$acc}\n";
		$orfCorSum += $cor * $alignLen;
		$orfResidues += $alignLen;
		if ($accIsRight{$acc}) {
			$orfCorSumC += $cor * $alignLen;
			$orfResiduesC += $alignLen;
		}
	} # foreach my $acc (@acc_arr) 

	if ($orfResidues > 0) {
		my $orfCor = $orfCorSum / $orfResidues;   # average correlation value for this gene
		print "$orf: average correlation=$orfCor\n";   
		print CORFILE "$orf: average correlation=$orfCor\n";   
		$totalCorSum += $orfCorSum;
		$totalResidues += $orfResidues;
		if ($orfResiduesC) {
			$orfCor = $orfCorSumC / $orfResiduesC;
			print "$orf: average correlation when correct=$orfCorSumC/$orfResiduesC=$orfCor\n";
			$totalCorSumC += $orfCorSumC;
			$totalResiduesC += $orfResiduesC;
		}

		if (@equivalogNames) {
			# @compare becomes a list with the score for each method in @accToPcntId
			my @compare = &compareMethods($orf, \%accIsRight, \%accToProt, \@accToPcntId, \@equivalogNames, @acc_arr);
			print "$orf compareMethods=(@compare)\n";
			if ($compare[0] > 0) {
				print "Equivalog names:\n";
			  foreach my $en (@equivalogNames) { print "  $en\n" }
			}
			push(@comparisons, @compare);
			for (my $i=0; $i<@compare; $i++) {
				for (my $j=0; $j<@compare; $j++) {
					$CompareStats[$i]->[$j]++ if $compare[$i] > ($compare[$j] + 1);
				}
			}
		}
	}
	close(CORFILE) or die "Can't close $orf_cor: $!";
	print "\n";
	&printCompareStats() if !($numOrfsTried % 5) || $DEBUG;
} # Closing of loop foreach my $orf (@orfs)

if($totalResidues > 0) {
	my $fDiso = $queryResiduesD / $queryResidues;
	print "Disordered residues in query proteins: $queryResiduesD / $queryResidues = $fDiso\n";
	#print "Fraction of genes with 11 disordered residues in a row = $numDisorderedGenes / $numWiEquiv\n";
	my $fGDiso = $numDisorderedGenes / @orfs;
	print "Fraction of genes with 11 disordered residues in a row = $numDisorderedGenes / " .
	       scalar(@orfs) . " = $fGDiso\n";

	my $aveCor = $totalCorSum / $totalResidues;
	my $msg = "Average blast vs max(anchor, 1-iupred) correlation, over $numOrfsTried genes, $totalResidues residues = $aveCor";
	print "$msg\n";
	print CORFILE "$msg\n\n";
	if($totalResiduesC > 0) {
		my $aveCorC = $totalCorSumC / $totalResiduesC;
		$msg = "Average correlation when correct over $totalResiduesC residues = $aveCorC";
		print "$msg\n\n";
		print CORFILE "$msg\n\n";
	}
	&printCompareStats();

	my $entropy = &entropy(values %ResidueCounts);
	print "Entropy of all residues=$entropy\n";
	$entropy = &entropy(values %ResidueCountsD);
	print "Entropy of disordered residues=$entropy\n";

	if ($Write) {
		&hashToTsv($ResidueFile, \%ResidueCounts);

		&hashToTsv("${ResidueFile}D", \%ResidueCountsD);
	}
}

$UniprotH->disconnect();
exit(0);

#################################################################################


# Read btab file, record %ID, get & record names, judge whether right or not
sub getNames {
	my ($btab_file, $accToProtRef, $accIsRightRef, $accToEvRefRef, $accToPcntIdRef, @equivalogNames) = @_;

	#Reading data from BTAB file
	open(my $BTABFILE, $btab_file)or die "Can't open btab file";
	my @btab_data = <$BTABFILE>;
	close($BTABFILE);

	#Parsing BTAB output file
	my $queryUniprot = "SELECT name FROM accession a, info i WHERE a.accession_db=? AND a.accession_id=? AND a.id = i.id";
	print "Getting names for " .scalar(@btab_data) ." BLAST hits\n" if $DEBUG;
	for(my $b=0; $b<scalar(@btab_data); $b++)
	{
	 	my $btab_line = $btab_data[$b];
		chomp($btab_line);
		my @btab_info = split('\t',$btab_line);
		my $pcntId = $btab_info[BT_PCNT_IDENT];
		my $line_mod = $btab_info[BT_COMMENT];
		my @prot = split(/\^\|\^/, $line_mod);
		my $btab_acc = $btab_info[BT_ACCESSION];
		#my $qlen = $btab_info[BT_QUERY_LEN];   # Wrong, due to praze
		#my $hitLen = $btab_info[BT_HIT_LEN];
		#my $hitAlignLen = $btab_info[BT_HIT_STOP] + 1 - $btab_info[BT_HIT_START];
		if($btab_acc=~/^([A-Z]+\|[^\|]+)\|/)
		{
			my $acc = $1;
			next if ($acc !~ m/^SP/);   # Only want SP hits, so annotation is reliable
			# Skip blast hits with very high %ID
			next if $pcntId > 50;   # Is this okay, or does it disadvantage straight-up %ID?
			$accToPcntIdRef->[0]->{$acc} = $pcntId;   # updated later
			# Creating object just to hold onto hit length info
			#my $evRef = new Anno::EvidenceBER($acc);
			#$evRef->initFromBtab(@btab_info);
			#$accToEvRefRef->{$acc} = $evRef;

			# Get protein name from Uniprot
			my ($accdb, $accid) = split(/\|/, $acc);
			my @uniAccs = &myDoSql($UniprotH, $queryUniprot, $Delimiter, $accdb, $accid);
			my $name;
			if (@uniAccs) {
				$name = $uniAccs[0];
			}
			else {
				# Get protein name from Panda
				$name = $prot[0];
				$name =~ s/\{.*\}.*$//g;
			}

			next if (&uselessName($name) > 65);  # Skip if bad name

			my $name2 = &thesaurusName($name, \%Dictionary, $GSHashRef, $PNUHashRef);
			my $name3 = $name2;    # in case !exists($StrippedToStripped{$name2})
			$name3 = $StrippedToStripped{$name2} if exists($StrippedToStripped{$name2});
			$name3 = &stripQualifiers4($name3);
			next if $name3 eq '';
			print "$acc: '$name' => '$name2' => '$name3'\n" if $DEBUG;
			$accToProtRef->{$acc}=$name3;
			print "ERR: $acc name=$name3\n" if $name3 =~ m/^\d+$/;
			# Judge hit $acc as right if its name matches an equivalog's name
			$accIsRightRef->{$acc} = &member($name3, @equivalogNames);
		}
		else {
			die "Bad btab accession: $btab_acc"
		}
	}
}


sub printCompareStats {
	print    "CompareStats:\n\n";
	print    " Loser:    0    1    2\n";
	print    "Winner\n";
	for (my $i=0; $i<@CompareStats; $i++) {
		print "   $i:  ";
		for (my $j=0; $j<@CompareStats; $j++) {
			my $wins = &rightJustify( $CompareStats[$i]->[$j], 5 );
			print $wins;
		}
		print "\n";
	}

	my @rightProb;
	foreach my $r (@RightCount) {
		my $rp = &precision($r / $RightTries, 4);
		push(@rightProb, $rp);
	}
	print "Probability of being right by method (out of $RightTries): (@rightProb)\n\n";
}

# Convert iupred short, iupred long, & anchor scores into an order score in [0 .. 1]
sub getOrder {
	my ($iuProb, $anchor) = @_;
	die "ERR: anchor=$anchor" if ($anchor < 0 || $anchor > 1);
	my $thresholdedAnchor = &minThresh($anchor, .0);
	#die "ERR: thresholdedAnchor=$thresholdedAnchor" if ($thresholdedAnchor < 0 || $thresholdedAnchor > 1);
	my $thresholdedIupred = 1 - &minThresh($iuProb, .0);
	#die "ERR: thresholdedIupred=$thresholdedIupred" if ($thresholdedIupred < 0 || $thresholdedIupred > 1);
	my $orderScore = &max($thresholdedAnchor, $thresholdedIupred);
	return $orderScore;
}

# See which %ID method scores higher for this orf
# @$hashesRef is a list of hash refs
# @$equivNamesRef is list of equivalog names from the array returned by parseHtab
sub compareMethods {
	my ($orf, $accIsRightRef, $accToProtRef, $hashesRef, $equivNamesRef, @acc_arr) = @_;
	my @equivalogNames = @$equivNamesRef;
	my @scores;   # $scores[$x] = scores using %accToPcntId[$x]

	my $i = 0;
	foreach my $hashRef (@$hashesRef) {
		# Construct a list of the blast hit info lines, sorted by %$hashRef in decreasing order
		my @accs = sort { $hashRef->{$b} <=> $hashRef->{$a} } @acc_arr;

		# Compare name of top hits to correct name
		my $score = 0;    # score for correct names
		foreach my $blastAcc (@accs[0 .. 19]) {
			$score *= 1.3;
			my $correct = $accIsRightRef->{$blastAcc};
			$score += $correct;
		}
		if ($score > 0) {
			print "compare $hashRef->{hashname}: $score\n";
			foreach my $acc (@accs[0 .. 19]) {
				my $id = &precision($hashRef->{$acc}, 3);
				print "  $acc	$id $accToProtRef->{$acc}\n"
			}
		}
		push(@scores, $score);
		# Also update global (sorry!) @RightCount whether top blast hit was right
		$RightCount[$i++] += $accIsRightRef->{$accs[0]};
	}
	$RightTries++;

	return @scores;    # evals to 0 if @equivalogs == 1
}