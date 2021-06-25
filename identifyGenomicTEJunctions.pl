use FileUtil;

# This file takes the raw read data mapped to a reference genome and identifies informative reads (aka reads that show signs of containing a junction meeting certain minimal criteria given a locally aligned read)
#
# The inputs are:
#   - Expected read length
#   - Max junction error (so account for fuzzy edges)
#   - Max read error (max number of mismatches allowed in read)
#   - log file
#   - junction type (used to name output file)
#   - input file 1: the R1 read mapped to the genome (either forward or reverse complement)
#   - input file 2: the R1 read mapped to the TEs (file 1 and 2 are paired)
#
# The goal is to see the read align to the genome with some soft clipping, and then the same read should also map to a TE with the rest of the read being soft-clipped.
# These two local hits should more or less be exclusive of each other. If these criteria and that the read maps somewhat uniquely - it is output with the relevant info.
#
# The output file contains the following information in a tab-delimited format:
#   - read id
#   - genome molecule id 
#   - position on genome
#   - q-score on genome
#   - cigar string on genome
#   - "TEs"
#   - position on TE
#   - qscore on TE
#   - cigar on TE
#   - errors
#   - total good reads for this TE
#   - total seen
#   - total with bad junction
#   - total with bad reference

my $expectedLen = shift;
my $maxJunctionError = shift;
my $maxReadError = shift;
my $logFile = shift;
my $type = shift;
my @input = @ARGV;

my $key = "TEs";

foreach my $input (@input) {
  ### this step is clunky - it looks at the input file name to determine what kind of TE it is working with and whether the file represents reads aligned to the genome or TE elements
	if ($input =~ /gz.(TEs).gz$/|| $input =~ /.*R1\.gz\.(dmel).gz/ || $input =~ /.*R1\.gz\.(dmel).rc.gz/||$input =~ /.*sam\.(\w+)\.R1\..*\.gz/||$input =~ /.*R1\.(\w+).rc.sam.gz/||
			$input =~ /.repeats.sam.\w+.R1.gz.(dmel).gz/ || $input =~ /.repeats.sam.\w+.R1.gz.(dmel).rc.gz/||
			$input =~ /.*R1\.(\w+).sam.gz/) {
		$key = $1;
		my $ifh = FileUtil::openFileHandle($input);
		while (<$ifh>) {
			my @d = split /\t/;
			$seen{$d[0]} = 1;
			if (/NM:i:(\d+)/) {
				my $readError = $1;
				if ($key eq "TEs") {
					next unless $d[2] eq $type;
				}
				next if $readError > $maxReadError;
				next if exists($mol{$d[0]}{$key});  # next if id-key pair has already been seen
				if ($key ne "TEs" && length($d[2]) > 3) {
					$badRef{$d[0]} = 1;
					next;				
				} else {
				}
				$mol{$d[0]}{$key} = $d[2];
				$pos{$d[0]}{$key} = $d[3];
				$qscore{$d[0]}{$key} = $d[4];
				$cigar{$d[0]}{$key} = $d[5];
			}
		}
		$ifh->close();
	} else {
		print STDERR "ERROR: file $input doesn't have a reference indicator\n";
		exit;
	}
}

foreach my $i (keys %mol) {
	next unless exists($mol{$i}{"dmel"}) && exists($mol{$i}{$key});
	next if $cigar{$i}{$key} eq "*";
	my $tmp = $cigar{$i}{"dmel"};
	my %dmel;
	while ($tmp =~ s/(\d+)(\w)//) {
		$dmel{$2} = $1;
	}
	my %cop;
	my $tmp = $cigar{$i}{"TEs"};
	while ($tmp =~ s/(\d+)(\w)//) {
		$cop{$2} = $1;
	}
	my $error = abs($cop{S}-$dmel{M}) + abs($cop{M}-$dmel{S});
	if ($error > $maxJunctionError) {
		$badJ{$i} = 1;
	} else {
		$good{$i} = 1;
	}
}

my $totalBadJ = keys %badJ;
my $totalBadR = keys %badRef;
my $totalSeen = keys %seen;
my $totalGood = keys %good;
foreach my $i (keys %mol) {
	next unless exists($mol{$i}{"dmel"}) && exists($mol{$i}{$key});
	next if $cigar{$i}{$key} eq "*";
	my $tmp = $cigar{$i}{"dmel"};
	my %dmel;
	while ($tmp =~ s/(\d+)(\w)//) {
		$dmel{$2} = $1;
	}
	my %cop;
	my $tmp = $cigar{$i}{$key};
	while ($tmp =~ s/(\d+)(\w)//) {
		$cop{$2} = $1;
	}
	my $error = abs($cop{S}-$dmel{M}) + abs($cop{M}-$dmel{S});
	next if $error > $maxJunctionError;
	my $l = join "\t",$i,$mol{$i}{"dmel"},$pos{$i}{"dmel"},$qscore{$i}{"dmel"},$cigar{$i}{"dmel"},$key,$pos{$i}{$key},$qscore{$i}{$key},$cigar{$i}{$key},$error,$totalGood,$totalSeen,$totalBadJ,$totalBadR;
	print $l,"\n";
}
