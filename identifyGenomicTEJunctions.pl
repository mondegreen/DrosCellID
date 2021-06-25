use FileUtil;

# This file takes the raw read data mapped to a reference genome and identifies informative reads (aka reads that show signs of containing a junction meeting certain minimal criteria given a locally aligned read)

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
