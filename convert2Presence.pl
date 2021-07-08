use FileUtil;

### Takes as input:
### 1) file containing a tab-delimited listing of the:
### - file names of the results
### - the sample name they should be assigned to
### - transposon class being surveyed
###
### along with 2 optional parameters:
### 2) the minimum number of positives per row
### 3) the maximum number of positives per row
###
### the files are read in, and junctions that are found in a sample are marked as found
### then a tab-delimited table is generated.
### the first row contains a header with all the sample names
### the subsequent rows contain a location (molecule.position) followed by a binary value
### indicating whether that particular sample was found or not
### 
### by default rows with transposons found in all samples or found only in less than 2 samples
### are excluded

my $mappingFile = shift;
my $min = shift;
my $max = shift;

$min = 2 unless defined($min);

my %sample;
my %type;
my $mfh = FileUtil::openFileHandle($mappingFile);
while (<$mfh>) {
	chomp;
	my @d = split /\t/;
	$sample{$d[0]} = $d[1];
	$type{$d[0]} = $d[2];
}
$mfh->close();

foreach my $file (keys %sample) {
	my $fh = FileUtil::openFileHandle($file);
	while (<$fh>) {
		chomp;
		my @d = split /\t/;
		my $pos = join ".",$d[0],sprintf "%9.9d",$d[1];
		$col{$sample{$file}} = 1;
		$m{$pos}{$sample{$file}} = 1;
		if (@d > 6) {
			$pos = join ".",$d[9],sprintf "%9.9d",$d[10];
			$col{$sample{$file}} = 1;
			$m{$pos}{$sample{$file}} = 1;
		}
	}
	$fh->close();
}

my @l;
foreach my $i (sort keys %col) {
	push @l,$i;
}

if (defined($max)) {
} else {
	$max = @l;
}

### print sample names in header
my $l = join "\t","",@l;
print $l,"\n";

foreach my $i (sort keys %m) {
	my @l;
	my $t = 0;
	foreach my $j (sort keys %col) {
		push @l,$m{$i}{$j}+0;
		$t += $m{$i}{$j};
	}
	next if $t < $min || $t >= $max;
	my $l = join "\t",$i,@l;
	print $l,"\n";
}
