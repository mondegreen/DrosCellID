use FileUtil;

### this file takes a set of candidate TE junctions and tries to determine if they represent intact transposons
###
### The input is:
### 1) the set of candidate TE junctions
### 2) the transposon-masked reference genome
### 
### transposons that are already in the genome should be masked so if a pair of junctions adjecent on the genome in the proper orientation
### is largely composed of N's, that identifies a known TE insertion
###
### conversely, if the transposon is not in the genome, then the junctions should be correctly oriented but lie very close to one another
### in their genomic coordinates (in this case we assume <=10 bps)
###
### the output for unpaired TE junctions is tab-delimited and looks like:
### 2L	309595	31	9	0
###
### the columns are:
### 1) molecule
### 2) coordinate (0 indexed interbase coordinate system)
### 3) number of reads identified
### 4) number of independent starts identified in the forward genome
### 5) number of independent starts identified in the reverse-complemented genome
###
### the output for paired TE junctions is tab-delimited and looks like:
### 2L	498332	236	0	29	x	7	0	x	2L	498339	30	9	0
### 1) molecule
### 2) coordinate (0 indexed interbase coordinate system)
### 3) number of reads identified
### 4) number of independent starts identified in the forward genome
### 5) number of independent starts identified in the reverse-complemented genome
### 6) spacer
### 7) distance between second element
### 8) number of N's between second element
### 9) spacer
### 10) molecule
### 11) coordinate (0 indexed interbase coordinate system)
### 12) number of reads identified
### 13) number of independent starts identified in the forward genome
### 14) number of independent starts identified in the reverse-complemented genome

my $file = shift;
my $fasta = shift;

my $ffh = FileUtil::openFileHandle("perl reformatFasta.pl $fasta 0 999999999|");
while (<$ffh>) {
	if (/^>(\S+)/) {
		my $id = $1;
		$_ = <$ffh>;
		chomp;
		$s{$id} = $_;
	}
}
$ffh->close();

my $ffh = FileUtil::openFileHandle($file);
while (<$ffh>) {
	chomp;
	my @d = split /\t/;
	my $r;
	if (@l) {
		my $delta = abs($d[1]-$l[1]);
		if ($d[0] eq $l[0]) {
			$r = predict(\@l,\@d);
		} else {
			$l = join "\t",@l;
			print $l,"\n";
		}
	}
	if ($r) {
		undef @l;
		undef @d;
	} else {
		@l = @d;
	}
}
my $delta = abs($d[1]-$l[1]);
if (@l) {
	$l = join "\t",@l;
	print $l,"\n";
}

sub predict {
	my ($a,$b) = @_;
	
	my $result = 0;
	my @a = @{$a};
	my @b = @{$b};
	my $sub = substr($s{$a[0]},$a[1],$b[1]-$a[1]);
	my ($ns) = $sub =~ tr/Nn/Nn/;
	my $delta = $b[1]-$a[1];
	if (($ns >= (0.9 * $delta) || $delta <= 10) && (($a[3] && $b[4] && !$a[4] && !$b[3]) || ($a[4] && $b[3] && !$a[3] && !$b[4]))) {
		my $l = join "\t",@a,"x",$delta,$ns,"x",@b;
		print $l,"\n";
		$result = 1;
	} else {
		my $l = join "\t",@a;
		print $l,"\n";
		$result = 0;
	}
	return ($result);
}
