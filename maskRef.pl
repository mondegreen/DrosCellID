use FileHandle;

my $fasta = shift; # genome to be masked
my $blast = shift; # blast tab-delimited output
my $ffh = new FileHandle($fasta);
my $bfh = new FileHandle($blast);
my %s;
my $id;
while (<$ffh>) {
	if (/^>(\S+)/) {
		$id = $1;
	} else {
		chomp;
		s/\s+//g;
		push @{$s{$id}},$_;
	}
}
$ffh->close();

my %ss;
foreach my $i (keys %s) {
	$ss{$i} = join "",@{$s{$i}};
}

while (<$bfh>) {
	chomp;
	my @d = split /\t/;
  ($d[8],$d[9]) = $d[8] < $d[9] ? ($d[8],$d[9]):($d[9],$d[8]);
	my $x = substr($ss{$d[1]},$d[8],$d[9]-$d[8]);
	my $n = "N" x length($x);
	substr($ss{$d[1]},$d[8],$d[9]-$d[8]) = $n;
}
$bfh->close();

foreach my $i (keys %s) {
	print ">$i\n$ss{$i}\n";
}
