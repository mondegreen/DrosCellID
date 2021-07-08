use FileUtil;

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
