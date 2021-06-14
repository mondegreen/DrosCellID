use FileHandle;

my $sam = shift;
my $r1 = shift;
my $r2 = shift;
my $outbase = shift;

my $fh;
if ($sam =~ /.gz$/) {
	$fh = new FileHandle("gunzip -c $sam|");
} else {
	$fh = new FileHandle($sam);
}
while (<$fh>) {
	chomp;
	my @d = split /\t/;
	next unless @d > 8;
	next if $d[1] == 4;
	$m{$d[2]}++;
	$mm{$d[0]} = $d[2];
}
$fh->close();

my $fhr1;
if ($r1 =~ /.gz$/) {
	$fhr1 = new FileHandle("gunzip -c $r1|");
} else {
	$fhr1 = new FileHandle($r1);
}
my %files;
foreach my $i (keys %m) {
	my $fhout = new FileHandle("| gzip -c > $outbase.$i.R1.fastq.gz");
	$files{$i} = $fhout;
}
while (<$fhr1>) {
	if (/^\@(\S+)/) {
		my $fho;
		if (exists($mm{$1})) {
			$fho = $files{$mm{$1}};
		} else {
			$fho = undef;
		}
		if (defined($fho)) {
			print $fho $_;
		}
		for (my $i = 0; $i < 3; $i++) {
			$_ = <$fhr1>;
			if (defined($fho)) {
				print $fho $_;
			}
		}
	}
}
foreach my $i (keys %m) {
	my $fhout = $files{$i};
	$fhout->close();
}
$fhr1->close();

my $fhr2;
if ($r2 =~ /.gz$/) {
	$fhr2 = new FileHandle("gunzip -c $r2|");
} else {
	$fhr2 = new FileHandle($r2);
}
my %files;
foreach my $i (keys %m) {
	my $fhout = new FileHandle("| gzip -c > $outbase.$i.R2.fastq.gz");
	$files{$i} = $fhout;
}
while (<$fhr2>) {
	if (/^\@(\S+)/) {
		my $fho;
		if (exists($mm{$1})) {
			$fho = $files{$mm{$1}};
		} else {
			$fho = undef;
		}
		if (defined($fho)) {
			print $fho $_;
		}
		for (my $i = 0; $i < 3; $i++) {
			$_ = <$fhr2>;
			if (defined($fho)) {
				print $fho $_;
			}
		}
	}
}
foreach my $i (keys %m) {
	my $fhout = $files{$i};
	$fhout->close();
	print STDERR "$sam\t$i\t$m{$i}\n";
}
$fhr2->close();
