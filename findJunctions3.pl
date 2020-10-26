use FileUtil;

my $counts = shift;
my $element = shift;
my $method = shift;
my $fwd1 = shift;
my $fwd2 = shift;
my $rc1 = shift;
my $rc2 = shift;
my $fasta = shift;
my $minCnt = shift;
my $minSide = shift;

$minCnt = 12 unless defined($minCnt);
$minSide = 4 unless defined($minSide);

print STDERR "$fwd1\n$rc1\n";

my $below = 0;
my $totalSam = 0;
my $minReads = 9999999999;
my $cfh = FileUtil::openFileHandle($counts);
while (<$cfh>) {
	my @d = split /\t/;
	$match{$d[0]} = $d[2];
	$reads{$d[0]} = $d[3];
#	print "$d[0] $match{$d[0]} $reads{$d[0]} $fwd1 $fwd2\n";
	if ($d[1] eq $element) {
		if ($method ne "minimum") {
			$below++ if $d[3] < $method;
		}
		$totalSam++;
		$minReads = $d[3] if $d[3] < $minReads;
	}
}
$cfh->close();

print STDERR $minReads,"\n";
if ($method eq "minimum") {
} else {
	if ($minReads < $method) {
		print STDERR "WARNING: For $element the minimum number of reads was $minReads but $below out of $totalSam samples were below the cutoff set at $method\n";
	}
	$minReads = $method;
}

my $ffh = FileUtil::openFileHandle($fasta);
while (<$ffh>) {
	if (/^>(\w+).*full_length=(\d+)/) {
		$len{$1} = $2;
	}
}
$ffh->close();

sub accumulateData {
	my ($f,$num) = @_;
#	print $num,"\n";
	my @r;
	my $fh = FileUtil::openFileHandle($f);
	while (<$fh>) {
		chomp;
		my @d = split /\t/;
		push @r,\@d;
	}
	$fh->close();
	my $seen = @r;
	if ($seen > $num) {
		my @result;
		my @seen;
		my $v = 0;
		while ($v < $num) {
			my $rnd = int(rand()*$seen);
			next if $seen[$rnd];
			$seen[$rnd] = 1;
			push @result,$r[$rnd];
			$v++;
		}
		return \@result;
	} else {
		return \@r;
	}
}

#print "$match{$fwd1} $minReads / $reads{$fwd1}\n";
my $data = accumulateData($fwd1,int($match{$fwd1}*$minReads/$reads{$fwd1}+0.5));
my $x = @{$data};
#print "xxx ",$x,"\n";
foreach my $i (@{$data}) {
	my @d = @{$i};
	my $cig = $d[4];
#	print $cig,"\n";
	$cig =~ s/^1S//;
	if ($cig =~ /^(\d+)M/) {
		$target = $d[2]+$1;
		$cand{$d[1]}{$target}++;
		$spread{$d[1]}{$target}{0}{$1}++;
	}
}

my $data = accumulateData($fwd2,int($match{$fwd2}*$minReads/$reads{$fwd2}+0.5));
foreach my $i (@{$data}) {
	my @d = @{$i};
	my $cig = $d[4];
#	print $cig,"\n";
	$cig =~ s/^1S//;
	if ($cig =~ /^(\d+)M/) {
		$target = $d[2]+$1;
		$cand{$d[1]}{$target}++;
		$spread{$d[1]}{$target}{0}{$1}++;
	}
}

my $data = accumulateData($rc1,int($match{$rc1}*$minReads/$reads{$rc1}+0.5));
foreach my $i (@{$data}) {
	my @d = @{$i};
	my $cig = $d[4];
	$cig =~ s/^1S//;
#	print $cig,"\n";
	if ($cig =~ /^(\d+)M/) {
		$target = $d[2]+$1;
		$target = $len{$d[1]} - $target;
		$cand{$d[1]}{$target}++;
		$spread{$d[1]}{$target}{1}{$1}++;
	}
}


my $data = accumulateData($rc2,int($match{$rc2}*$minReads/$reads{$rc2}+0.5));
foreach my $i (@{$data}) {
	my @d = @{$i};
	my $cig = $d[4];
	$cig =~ s/^1S//;
#	print $cig,"\n";
	if ($cig =~ /^(\d+)M/) {
		$target = $d[2]+$1;
		$target = $len{$d[1]} - $target;
		$cand{$d[1]}{$target}++;
		$spread{$d[1]}{$target}{1}{$1}++;
	}
}

foreach my $i (keys %cand) {
	next if length($i) > 3;
	foreach my $j (sort { $a <=>$b;} keys %{$cand{$i}}) {
		next if $cand{$i}{$j} < $minCnt;
		my $sprdLeft = scalar(keys %{$spread{$i}{$j}{0}});
		my $sprdRight = scalar(keys %{$spread{$i}{$j}{1}});
		next unless $sprdLeft >= $minSide || $sprdRight >= $minSide;
		my $pos = $j;
		print "$i\t$pos\t$cand{$i}{$j}\t$sprdLeft\t$sprdRight\n";
#		for (my $k = 60; $k > 20; $k--) {
#			my ($x,$y) = ($spread{$i}{$j}{0}{$k}+0,$spread{$i}{$j}{1}{$k}+0);
#			print "=\t$k\t$x\t$y\n";
#		}
	}
}
