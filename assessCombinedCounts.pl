use FileUtil;

### Takes as input of samples (columns) and transposon,position (rows) each containing the number of reads associated with that sample/tranposon,position
###
### Output is an anologous table with the values replaced with 0 or 1 depending on the z-score associated with each row (indicating absence/presence)
### 
### Note: Because a particular transposon can occur in every sample, during z-score analysis the row values are supplemented with an equal number of zero values
### before z-scores are calculated

my $f = shift;

my $fh = FileUtil::openFileHandle($f);
$_ = <$fh>;
chomp;
my @header = split /\t/,$_;
shift @header;
while (<$fh>) {
	chomp;
	my @d = split /\t/;
	my $id = shift @d;
	my $sum = 0;
	for (my $i = 0; $i < @d; $i++) {
		$sum += $d[$i];
	}
	next unless $sum > 0;
	$m{$id} = \@d;
}
$fh->close();

my @totals;
for (my $i = 0; $i < @header; $i++) {
	my $tot = 0;
	foreach my $j (keys %m) {
		$tot += ${$m{$j}}[$i];
	}
	$totals[$i] = $tot;
}

### normalize each samples to each other
my $mean = getMean(\@totals);
for (my $i = 0; $i < @header; $i++) {
	my $xfactor = $totals[$i]/$mean;
	foreach my $j (keys %m) {
		${$mm{$j}}[$i] = ${$m{$j}}[$i]/$xfactor;
	}
}

my @balance;
for (my $j = 0; $j < @header; $j++) {
	push @balance,0;
}
### z-score each row and convert to binary based on whether the score is positive or negative
foreach my $i (keys %m) {
	my @temp = @{$mm{$i}};
	push @temp,@balance;
	$mean{$i} = getMean(\@temp);
	$std{$i} = getMeanStandardDeviation($mm{$i});
	for (my $j = 0; $j < @header; $j++) {
		${$mmm{$i}}[$j] = (${$mm{$i}}[$j] - $mean{$i})/$std{$i};
		${$mmm{$i}}[$j] = ${$mmm{$i}}[$j] < 0 ? 0:1;
	}
}

### print binary table
my $l = join "\t","",@header;
print $l,"\n";
foreach my $i (keys %m) {
	my $l = join "\t",$i,@{$mmm{$i}};
	print $l,"\n";
}

sub getMean {
  my ($values) = @_;
  
  my $count = 0;
  my $sum = 0;
  foreach my $value (@{$values}) {
    $count++;
    $sum+= $value;
  }
  my $avg;
  if ($count) {
    $avg = $sum / $count;
  } else {
    
  }
  return ($avg);
}

sub getMeanStandardDeviation {
  my ($values,$mean) = @_;
  
  if (!defined($mean)) {
    $mean = getMean($values);
  }
  
  my $sumOfSquares = 0;
  my $count = 0;
  foreach my $value (@{$values}) {
    $sumOfSquares += ($value - $mean) * ($value - $mean);
    $count++;
  }
  my $std = 0;
  if ($count > 1) {
    $std = sqrt($sumOfSquares/($count - 1));
  }
  
  return ($std);
}
