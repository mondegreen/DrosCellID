use FileUtil;

### Takes as input:
### 1) file containing a tab-delimited listing of the:
### - file names of the results
### - the sample name they should be assigned to
### - transposon class being surveyed
###
### data is counted and converted into a tab-delimited counts files with samples in the header and rows indication the type of transposon, the molecule it is associated with, and the bp positions of the window that was hit
###
### for example a row id consists of 4 period-separated values:
### 1731.2L.11581310.11581902
### 
### 1) 1731 is the transposon id
### 2) 2L is the molecule or chromosome
### 3) 11581310 is the begin coordinate
### 4) 11581902 is the end coordinate

my $list = shift;
my $lfh = FileUtil::openFileHandle($list);
while (<$lfh>) {
	chomp;
	my @d = split /\t/;
	my $file = $d[0];
	my $id = $d[1];
	my $te = $d[2];
	my $fh = FileUtil::openFileHandle($file);
	while (<$fh>) {
		chomp;
		my @d = split /\t/;
		next if $d[1] == 4;
		my @e = split /\./,$d[2];
		next unless $e[0] eq $te;
		$m{$id}{$d[2]}++;
		$mm{$d[2]}++;
	}
	$fh->close();
}

my @l;
foreach my $i (sort keys %m) {
	push @l,$i;
}
my $l = join "\t","",@l;
print $l,"\n";

foreach my $i (sort keys %mm) {
	my @l;
	foreach my $j (sort keys %m) {
		push @l,$m{$j}{$i}+0;
	}
	my $l = join "\t",$i,@l;
	print $l,"\n";
}
