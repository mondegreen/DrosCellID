use FileUtil;

# input is a list of processed read alignments. Files should be grepped such that the filename is pre-pended to every line of data that is read.
# Example usage:
#   grep TEs outputDirectory/*txt | perl generateReadCountSummaryFile.pl - > summary.file.txt

my $list = shift;
my $lfh = FileUtil::openFileHandle($list);

my %seen;
while (<$lfh>) {
	chomp;
	s/.repeats//;
	s/:/\t/;
	my @d = split /\t/;
	my $x = $d[0];
	$d[0] =~ s/.*GSF\d+-//;
	my @e = split /\./,$d[0];
	my @f;
	push @f,$x,$e[0],$e[1],@d[1..13];
	my @g;
	push @g,$f[0],$f[2],$f[13],$f[14],$f[15];
	my $l = join "\t",@g;
	push @data,\@g unless exists($seen{$l});
	$seen{$l}++;
}
$lfh->close();

foreach my $i (sort { ${$a}[0] cmp ${$b}[0] || ${$a}[2] <=> ${$b}[2];} @data) {
	my $l = join "\t",@{$i};
	print $l,"\n";
}
