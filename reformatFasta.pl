use FileUtil;

### This script takes a fasta file and reformats it in a variety of ways for easy viewing/processing.
### It automatically corrects for some common artifacts in fasta files like control-M's
###
### the inputs are:
### 1) the fasta file
### 2) optional: boolean indicating whether to reverse complement the fasta file (default if 0 or false)
### 3) optional: line length to format the fasta file (default 80)
### 4) optional: convert to uppercase (if true upper cases the fasta sequence; default false)

my @seq;
my $def;
my $file = shift;
my $revcomp = shift;
my $lineLength = shift;
my $upperCase = shift;
$revcomp = 0 unless defined($revcomp);
$lineLength = 80 unless defined($lineLength);
my $zero = 0;
my $fh = FileUtil::openFileHandle($file);
while (<$fh>) {
  chomp;
  if (/^(\#?)>(.*)/) {
    my $newMark = $1;
    my $newDef = $2;
    $newDef =~ s/\cM//g;
    if (defined($def)) {
      my $seq;
      if ($revcomp) {
        $seq = join "",@seq;
        $seq =~ tr/ACTGactg/TGACtgac/;
        $seq = reverse($seq);
      } else {
        $seq = join "",@seq;
      }
      $seq =~ s/\cM//g;
      my $len = length($seq);
      if ($len == 0) {
        if ($def =~ /\/full_length=(\S+)/) {
          $len = $1;
        }
      }
      my ($nlen) = $seq =~ tr/Nn/Nn/;
      my ($gc) = $seq =~ tr/gcGC/gcGC/;
      my $gcNum= $gc;
      if (defined($len) && $len) {
        if ($len-$nlen) {
          $gc = int(1000*$gc/($len-$nlen))/10;
        } else {
          $gc = 0;
        }
      }
      if (defined($len) && $len) {
        if ($def =~ s/(\/gc=)\S+/$1$gc/) {
        } else {
          $def .= " /gc=$gc";
        }
      }
      if ($def =~ s/(\/offset=)\S+/$1$zero/) {
      } else {
        $def .= " /offset=$zero";
      }
      if ($def =~ s/(\/full_length=)\S+/$1$len/) {
      } else {
        $def .= " /full_length=$len";
      }
      if ($def =~ s/(\/nlength=)\S+/$1$nlen/) {
      } else {
        $def .= " /nlength=$nlen";
      }
      if ($len-$nlen) {
        print "$mark>$def\n";
        if (length($seq)) {
          for (my $s = 0; $s < $len; $s += $lineLength){
            print substr($seq, $s, $lineLength),"\n";
          }
        }
      }
    }
    undef @seq;
    $def = $newDef;
    $mark = $newMark;
    next;
  }
  if ($upperCase) {
    $_ = uc($_);
  }
  push @seq,$_;
}
$fh->close();
if (defined($def)) {
  my $seq;
  if ($revcomp) {
    $seq = join "",@seq;
    $seq =~ tr/ACTGactg/TGACtgac/;
    $seq = reverse($seq);
  } else {
    $seq = join "",@seq;
  }
  $seq =~ s/\cM//g;
  my $len = length($seq);
  if ($len == 0) {
    if ($def =~ /\/full_length=(\S+)/) {
      $len = $1;
    }
  }
  my ($gc) = $seq =~ tr/gcGC/gcGC/;
  if (defined($len) && $len) {
    $gc = int(1000*$gc/($len-$nlen))/10;
  }
  my $nlen = $seq =~ tr/Nn/Nn/;
  if (defined($len) && $len) {
    if ($def =~ s/(\/gc=)\S+/$1$gc/) {
    } else {
      $def .= " /gc=$gc";
    }
  }
  if ($def =~ s/(\/offset=)\S+/$1$zero/) {
  } else {
    $def .= " /offset=$zero";
  }
  if ($def =~ s/(\/full_length=)\S+/$1$len/) {
  } else {
    $def .= " /full_length=$len";
  }
  if ($def =~ s/(\/nlength=)\S+/$1$nlen/) {
  } else {
    $def .= " /nlength=$nlen";
  }
  print "$mark>$def\n";
  if (length($seq)) {
    for (my $s = 0; $s < $len; $s += $lineLength){
      print substr($seq, $s, $lineLength),"\n";
    }
  }
}
