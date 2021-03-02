#!/usr/bin/perl -w
# script to automatically unlock features against each other
#
# Usage: perl unlock_features.pl --feature <myfeature> src/io/parameter/check_parameter.cc
#
# if --feature is given, it unlocks all features given in the file against the given feature.
# if --feature is omitted, it unlocks all features against each other.
#
# by Nathan Schmid, May 2008

use strict;
use Getopt::Long;

my $feature = '';
if (!GetOptions('feature=s' => \$feature)) {
  print STDERR "Usage: perl $0 --feature <f> src/io/parameter/check_parameter.cc\n";
  exit(1);
}

my @available_features = ();

# read and parse the file
while(my $line = <>) {
  while($line =~ m/(?:^|\s)add\(\w*"([^"]+)",/g) {
    push @available_features, $1;
  }
}

if ($feature) {
  if (not grep { $_ eq $feature; } @available_features) {
    print STDERR "Feature $feature is not available\n";
    exit(1);
  }

  foreach my $f (@available_features) {
    next if ($f eq $feature);
    print <<EOF;
  fc.unlock("${feature}", "${f}");
EOF
  }
} else {
  for(my $i = 0; $i <= $#available_features; $i++) {
    for(my $j = $i + 1; $j <= $#available_features; $j++) {
      my $f1 = $available_features[$i];
      my $f2 = $available_features[$j];
      print <<EOF;
  fc.unlock("${f1}", "${f2}");
EOF
    }
  }
}
