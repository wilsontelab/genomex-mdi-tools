use strict;
use warnings;

# list canonical chromosomes

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);

# initialize the genome
use vars qw(@canonicalChroms);
setCanonicalChroms();
print join(" ", @canonicalChroms);
