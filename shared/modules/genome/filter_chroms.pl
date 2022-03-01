use strict;
use warnings;

# action:
#     remove all chromosome entries from input except chr1..chr22,chrX,chrY[,chrM]
# optional arguments:
#     chromCol useChrM
# usage:
#     call in a stream as ... | perl $MODULES_DIR/genome/filter_chroms.pl | ...
#     to use a column other than column 1 for chromosome:
#        perl $MODULES_DIR/genome/filter_chroms.pl 3
#     to include chrM, include chromosome column and set useChrM:
#        perl $MODULES_DIR/genome/filter_chroms.pl 1 TRUE

# optional boolean to include chrM (omitted by default)
my ($chromCol, $useChrM) = @ARGV;
$chromCol or $chromCol = 1; # default to BED-type format
$chromCol--;
my $nSplit = $chromCol + 2;

# working variables
my @useChroms = map { "chr$_" } (1..22, "X", "Y", $useChrM ? "M" : ());
my %useChroms = map { $_ => 1 } @useChroms;

# filter the stream
while(my $line = <STDIN>){
    chomp $line;
    my @f = split("\t", $line, $nSplit);
    $useChroms{$f[$chromCol]} or next;
    print "$line\n";
}

