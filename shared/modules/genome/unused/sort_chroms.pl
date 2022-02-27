use strict;
use warnings;

# action:
#     ensure all chromosome entries are sorted as chr1..chr22,chrX,chrY[,chrM]
#     implicitly filters by chromosome at the same time
#     does NOT sort by position, only chrom; position sort is the same as the input stream
# optional arguments:
#     chromCol useChrM
# usage:
#     call in a stream as ... | perl $MODULES_DIR/genome/sort_chroms.pl | ...
#     to use a column index other than column 1 for chromosome:
#        perl $MODULES_DIR/genome/sort_chroms.pl 3
#     to include chrM, include chromosome column and set useChrM:
#        perl $MODULES_DIR/genome/sort_chroms.pl 1 TRUE
# caution:
#     this script will take a lot of memory for large input streams!

# optional boolean to include chrM (omitted by default)
my ($chromCol, $useChrM) = @ARGV;
$chromCol or $chromCol = 1; # default to BED-type format
$chromCol--;
my $nSplit = $chromCol + 2;

# working variables
my @useChroms = map { "chr$_" } (1..22, "X", "Y", $useChrM ? "M" : ());
my %useChroms = map { $_ => 1 } @useChroms;

# filter the stream, collect all lines in memory
my %lines;
while(my $line = <STDIN>){
    chomp $line;
    my @f = split("\t", $line, $nSplit);
    $useChroms{$f[$chromCol]} or next;
    push @{$lines{$f[$chromCol]}}, "$line\n";
}

# print back out in natural sort order
foreach my $chrom(@useChroms){
    $lines{$chrom} and print join("", @{$lines{$chrom}}); # retains newlines
}
