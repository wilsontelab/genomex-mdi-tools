use strict;
use warnings;

# action:
#     expand a bedgraph stream to include missing 0-value runs
#     some of these will typically be in genome gaps; filter those later as needed
#     only adds runs between existing runs (not at the beginning and end of chromosome)
# optional arguments:
#     nullValue [default; 0] 
# usage:
#     call in a stream as ...| perl $MODULES_DIR/genome/add_missing_runs.bedgraph.pl | ...

# working variables
my ($nullValue) = @ARGV; # set the value used for missing runs (typically 0)
$nullValue or $nullValue = 0;
my ($prevChrom, $prevEnd) = ('');

# run the bedgraph lines
while (my $line = <STDIN>) {
    $line =~ m/^track/ and next;
    my ($chrom, $start, $end) = split("\t", $line);
    $chrom ne $prevChrom and ($prevChrom, $prevEnd) = ($chrom);
    defined $prevEnd and $start > $prevEnd and print join("\t", $chrom, $prevEnd, $start, $nullValue), "\n";
    print $line;
    $prevEnd = $end;    
}
