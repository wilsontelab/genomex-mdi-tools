use strict;
use warnings;

# convert UCSC gc5Base.wigVarStep.gz to binned bed
# not all possible genome bins will be present
# no guarantee how much data actually contributed to a bin based on this output alone

#variableStep chrom=chr1 span=5
#10001   40
#10006   40
#10011   40
#10016   60
#10021   60
#10026   60
#10031   40
#10036   40
#10041   40

# working variables
my $binSize = $ENV{BIN_SIZE};
my ($chrom, $n5s, $sum5s, $prevBin) = ('', 0, 0);

# run the gc5base varStep lines
while (my $line = <STDIN>) {
    chomp $line;    
    
    # set the working chromosome
    if ($line =~ m/^v/){
        $line =~ m/chrom=(\S+)/;
        $chrom = $1;
        $prevBin = undef;
        next;
    }
    
    # get GC data over 5 genome bases
    my ($pos0, $val) = split("\t", $line);
    $pos0 = $pos0 - 1;      
    my $bin = int($pos0 / $binSize) * $binSize;

    # set the first bin on the chromosome
    if (!defined $prevBin) {
        $prevBin = $bin;
    }
    
    # print average mappability per bin
    if ($bin != $prevBin) {
        my $gc = $sum5s / $n5s / 100;
        print join("\t", $chrom, $prevBin, $prevBin + $binSize, ".",
                         int($gc * 1000 + 0.5) / 1000, "."), "\n";
        ($n5s, $sum5s, $prevBin) = (0, 0, $bin);
    }
    
    # aggregate over all 5-base chunks in a bin
    $n5s++;
    $sum5s += $val;

}

# print last bin
my $gc = $sum5s / $n5s / 100;
print join("\t", $chrom, $prevBin, $prevBin + $binSize, ".",
                 int($gc * 1000 + 0.5) / 1000, "."), "\n";
