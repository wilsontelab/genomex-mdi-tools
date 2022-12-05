use strict;
use warnings;

# calculate the weighted average run score over all bases in each bin

# working variables
my $chrom   = $ENV{CHROMOSOME};
my $binSize = $ENV{BIN_SIZE};
my ($binStart, $binEnd,      $weightedSum, $sumWeights) =
   (0,         $binSize - 1, 0,            0); # all working coordinates 0-based
my $binN = 1; # incremented over all bins

# run the input
while (my $line = <STDIN>) {
    chomp $line;    
    my ($runStart, $runEnd, $score) = split("\t", $line);
    $runEnd--; # all working coordinates 0-based
    
    # finish previous bin and add empty full-gap bins
    while($runStart > $binEnd) { printBin() }

    # handle runs that start in the working bin
    HANDLE_RUN: if ($runEnd <= $binEnd) { # bin entirely contains the run
        my $weight = $runEnd - $runStart + 1;
        $weightedSum += ($score * $weight);
        $sumWeights += $weight;
    } else { # run overruns the bin
        my $weight = $binEnd - $runStart + 1;
        $weightedSum += ($score * $weight);
        $sumWeights += $weight;
        printBin();
        $runStart = $binStart;
        goto HANDLE_RUN;
    }
}

# finish the last bin on the chromosome
$sumWeights and printBin();

# commit a bin
sub printBin {
    my $score = $sumWeights ? $weightedSum / $sumWeights : 0;
    print join("\t", $chrom, $binStart, $binEnd + 1, $binN, $score, "."), "\n";
    $binN++;
    $binStart += $binSize;
    $binEnd   += $binSize;
    ($weightedSum, $sumWeights) = (0, 0);
}
