use strict;
use warnings;

# for parabricks gpu acceleration, split interleaved and merged reads to three temporary files
#   - unmerged_1.fastq
#   - unmerged_2.fastq
#   - merged.fastq (includes orphaned reads)

print STDERR "pre-filtering, trimming and merging reads\n";

open my $unmergedH_1, "|-", "pigz -c -p $ENV{N_CPU} > $ENV{UNMERGED_FILE_1}" or die "could not open: $ENV{UNMERGED_FILE_1}: \n";
open my $unmergedH_2, "|-", "pigz -c -p $ENV{N_CPU} > $ENV{UNMERGED_FILE_2}" or die "could not open: $ENV{UNMERGED_FILE_2}: \n";
open my $mergedH,     "|-", "pigz -c -p $ENV{N_CPU} > $ENV{MERGED_FILE}"     or die "could not open: $ENV{MERGED_FILE}: \n";

# run the interleaved pairs
my ($lineN, $readN, $prevQName, @buffer) = (0, 0);
while(my $line = <STDIN>){
    unless($lineN % 4){
        if($prevQName){
            if($prevQName ne $line){        
                if(@buffer == 1){
                    print $mergedH join("", @{$buffer[0]});
                } else {
                    print $unmergedH_1 join("", @{$buffer[0]});
                    print $unmergedH_2 join("", @{$buffer[1]});
                }
                ($readN, @buffer) = (0);
                $prevQName = $line;
            } else {
                $readN = 1;
            }
        } else {
            $prevQName = $line;  
        }
    }
    push @{$buffer[$readN]}, $line;
    $lineN++;
}

close $unmergedH_1;
close $unmergedH_2;
close $mergedH;

# @1:1:1:0
# AGGAGAATCACCAATCCCAGCAGTCTGAACTTGAGTTCTGAAAAACCTCGCTACCAAAGGCCAAAGTGCTCTGGGTCTCTAAGTAAACTTGAAAGGCAATCTAGGCCACAAGGACTGCAGCTCCGATGTGCTGATGTGAGTCTTAGGGCTG
# +
# FFFFFFF-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF8FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF-FFFFFFFFF-FF-FFFFFFFFFFFFF-FFFFFFFFFFFFFFFFFFFFFFFFFF8FF-FFFFF-FF
# @1:1:1:0
# TACTCAACTGAAGTGAGCAAGTCTCAGAGGTTCAGCCAAGGCCCTTGACATAGTACCTGTTTATTGCTGCTGATTGTTCAGGTCCCAAGGGCTCTTCAATTAGAAGGGGATAACTGCTGCCTGCCAGGGCTATGTCCTACTCTTCAAGGTG
# +
# FFF-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF-FFFFFFFFFFFFFF-FFF--F-FFFFFF-F-FFFFFFFFFFFFFFFFFF-FF8FF--FFFFFF8FFFFFFF--F-FFFFFFFFFFFFF-FFFFFFFFF8-F88FFF---FFF

# @2:1:1:0
# CTGACTACATTGAATAGTCTCCCAACTCACTAATAGAATCTGACCCTGTCTCTAAGATGCTTCTCTAACTCAGGCTGAAGAGCTAATGGTTAATAATAATGTGCTTCTTAAAGAGATATTTGTGTTGTAAAGGAGCACATTTGAAGGAGTG
# +
# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF-FFFFFFFFFF-FFF-FFFFFFFFFFFFFFFF-FFFFFFFFFFFFFFFFFFF-FFFFFFFFFFFFFFFFFFFF
# @2:1:1:0
# GTCCTTAGGCAAGTTAGTAGGAGCAAGGTTACATTGTAGTTGGCCTCAGTGGTCAGTCAGACTGGAGTTTAAATCTTTGCTTTACTAATGTACAACTCAGTACTGTGAGCATCAGTTTCCTCATCTATAAAAATAGCTACAATAGCTACCT
# +
# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF8F-FF8FFFFFF8FFFFFFFFFF-FFFFFFFFF-FFF
