use strict;
use warnings;

# change the fastp merged tag to something short that is carried in QNAME
#   QNAME trailing :2 = merged by fastp, :0 = unmerged (later, :1 = merged with alignment guidance)
#   value of 2 allows fastp highest quality merges to sort first
#   if :0, 2 reads = unmerged pair, 1 read = orphaned read (where mate failed filters)
# unmerged pairs may still be merged later with alignment guidance
# if requested, append an additional trailing :int with the number of bases of read overlap as set by fastp

# support single-read alignment of paired-end reads
my $SUPPRESS_SMART_PAIRING = $ENV{SUPPRESS_SMART_PAIRING};
my $prevName = "";

# for legacy consistency, require callers to request processing of fastq overlap lengths
# QNAME merged_150_15 means that 150bp are from read1, and 15bp are from read2
my $APPEND_FASTP_OVERLAP = $ENV{APPEND_FASTP_OVERLAP};
my $READ_LEN_x2 = 2 * $ENV{READ_LEN};
sub getMergeOverlap {
    $_[0] or return 0;
    my @x = split("_", $_[0]);
    $READ_LEN_x2 - ($x[1] + $x[2]);
}

# run the interleaved pairs
my $lineN = 0;
while(my $line = <STDIN>){
    unless($lineN % 4){
        my @f = split(" ", $line); # split drops trailing whitespace when splitting on whitespace, i.e., chomp occurs implicitly
        
        # single-read alignment, QNAME = readPairN:umi1:umi2:mergeLevel[:overlap]:readN
        if($SUPPRESS_SMART_PAIRING){ 
            $line = join(":", 
                $f[0], 
                $f[1] ? 2 : 0, 
                $APPEND_FASTP_OVERLAP ? getMergeOverlap($f[1]) : (),
                $f[0] eq $prevName ? 2 : 1
            )."\n"; 
            $prevName = $f[0];

        # paired alignment, QNAME = readPairN:umi1:umi2:mergeLevel[:overlap] 
        } else { 
            $line = join(":", 
                $f[0], 
                $f[1] ? 2 : 0, 
                $APPEND_FASTP_OVERLAP ? getMergeOverlap($f[1]) : ()
            )."\n";
        }    
    }
    print $line;
    $lineN++;
}
