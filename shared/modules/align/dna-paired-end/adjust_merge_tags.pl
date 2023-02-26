use strict;
use warnings;

# change the fastp merged tag to something short that is carried in QNAME
#   QNAME trailing :2 = merged by fastp, :0 = unmerged (later, :1 = merged with alignment guidance)
#   value of 2 allows fastp highest quality merges to sort first
#   if :0, 2 reads = unmerged pair, 1 read = orphaned read (where mate failed filters)
# unmerged pairs may still be merged later with alignment guidance

# support single-read alignment of paired-end reads
my $SUPPRESS_SMART_PAIRING = $ENV{SUPPRESS_SMART_PAIRING};
my $prevName = "";

while(my $line = <STDIN>){
    if($line =~ m/^\@/){
        my @f = split(" ", $line);

        # single-read alignment, QNAME = readPairN:umi1:umi2:mergeLevel:readN
        if($SUPPRESS_SMART_PAIRING){ 
            $line = join(":", $f[0], $f[1] ? 2 : 0, $f[0] eq $prevName ? 2 : 1)."\n"; 
            $prevName = $f[0];

        # paired alignment, QNAME = readPairN:umi1:umi2:mergeLevel    
        } else { 
            $line = join(":", $f[0], $f[1] ? 2 : 0)."\n";
        }  
    }
    print $line;
}
