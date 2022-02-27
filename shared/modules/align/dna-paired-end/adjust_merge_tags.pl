use strict;
use warnings;

# change the fastp merged tag to something short that is carried in QNAME
#   QNAME trailing :2 = merged by fastp, :0 = unmerged (later, :1 = merged with alignment guidance)
#   if :0, 2 reads = unmerged pair, 1 read = orphaned read (where mate failed filters)
# unmerged pairs may still be merged later with alignment guidance

while(my $line = <STDIN>){
    if($line =~ m/^\@/){
        my @f = split(" ", $line);
        $line = join(":", $f[0], $f[1] ? 2 : 0)."\n"; # 2 allow fastp highest quality merges to sort first
    }
    print $line;
}
