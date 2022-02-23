use strict;
use warnings;
use File::Copy;

my @EXPERIMENT_IDS   = split(" ", $ENV{EXPERIMENT_IDS});
my @EXPERIMENT_NAMES = split(" ", $ENV{EXPERIMENT_NAMES});
@EXPERIMENT_IDS == @EXPERIMENT_NAMES or die "different number of experiment-ids and experiment-names\n";

foreach my $i(0..$#EXPERIMENT_IDS){
    my $source = "$ENV{TASK_DIR}/$EXPERIMENT_IDS[$i]";
    my $dest   = "$ENV{TASK_DIR}/$EXPERIMENT_NAMES[$i]";
    move($source, $dest);    
}
