use strict;
use warnings;

#----------------------------------------------------------
# general sequence manipulations
#----------------------------------------------------------

# regular expressions to process read clips
our $leftClip_  = qr/^(\d+)S/;
our $rightClip_ = qr/(\d+)S$/;

# reverse complement reads in place (i.e., by reference)
sub rc { 
    my ($seq) = @_;
    $$seq = reverse $$seq;
    $$seq =~ tr/ATGC/TACG/;
}

# get rightmost mapped read position in reference genome from POS and CIGAR
sub getEnd { 
    my ($start, $cigar) = @_;
    $cigar =~ s/\d+S//g;
    my $end = $start - 1;
    while ($cigar =~ (m/(\d+)(\w)/g)) {
        $2 eq "I" or $end += $1;
    }
    return $end;
}

1;
