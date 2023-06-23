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

# get the average per-base Phred QUAL for a single read or segment of a read
sub getAvgQual {
    ($_[0] and length($_[0])) or return 0;
    my $sum = 0;
    map{ $sum += ord($_) } split("", $_[0]);
    $sum / length($_[0]) - 33;
}

# parse a CIGAR string to match a query SEQ to its reference
sub getQryOnRef {
    my ($qry, $cigar) = @_;
    my @qry = split("", $qry);
    my @qryOnRef;
    my $index = 0;
    my $nDeleted = 0;
    my @insPos;
    while ($cigar =~ (m/(\d+)(\w)/g)) { 
        my ($size, $operation) = ($1, $2);
        if($operation eq 'D'){
            $nDeleted += $size;
            push @qryOnRef, (("-") x $size);
        } elsif($operation eq 'I'){
            push @insPos, $index - 1 + $nDeleted;
            $index += $size; 
        } else {
            push @qryOnRef, @qry[$index..($index + $size - 1)];
            $index += $size;
        } 
    }
    foreach my $i(@insPos){ $qryOnRef[$i] = "+" } # mark the position to the left of each novel insertion
    return \@qryOnRef
}

1;
