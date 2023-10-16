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
    my ($qry, $cigar, $noClip) = @_;
    my @qry = split("", $qry);
    my @qryOnRef;
    my $qryI = 0;
    my $nDeleted = 0;
    my $nInserted = 0;
    my @insI;
    if($noClip){
        $cigar =~ s/^(\d+)[S|H]//g and @qry = @qry[$1..$#qry];
        $cigar =~ s/(\d+)[S|H]$//g and @qry = @qry[0..($#qry - $1)];
    }
    while ($cigar =~ (m/(\d+)(\w)/g)) { 
        my ($size, $operation) = ($1, $2);
        if($operation eq 'D'){
            push @qryOnRef, (("-") x $size);            
            $nDeleted += $size;
        } elsif($operation eq 'I'){
            push @insI, $qryI - 1 + $nDeleted - $nInserted;
            $nInserted += $size;
            $qryI += $size; 
        } else {
            push @qryOnRef, @qry[$qryI..($qryI + $size - 1)];
            $qryI += $size;
        } 
    }
    foreach my $i(@insI){ $qryOnRef[$i] = "+" } # mark the position to the left of each novel insertion
    return \@qryOnRef;
}

sub getN50 {
    my (@sizes) = @_;
    my $sumSizes = 0;
    map { $sumSizes += $_ } @sizes;
    my $halfSumSizes = $sumSizes / 2;
    my $runningSum = 0;
    foreach my $size(sort {$b <=> $a} @sizes){
        $runningSum += $size;
        $runningSum >= $halfSumSizes and return $size;
    }
}

1;
