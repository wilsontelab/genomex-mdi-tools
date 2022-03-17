use strict;
use warnings;

#----------------------------------------------------------
# manipulations related to capture or amplicon target regions
#----------------------------------------------------------

# output variables
our ($nRegions, $sumTargetLens, $sumPaddedTargetLens) = (0, 0, 0);
our %targetRegions;

# constants
use constant {
    ON_TARGET   => 'T', # T for 'target'
    NEAR_TARGET => 'A', # A for 'adjacent'
    OFF_TARGET  => '-',
};

# operating parameters
use vars qw($TARGETS_BED $REGION_PADDING $TARGET_SCALAR);
defined $REGION_PADDING or $REGION_PADDING = 0;
defined $TARGET_SCALAR  or $TARGET_SCALAR  = 1;

# load the target regions
sub loadTargetRegions {
    my ($quiet) = @_;
    ($TARGETS_BED and $TARGETS_BED ne "null") or return;

    # first pass to record target regions (T=target type)
    $nRegions = loadTargetRegions_(\$sumTargetLens, 0, ON_TARGET);

    # second pass to record padded regions (A=adjacent type)
    $REGION_PADDING and loadTargetRegions_(\$sumPaddedTargetLens, $REGION_PADDING, NEAR_TARGET);

    # report target summary to log
    if(!$quiet){
        printCount($nRegions,            'nRegions',            'target regions');
        printCount($sumTargetLens,       'sumTargetLens',       'total bp covered by target regions');
        $REGION_PADDING and 
        printCount($sumPaddedTargetLens, 'sumPaddedTargetLens', 'total bp covered by padded target regions');
    }
}
sub loadTargetRegions_ {
    my ($sumLens, $regPad, $type) = @_;
    open my $inH, "<", $TARGETS_BED or die "could not open $TARGETS_BED: $!\n";
    my $regId = 0;
    while(my $line = <$inH>){
        $regId++;
        chomp $line;
        $line =~ s/\r//g;        
        my ($chr, $start, $end) = split("\t", $line);    
        for(my $pos =  int(($start - $regPad) / $TARGET_SCALAR);
               $pos <= int(($end   + $regPad) / $TARGET_SCALAR);
               $pos++){    
            $targetRegions{$chr} and $targetRegions{$chr}{$pos} and next;    
            $targetRegions{$chr}{$pos} = [$regId, $type]; # lookup for the target region associated with a coordinate
        }
        $$sumLens += ($end + $regPad) - ($start - $regPad);
    }
    close $inH;
    return $regId;
}

# codified form of the relationship of two coordinate pairs with respect to target regions
sub getTargetClass { 
    my ($chr1, $chr2, $pos1, $pos2) = @_;
    my $ct1 = $targetRegions{$chr1} ? ($targetRegions{$chr1}{int($pos1 / $TARGET_SCALAR)} || 0) : 0;
    my $ct2 = $targetRegions{$chr2} ? ($targetRegions{$chr2}{int($pos2 / $TARGET_SCALAR)} || 0) : 0;
    if($ct1 and $ct2){
        my ($target1, $type1) = @$ct1;
        my ($target2, $type2) = @$ct2;    
        my $type = join("", reverse sort {$a cmp $b} $type1, $type2); # yields TT, TA or AA
        return $target1 == $target2 ?
            $type :     # both outer molecule endpoints in same target region,  uppercase
            lc($type);  # outer molecule endpoints in different target regions, lowercase
    } elsif($ct1 or $ct2){
        my ($target, $type) = $ct1 ? @$ct1 : @$ct2;
        return lc($type).OFF_TARGET; # one endpoint is/near in a target region   
    } else {
        return OFF_TARGET.OFF_TARGET; # neither endpoint is in a target region
    }   
}

1;
