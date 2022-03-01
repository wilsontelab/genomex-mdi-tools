use strict;
use warnings;

# utility functions common to svCapture scripts

# constants
use constant {
    ON_TARGET => 'T',
    NEAR_TARGET => 'A', # A for adjacent
    OFF_TARGET => '-',
};

# working variables
use vars qw(%ACGTNMatches %noIndelMatches @maxMergeMismatch);

#----------------------------------------------------------
# sequence manipulations
#----------------------------------------------------------
sub noIndelMatch { # compare two equal length sequences without allowing indels
    my ($seq1, $seq2) = @_; # expects only ACGTN in input, and returns all IUPAC
    my $length = length($seq1); # caller must ensure that seq lengths are the same
    my ($seq, $nMismatch, $maxNMismatch)
     = ('',   0,          $maxMergeMismatch[$length]);
    foreach my $i(0..$length-1){  
        my $pairedBases = substr($seq1, $i, 1).substr($seq2, $i, 1);
        my $base = $ACGTNMatches{$pairedBases};
        if($base){          
            $seq .= $base;
        } else {
            $nMismatch++; # abort early when can't succeed, estimate overall match count based on density
            $nMismatch > $maxNMismatch and return (int($length * ($i + 1 - $nMismatch) / ($i + 1)));
            $seq .= $noIndelMatches{$pairedBases} || 'N';
        } 
    }
    return ($length - $nMismatch, $seq);
}

#----------------------------------------------------------
# load sample-specific data resources
#----------------------------------------------------------
# the capture target regions
our ($nRegions, $targetScalar, $sumTargetLens, $sumPaddedTargetLens)  = (0, 10, 0, 0); # use 10 bp resolution
our %captureTargets;
my $regionPadding = $ENV{REGION_PADDING} || 0;
#sub loadTargetRegions { 
#    if($ENV{CAPTURE_BED}){
#        open my $inH, "<", $ENV{CAPTURE_BED} or die "could not open $ENV{CAPTURE_BED}: $!\n";
#        while(my $line = <$inH>){
#            $nRegions++;
#            chomp $line;
#            $line =~ s/\r//g;
#            my ($chr, $start, $end) = split("\t", $line);
#            for(my $pos=int(($start-$regionPadding)/$targetScalar);
#                  $pos<=int(($end  +$regionPadding)/$targetScalar);
#                  $pos++){
#                $captureTargets{$chr} and $captureTargets{$chr}{$pos} and next;
#                $captureTargets{$chr}{$pos} = $nRegions; # lookup for the target region associated with a coordinate
#            }
#            $sumTargetLens += $end - $start;
#        }
#        close $inH;    
#    }
#    printCount($nRegions, 'nRegions', 'capture target regions');
#    printCount($sumTargetLens, 'sumTargetLens', 'total bp covered by target regions');    
#}
sub loadTargetRegions {
    my ($quiet) = @_;
    if($ENV{CAPTURE_BED}){
        # first pass to record target regions (T=target type)
        $nRegions = loadTargetRegions_(\$sumTargetLens, 0, ON_TARGET);
        # second pass to record padded regions (A=adjacent type)
        $regionPadding and loadTargetRegions_(\$sumPaddedTargetLens, $regionPadding, NEAR_TARGET);
    }
    if (!$quiet) {
        printCount($nRegions, 'nRegions', 'capture target regions');
        printCount($sumTargetLens, 'sumTargetLens', 'total bp covered by target regions');
        printCount($sumPaddedTargetLens, 'sumPaddedTargetLens', 'total bp covered by padded target regions');
    }
}
sub loadTargetRegions_ {
    my ($sumLens, $regPad, $type) = @_;
    open my $inH, "<", $ENV{CAPTURE_BED} or die "could not open $ENV{CAPTURE_BED}: $!\n";
    my $regId = 0;
    while(my $line = <$inH>){
        $regId++;
        chomp $line;
        $line =~ s/\r//g;        
        my ($chr, $start, $end) = split("\t", $line);    
        for(my $pos =  int(($start - $regPad) / $targetScalar);
               $pos <= int(($end   + $regPad) / $targetScalar);
               $pos++){    
            $captureTargets{$chr} and $captureTargets{$chr}{$pos} and next;    
            $captureTargets{$chr}{$pos} = [$regId, $type]; # lookup for the target region associated with a coordinate
        }
        $$sumLens += ($end + $regPad) - ($start - $regPad);
    }
    close $inH;
    return $regId;
}
sub getTargetClass { # codified form of the relationship of two coordinate pairs with respect to capture targets
    my ($chr1, $chr2, $pos1, $pos2) = @_;
    my $ct1 = $captureTargets{$chr1} ? ($captureTargets{$chr1}{int($pos1/$targetScalar)} || 0) : 0;
    my $ct2 = $captureTargets{$chr2} ? ($captureTargets{$chr2}{int($pos2/$targetScalar)} || 0) : 0;
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

# fragment size distribution to aid in SV filtering and calling
our (%insertSizes, $maxTLen, $minSvSize);
sub loadInsertSizes { 
    my ($svSizeFactor, $percentile, $insSizeFile) = @_;
    $svSizeFactor or $svSizeFactor = 2;    
    $percentile or $percentile = 0.99; # maxTLen defaults to 99%th percentile of insert size distribution
    $insSizeFile or $insSizeFile = "$ENV{EXTRACT_PREFIX}.insertSizeDistribution.txt";
    open my $inH, "<", $insSizeFile or die "could not open $insSizeFile: $!\n";
    $maxTLen = undef; # for compare
    while(my $line = <$inH>){
        chomp $line;
        my ($insertSize, $freq, $cumFreq) = split("\t", $line);
        $insertSizes{$insertSize} = $freq;
        !$maxTLen and $cumFreq >= $percentile and $maxTLen = $insertSize;
        $maxTLen and !$minSvSize and $minSvSize = $maxTLen * $svSizeFactor;
    }
    close $inH;    
}

1;

