use strict;
use warnings;

#----------------------------------------------------------
# genome manipulations
#----------------------------------------------------------

# working variables
our (@canonicalChroms, %chromIndex, %revChromIndex);

# parsers for restricting work to properly ordered canonical chromosomes
sub setCanonicalChroms { 
    my %canonicalChroms = map { $_ => 1 } split(/\s+/, $ENV{GENOME_CHROMS});
    sub getPushValue{
        my ($chr, $canonicalChroms) = @_;
        my $chrom = "chr$chr";
        $$canonicalChroms{$chrom} ? $chrom : ();
    } 
    @canonicalChroms = map { getPushValue($_, \%canonicalChroms) } (1..90, 'X');
    $ENV{SUPPRESS_CHR_Y} or push @canonicalChroms, getPushValue('Y', \%canonicalChroms); # chrY included unless specifically excluded
    $ENV{USE_CHR_M}     and push @canonicalChroms, getPushValue('M', \%canonicalChroms); # chrM excluded unless specifically included
    %chromIndex    = map { $canonicalChroms[$_] => $_ + 1 } 0..$#canonicalChroms; # 1-referenced chrom indices, i.e., chr3 => 3
    %revChromIndex = map { $_ + 1 => $canonicalChroms[$_] } 0..$#canonicalChroms;
    $chromIndex{'*'} = 99; # special handling of unmapped reads
    $revChromIndex{99} = '*';
}

1;
