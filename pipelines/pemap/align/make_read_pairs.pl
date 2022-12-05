use strict;
use warnings;
use List::Util qw(shuffle); 

# create a set of pseudo-randomized fragment endpoints that samples ~every genome base
# exactly twice, once as a fragment start, once as a fragment end
# output paired reads for each fragment as two sequential regions for samtools faidx to retrieve

# initialize reporting
our $script = "make_read_pairs";
our $error  = "$script error";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);
map { require "$perlUtilDir/sequence/$_.pl" } qw(faidx);

# environment variables
fillEnvVar(\our $CHROMOSOME,       'CHROMOSOME');
fillEnvVar(\our $GENOME_FASTA_SHM, 'GENOME_FASTA_SHM');
fillEnvVar(\our $READ_LENGTH,      'READ_LENGTH');

# initialize the genome
loadFaidx($GENOME_FASTA_SHM);

# initialize a set of distributions of genomic fragment insert sizes
my $minFragStart = 0;
my $maxFragStart = 1 * $READ_LENGTH - 1;
my $minFragEnd   = 2 * $READ_LENGTH - 2; # prevents read overruns in too-small molecules
my $maxFragEnd   = 3 * $READ_LENGTH - 3;
my @fragStarts   = $minFragStart..$maxFragStart; # 0-referenced
my @fragEnds     = $minFragEnd..$maxFragEnd;     
my @read1Starts  = @fragStarts;
my @read1Ends    = map { $_ + $READ_LENGTH - 1 } @read1Starts;
my $nDist = 1000;
my @dists = map {
    my @read2Ends = shuffle(@fragEnds);
    {
        read2Starts => [map { $_ - $READ_LENGTH + 1 } @read2Ends],
        read2Ends   => \@read2Ends        
    }
} 1..$nDist;

# print all matching read pair mimics on chromosome
my $chromLength = getChromLength($CHROMOSOME);
my @blockStarts; # 1-referenced
my $bs = 1;
while($bs < $chromLength - $maxFragEnd){
    push @blockStarts, $bs;
    $bs += $READ_LENGTH;
}
@blockStarts = shuffle(@blockStarts); # jump around the chromosome to help unsure mappability in bwa stream
foreach my $blockStart(@blockStarts){
    my $dist = $dists[int(rand($nDist))];
    foreach my $i(0..$#read1Starts){
        my $start = $blockStart + $read1Starts[$i];
        my $end   = $blockStart + $read1Ends[$i];
        print "$CHROMOSOME:$start-$end\n";
        $start = $blockStart + $$dist{read2Starts}[$i];
        $end   = $blockStart + $$dist{read2Ends}[$i];
        print "$CHROMOSOME:$start-$end\n";
    }
}
