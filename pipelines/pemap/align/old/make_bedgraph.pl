use strict;
use warnings;

# collapse endpoint hits to a bedgraph of coverage runs
# values at each base are either 0, 0.5, or 1, depending on whether
# it matched as a fragment start, end, both or neither

# initialize reporting
our $script = "make_bedgraph";
our $error  = "$script error";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$perlUtilDir/sequence/$_.pl" } qw(faidx);

# environment variables
fillEnvVar(\our $GENOME_FASTA_SHM, 'GENOME_FASTA_SHM');

# initialize the genome
use vars qw(@canonicalChroms %chromIndex %revChromIndex);
setCanonicalChroms();
loadFaidx($GENOME_FASTA_SHM);

# constants
use constant {
    CHROM => 0, # SAM fields
    POS => 1,
    COUNT => 2,
};

# working variables 
my ($prevChrom, $prevPos, $prevCount, $startPos) = ("", 0, -1, 1);

# run the base position
while(my $line = <STDIN>){
    chomp $line;
    my @line = split("\t", $line);
    if($prevChrom){
        if($prevChrom != $line[CHROM]){
            finishChromosome();
            startChromosome($line[CHROM], $line[POS] - 1);
        } elsif($prevPos != $line[POS] - 1){
            printRun($prevPos, $prevCount);
            printRun($line[POS] - 1, 0);
        } elsif($prevCount != $line[COUNT]){
            printRun($prevPos, $prevCount);
        }
    } else { # start the first chromosome
        startChromosome($line[CHROM], $line[POS] - 1);
    }
    ($prevChrom, $prevPos, $prevCount) = @line;
}
finishChromosome();

sub startChromosome {
    my ($chrom, $endPos) = @_;
    $prevChrom = $chrom;
    $endPos > 0 and printRun($endPos, 0);    
}
sub finishChromosome {
    printRun($prevPos, $prevCount);
    my $chromEndPos = getChromLength($revChromIndex{$prevChrom});
    $startPos <= $chromEndPos and printRun($chromEndPos, 0);
    ($prevChrom, $prevPos, $prevCount, $startPos) = ("", 0, -1, 1);
}
sub printRun {
    my ($endPos, $count) = @_;
    print join("\t", $revChromIndex{$prevChrom}, $startPos - 1, $endPos, $count / 2), "\n";
    $startPos = $endPos + 1;
}
