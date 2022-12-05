use strict;
use warnings;

# convert read alignments to a vector of mapped endpoint counts
# output as bedgraph-like, but without the chrom column

# initialize reporting
our $script = "parse_bwa";
our $error  = "$script error";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general faidx);

# environment variables
fillEnvVar(\our $CHROMOSOME,       'CHROMOSOME');
fillEnvVar(\our $MIN_MAPQ,         'MIN_MAPQ');
fillEnvVar(\our $GENOME_FASTA_SHM, 'GENOME_FASTA_SHM');

# initialize the genome
loadFaidx($GENOME_FASTA_SHM);
my $chromLength = getChromLength($CHROMOSOME);

# constants
use constant {
    QNAME => 0, # SAM fields
    FLAG => 1,
    RNAME => 2,
    POS => 3,
    MAPQ => 4,
    CIGAR => 5,
    RNEXT => 6,
    PNEXT => 7,
    TLEN => 8,
    SEQ => 9,
    #-------------
    _IS_PAIRED => 1, # SAM FLAG bits
    _PROPER => 2,
    _UNMAPPED => 4,
    _REVERSE => 16,
    _FIRST_IN_PAIR => 64,
    _SECOND_IN_PAIR => 128,
    _SUPPLEMENTARY => 2048,  
};

# working variables
my @counts = (0) x $chromLength;
my $rejectFlags = _UNMAPPED + _REVERSE + _SUPPLEMENTARY; # we know all read pairs are FF 

# run the alignments
# NB: do NOT enforce proper fragment assessments to allow one half of a read pair
# to yield a high quality alignment, even if the other cannot be mapped 
while(my $aln = <STDIN>){
    $aln =~ m/^\@/ and next;
    my @aln = (split("\t", $aln, 7));
    ($aln[FLAG] & $rejectFlags) and next;    
    $aln[RNAME] eq $CHROMOSOME or next;
    $aln[MAPQ] >= $MIN_MAPQ or next;
    my $pos = ($aln[FLAG] & _FIRST_IN_PAIR) ? $aln[POS]: getEnd($aln[POS], $aln[CIGAR]);
    $counts[$pos - 1]++;
}

# print the base-level map
my ($prevVal, $startPos) = (-1, 1);
foreach my $pos(1..$chromLength){
    my $val = min(2, $counts[$pos - 1]) / 2;
    if($prevVal >= 0 and $prevVal != $val){
        printRun($pos - 1);
    }
    $prevVal = $val;
}
printRun($chromLength);
sub printRun {
    my ($endPos) = @_;
    print join("\t", $startPos - 1, $endPos, $prevVal), "\n";
    $startPos = $endPos + 1;
}
