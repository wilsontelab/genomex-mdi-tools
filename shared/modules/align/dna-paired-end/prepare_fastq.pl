use strict;
use warnings;

# action:
#     prepare interleaved FASTQ from different input types, including SRA files
#     extract UMI sequences and skip base(s) in preparation for read-pair alignment to genome
#         if UMIs not in use, append dummy UMI index of 1 for all read in all pairs
#     apply read quality filtering as requested (any failed read fails both reads of the pair)
# expects:
#     source $MODULES_DIR/source/set_read_file_vars.sh (sets FASTQ_FILE1, FASTQ_FILE2, SRA_FILES)
#     input as either paired fastq.gz files or a set of .sra files

# initialize reporting
our $action  = "prepare_fastq";
my ($nInputPairs, $nOutputPairs, $nInputReads,
    $nUmisExpected, $nUmisInferred, $nUmisFailed, $nQualFailed) = (0) x 10;

# load dependencies
my $perlUtilDir = "$ENV{SHARED_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general umi);
resetCountFile();

# environment variables
fillEnvVar(\my $FASTQ_FILE1,        'FASTQ_FILE1');
fillEnvVar(\my $FASTQ_FILE2,        'FASTQ_FILE2');
fillEnvVar(\my $SRA_FILES,          'SRA_FILES');
fillEnvVar(\my $UMI_FILE,           'UMI_FILE',       1, "");
fillEnvVar(\my $UMI_SKIP_BASES,     'UMI_SKIP_BASES', 1, 1); # the single-base A addition used in library prep/1-base primer tail
fillEnvVar(\my $MIN_QUAL,           'MIN_QUAL', 1, "NA");
fillEnvVar(\my $N_TERMINAL_BASES,   'N_TERMINAL_BASES', 1, 30);
fillEnvVar(\my $MAX_HOMOPOLYMER,    'MAX_HOMOPOLYMER', 1, 0);

# constants
use constant {
    SEQ   => 0,
    QUAL  => 1,
    UMI   => 2
};

# load the UMI library
my $seqRegex = loadFixedUmis($UMI_FILE, $UMI_SKIP_BASES);
my $isFixedUmi = $seqRegex ? 1 : 0;
use vars qw(%expectedUmis %allowedUmis);

# prepare for quality filtering
my ($areQualFiltering, $MIN_QUAL_READ, $MIN_QUAL_FIRST, $MIN_QUAL_LAST) = (0, 0, 0, 0);
if($MIN_QUAL and $MIN_QUAL ne "NA"){
    ($MIN_QUAL_READ, $MIN_QUAL_FIRST, $MIN_QUAL_LAST) = split(":", $MIN_QUAL);
    $areQualFiltering = ($MIN_QUAL_READ || $MIN_QUAL_FIRST || $MIN_QUAL_LAST);
}
my $BAD_HOMOPOLYMER = $MAX_HOMOPOLYMER + 1;
my $homopolymer = qr/([A,N]{$BAD_HOMOPOLYMER}|[C,N]{$BAD_HOMOPOLYMER}|[G,N]{$BAD_HOMOPOLYMER}|[T,N]{$BAD_HOMOPOLYMER})/;

# set the file input handles
if($FASTQ_FILE1){
    open my $inH1, "-|", "zcat $FASTQ_FILE1" or throwError("could not open $FASTQ_FILE1: $!");
    open my $inH2, "-|", "zcat $FASTQ_FILE2" or throwError("could not open $FASTQ_FILE2: $!");
    runReadPairs($inH1, $inH2);
    close $inH1;
    close $inH2;    
} else {
    foreach my $sraFile(split(" ", $SRA_FILES)){
        open my $inH, "-|", "fastq-dump --stdout --split-files $sraFile" or throwError("could not open $sraFile: $!");
        runReadPairs($inH, $inH);
        close $inH;        
    }
}

# run the paired reads
# no need to multi-thread since the downstream alignment process is rate limiting
sub runReadPairs {
    my ($inH1, $inH2) = @_;
    while (my $read1 = getRead($inH1)){
        my $read2 = getRead($inH2);
        $nInputPairs++;
        if($$read1[UMI] and $$read2[UMI]){ # both UMIs must be known to proceed; discard all other read pairs
            $nOutputPairs++;
            my $name = "@".join(":", $nOutputPairs, $$read1[UMI], $$read2[UMI]); # append paired UMIs to numeric read pair id
            print join("\n", $name, $$read1[SEQ], '+', $$read1[QUAL]), "\n"; # print interleaved read pairs
            print join("\n", $name, $$read2[SEQ], '+', $$read2[QUAL]), "\n";        
        }
    }
}

# parse a FASTQ set of 4 lines
sub failReadQuality {
    $nQualFailed++;
    return [];
}
sub getRead {
    my ($inH) = @_;

    # get four fastq lines
    my $name    = <$inH>; # discarded in favor over numeric read pair id
    $name or return; # EOF
    my $seq     = <$inH>;
    my $discard = <$inH>;
    my $qual    = <$inH>;

    # # process read pair name
    # chomp $name;
    # my @name = split(/\s/, $name, 2); 
    
    # process seq and qual  
    my $umi;    
    if($isFixedUmi){
        $nInputReads++;
        $seq =~ m/$seqRegex/;
        $seq = $2;        
        $umi = $allowedUmis{$1} || 0; # revert observed UMIs to the inferred parent UMI value
        if($expectedUmis{$1}){
            $nUmisExpected++;
        } elsif($umi == 0){
            $nUmisFailed++;
        } else {
            $nUmisInferred++;
        }       
        $qual =~ m/$seqRegex/;
        $qual = $2;        
    } else {
        chomp $seq;
        chomp $qual;        
        $umi = 1; # use dummy UMI values when UMIs not present; seq and qual passed as is
    }

    # apply qual filtering upstream of alignment and merging
    # use with discretion, these are computationally more costly
    if($areQualFiltering){
        !$MIN_QUAL_READ  or getAvgQual($qual) >= $MIN_QUAL_READ  or return failReadQuality();
        !$MIN_QUAL_FIRST or getAvgQual(substr($qual, 0, $N_TERMINAL_BASES)) >= $MIN_QUAL_FIRST or return failReadQuality();
        !$MIN_QUAL_LAST  or getAvgQual(substr($qual,   -$N_TERMINAL_BASES)) >= $MIN_QUAL_LAST  or return failReadQuality(); # fast, usually sufficient
    }
    if($MAX_HOMOPOLYMER){
        $seq =~ m/$homopolymer/;
        $1 and return failReadQuality();
    }
    return [$seq, $qual, $umi]; # $name[0]
}

# print summary information
printCount($nInputReads,   'nInputReads',   'input reads');
printCount($nInputPairs,   'nInputPairs',   'input read pairs');
if($isFixedUmi){
    printCount($nUmisExpected, 'nUmisExpected', 'reads matched an expected UMI exactly');
    printCount($nUmisInferred, 'nUmisInferred', 'reads had a 1 base mismatch from an expected UMI');
    printCount($nUmisFailed,   'nUmisFailed',   'reads failed UMI lookup');    
}
if($areQualFiltering or $MAX_HOMOPOLYMER){
    printCount($nQualFailed,   'nQualFailed',   'reads failed one or more quality filters');
}
printCount($nOutputPairs,  'nOutputPairs',  'output read pairs');
