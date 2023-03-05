use strict;
use warnings;

# action:
#   pull up to $N_READS QUAL values from the head of each FASTQ file
# expects:
#     source $MODULES_DIR/source/set_read_file_vars.sh (sets FASTQ_FILE1, FASTQ_FILE2, SRA_FILES)
#     input as either paired fastq.gz files or a set of .sra files

# initialize reporting
our $action  = "pull-quals";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);

# environment variables
fillEnvVar(\my $FASTQ_FILE1,        'FASTQ_FILE1');
fillEnvVar(\my $FASTQ_FILE2,        'FASTQ_FILE2');
fillEnvVar(\my $SRA_FILES,          'SRA_FILES');
fillEnvVar(\my $N_READS,            'N_READS');
fillEnvVar(\my $N_TERMINAL_BASES,   'N_TERMINAL_BASES');

# set the file input handles
if($FASTQ_FILE1){
    open my $inH1, "-|", "zcat $FASTQ_FILE1 2>/dev/null" or throwError("could not open $FASTQ_FILE1: $!");
    my $inH2;
    if($FASTQ_FILE2){
        open $inH2, "-|", "zcat $FASTQ_FILE2 2>/dev/null" or throwError("could not open $FASTQ_FILE2: $!");
    } else {
        $inH2 = $inH1;
    }
    runReadPairs($inH1, $inH2);
} else {
    foreach my $sraFile(split(" ", $SRA_FILES)){
        open my $inH, "-|", "fastq-dump --stdout --split-files $sraFile" or throwError("could not open $sraFile: $!");
        runReadPairs($inH, $inH);       
    }
}

# run the reqested number of paired reads
sub runReadPairs {
    my ($inH1, $inH2) = @_;
    my $nReads = 0;
    print join("\t", "length", "read", "first", "last"), "\n";
    while (my $read1 = printRead($inH1)){
        my $read2 = printRead($inH2) or return;
        $nReads += 2;
        $nReads >= $N_READS and return;       
    }
}

# parse a FASTQ set of 4 lines
sub printRead {
    my ($inH) = @_;

    # get four fastq lines
    my $name    = <$inH>;
    $name or return; # EOF
    my $seq     = <$inH>;
    my $discard = <$inH>;
    my $qual    = <$inH>;

    # print the read length and three different average QUALs
    chomp $qual;
    print join("\t",
        length($qual),
        getAvgQual($qual),
        getAvgQual(substr($qual, 0, $N_TERMINAL_BASES)),
        getAvgQual(substr($qual, -$N_TERMINAL_BASES))
    ), "\n";
}
