use strict;
use warnings;

# action:
#     prepare interleaved FASTQ from different input types, including SRA files
#     extract UMI sequences and skip base(s) in preparation for read alignment to genome
#     if UMIs not in use, use dummy UMI index of 1 for all reads
# expects:
#     source $MODULES_DIR/align/set_read_file_vars.sh (sets FASTQ_FILE1, FASTQ_FILE2, SRA_FILES)
#     input as either paired fastq.gz files or a set of .sra files

# initialize reporting
our $action  = "prepare_reads";
my ($nInputPairs, $nOutputPairs, $nInputReads,
    $nKnownUmis, $nInferredUmis, $nFailedUmis) = (0) x 10;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
require "$perlUtilDir/workflow.pl";
resetCountFile();

# environment variables
fillEnvVar(\my $FASTQ_FILE1, 'FASTQ_FILE1');
fillEnvVar(\my $FASTQ_FILE2, 'FASTQ_FILE2');
fillEnvVar(\my $SRA_FILES,   'SRA_FILES');
fillEnvVar(\my $UMI_FILE,    'UMI_FILE', 1, "");

# constants
use constant {
    QNAME => 0,
    SEQ   => 1,
    QUAL  => 2,
    UMI   => 3
};

# load the UMI library
my (%knownUmis, %umis, $skipLen, $seqRegex);
my ($umiLen, $nUmi, $nObsUmi) = (0) x 10;
my $isFixedUmi = ($UMI_FILE ne "");
if($isFixedUmi){
    open my $inH, "<", $UMI_FILE or throwError("could not open $UMI_FILE: $!");
    my $discardHeader = <$inH>;
    while (my $line = <$inH>){
        $nUmi++;    
        chomp $line;
        $line =~ s/\r//g;
        my @line = split("\t", $line);
        my $umi = $line[0]; # the known, i.e., expected UMI value, must be in first column
        $knownUmis{$umi}++;
        my @umi = split("", $umi);
        $umiLen = @umi;
        foreach my $i(0..$#umi){ # allow one base mismatch from known UMIs
            my @obs = @umi;
            foreach my $base(qw(A C G T)){
                $obs[$i] = $base;
                $umis{join("", @obs)} = $nUmi; # key=UMI as observed, value=inferred source UMI as 1-referenced index
                $nObsUmi++;
            }
        }
    }
    close $inH;

    # report fixed UMI metadata
    printCount($umiLen,  'umiLen',  'UMI length');
    printCount($nUmi,    'nUmi',    'expected UMI values');
    printCount($nObsUmi, 'nObsUmi', 'allowed UMI values (up to 1 base mismatch from expected)');  

    # declare the read pattern
    $skipLen = 1; # the single-base A addition used in library prep/1-base primer tail
    $seqRegex = qr/^(.{$umiLen}).{$skipLen}(.+)/;
}

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
            my $name = join(":", $$read1[QNAME], $$read1[UMI], $$read2[UMI]); # append paired UMIs to read pair id
            print join("\n", $name, $$read1[SEQ], '+', $$read1[QUAL]), "\n";  # print interleaved read pairs
            print join("\n", $name, $$read2[SEQ], '+', $$read2[QUAL]), "\n";        
        }
    }
}

# parse a FASTQ set of 4 lines
sub getRead {
    my ($inH) = @_;

    # get four fastq lines
    my $name    = <$inH>;
    $name or return;    
    my $seq     = <$inH>;
    my $discard = <$inH>;
    my $qual    = <$inH>;

    # process read pair name
    chomp $name;
    my @name = split(/\s/, $name, 2); 
    
    # process seq and qual  
    my $umi;    
    if($isFixedUmi){
        $nInputReads++;
        $seq =~ m/$seqRegex/;
        $umi = $umis{$1} || 0; # revert observed UMIs to the inferred parent UMI value
        $seq = $2;
        if($knownUmis{$1}){
            $nKnownUmis++;
        } elsif($umi == 0){
             $nFailedUmis++;
        } else {
            $nInferredUmis++;
        }       
        $qual =~ m/$seqRegex/;
        $qual = $2;        
    } else {
        chomp $seq;
        chomp $qual;        
        $umi = 1; # use dummy UMI values when UMIs not present; seq and qual passed as is
    }
    return [$name[0], $seq, $qual, $umi];
}

# print summary information
printCount($nInputPairs,   'nInputPairs',   'total input read pairs');
if($isFixedUmi){
    printCount($nOutputPairs,  'nOutputPairs',  'output read pairs');
    printCount($nInputReads,   'nInputReads',   'input reads');
    printCount($nKnownUmis,    'nKnownUmis',    'reads matched an expected UMI exactly');
    printCount($nInferredUmis, 'nInferredUmis', 'reads had a 1 base mismatch from an expected UMI');
    printCount($nFailedUmis,   'nFailedUmis',   'reads failed UMI lookup');
}
