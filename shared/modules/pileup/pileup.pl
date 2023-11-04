use strict;
use warnings;

# create a base-level map of all bam alignments, like samtools mpileup but more focused in purpose
# output is runs of unique base coverage as chrom,start,end,bases(M:nRefBases,A:nAltBases_A,etc.)

# initialize reporting
our $script = "pileup";
our $error  = "$script error";

# load dependencies
my $perlUtilDir = "$ENV{SHARED_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general faidx);

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
    ALN_N => 10,
    RNAME_INDEX => 11,
};

# environment variables
fillEnvVar(\our $N_CPU,                 'N_CPU'); 
fillEnvVar(\our $GENOME_FASTA,          'GENOME_FASTA');
fillEnvVar(\our $COORDINATE_BAM_FILE,   'COORDINATE_BAM_FILE');
fillEnvVar(\our $SHM_DIR_WRK,           'SHM_DIR_WRK');
fillEnvVar(\our $MIN_MAPQ,              'MIN_MAPQ');

# initialize the genome
use vars qw(@canonicalChroms %chromIndex);
setCanonicalChroms();
loadFaidx($GENOME_FASTA);

# process data by chromosome over multiple parallel threads
launchChildThreads(\&parseChromosome);
use vars qw(@readH @writeH);
my $writeH = $writeH[1];
foreach my $i(0..$#canonicalChroms){
    $writeH = $writeH[$i % $N_CPU + 1];   
    print $writeH "$canonicalChroms[$i]\n"; # commit to worker thread
}
finishChildThreads();

# child process to parse bam read pairs
sub parseChromosome {
    my ($childN) = @_;
    my $readH = $readH[$childN];    
    open our $faH, "<", $GENOME_FASTA or die "$error: could not find genome: $!\n";
    while(my $chrom = <$readH>){

        # initialize the chromosome
        chomp $chrom;
        our $shortChrom = $chrom;
        $shortChrom =~ s/chr//;
        my $chromSeq = getChromSeq($chrom);
        open our $outH, "|-", "gzip -c > $SHM_DIR_WRK/pileup.$chrom.txt.gz" 
            or die "$error: $chrom output stream failed: $!";
        open my $inH, "-|", "slurp -s 50M samtools view -F 1796 -q $MIN_MAPQ $COORDINATE_BAM_FILE $chrom" 
            or die "$error: $chrom input stream failed: $!";
        our ($pileupStartPos, $lastProcessedI, @pileup) = (1, -1);

        # thread alignment through a transient buffer for memory efficiency
        sub finishBases {
            my ($lastPos) = @_;
            my $lastI = $lastPos - $pileupStartPos;  
            my $hasBuffer = @pileup; # false on first process request on chromosome 
            my $workingRun = $hasBuffer ? $pileup[$lastProcessedI] : ""; # the one run always present in the pileup buffer, that might be extended by the next processed position
            my $hasBreak;

            # process the collection of nominated base values into a final, streamlined pileup at each base
            foreach my $i(($lastProcessedI + 1)..$lastI){
                if($pileup[$i]){
                    foreach my $base(keys %{$pileup[$i]}){ # suppress singleton variants as likely errors
                        if($base ne "M" and $pileup[$i]{$base} == 1){
                            $pileup[$i]{"M"} += 1;
                            $pileup[$i]{$base} = 0;
                        }
                    }
                }
                $pileup[$i] = $pileup[$i] ? 
                    join("", map { $pileup[$i]{$_} ? "$pileup[$i]{$_}$_" : () } keys %{$pileup[$i]}) : 
                    "0M";
                $pileup[$i] ne $workingRun and $hasBreak = 1; # if never hits, we are just extended the current run
            }
            $lastProcessedI = $lastI;
            ($hasBuffer and $hasBreak) or return;

            # further collapse runs of identical base coverage (e.g., M12) into single output lines
            TRY_NEXT_RUN: my ($runBases, $runEndI) = ($pileup[0], 0);
            while($runEndI < $lastProcessedI and $pileup[$runEndI + 1] eq $runBases){
                $runEndI++;
            }    
            if($runEndI < $lastProcessedI){ # thus, always keep one run in the pileup buffer
                my $nPos = $runEndI + 1;            
                print $outH join("\t", $shortChrom, $pileupStartPos, $nPos, $runBases), "\n";
                splice(@pileup, 0, $nPos);
                $pileupStartPos += $nPos;
                $lastProcessedI -= $nPos;
                goto TRY_NEXT_RUN;
            }
        }

        # run the alignments
        while(my $line = <$inH>){
            my @aln = split("\t", $line, 11);
            $aln[POS] > $pileupStartPos and finishBases($aln[POS] - 1);
            my $qryOnRef = getQryOnRef($aln[SEQ], $aln[CIGAR], 1); 
            foreach my $qI(0..$#$qryOnRef){
                my $base = ($$qryOnRef[$qI] eq substr($chromSeq, $aln[POS] - 1 + $qI, 1) or $$qryOnRef[$qI] eq "N") ? "M" : $$qryOnRef[$qI];
                $pileup[$aln[POS] - $pileupStartPos + $qI]{$base}++;        
            }
        }
        finishBases($pileupStartPos + @pileup - 1);
        print $outH join("\t", $shortChrom, $pileupStartPos, scalar(@pileup), $pileup[0]), "\n";
        close $inH;
        close $outH;

        print STDERR "  $chrom done\n";
    }
    close $faH;
}
