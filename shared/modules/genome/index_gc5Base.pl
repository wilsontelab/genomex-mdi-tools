use strict;
use warnings;

# convert UCSC gc5Base.wigVarStep.gz to bed and index for
# rapid retrieval of GC base content in an arbitrary genome region

#variableStep chrom=chr1 span=5
#10001   40
#10006   40
#10011   40
#10016   60
#10021   60
#10026   60
#10031   40
#10036   40
#10041   40

# check for something to do
my $gc5BaseWigFile = "$ENV{GENOME_UCSC_DIR}/$ENV{GENOME}.gc5Base.wigVarStep.gz";
my $gc5BaseBgzFile = "$ENV{GENOME_UCSC_DIR}/$ENV{GENOME}.gc5Base.txt.bgz";
-f $gc5BaseBgzFile and exit;
print STDERR "\nindexing $ENV{GENOME} gc5Base\n";

# working variables
my ($chrom);
open my $inH,  "-|", "zcat $gc5BaseWigFile" or die "could not open $gc5BaseWigFile: $!";
open my $outH, "|-", "bgzip -c > $gc5BaseBgzFile" or die "could not open $gc5BaseBgzFile: $!";

# run the gc5base varStep lines
while (my $line = <$inH>) {   
    if ($line =~ m/^v/){
        $line =~ m/chrom=(\S+)/;
        $chrom = $1;
        next;
    }
    print $outH join("\t", $chrom, $line);
}
close $inH;
close $outH;

# create the tabix index
system("tabix -s 1 -b 2 -e 2 $gc5BaseBgzFile");
print STDERR "done\n\n";
