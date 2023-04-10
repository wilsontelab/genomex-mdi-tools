use strict;
use warnings;

# use a previously constructed fastq index to give random access to a read's SEQ and QUAL

# ba59c9b1-7fbd-43c3-b7b0-bc662bf9770a    PAK57679_pass_909c8bfb_7f82fa6d_0       0
# ccc98ece-eeab-4b5d-bed7-7e86552cc6ec    PAK57679_pass_909c8bfb_7f82fa6d_0       1137

# load the index
my (%index, $prevQName, $prevFile, $prevOffset);
my $indexFile = "$ENV{DATA_FILE_PREFIX}.align.fastq_index.txt.gz";
open my $indexH, "-|", "zcat $indexFile" or die "could not open: $indexFile\n";
while (my $line = <$indexH>){
    chomp $line;
    my ($qName, $file, $offset) = split("\t", $line);
    $index{$qName} = [$file, $offset];
    if($prevQName){
        $index{$prevQName}[2] = $file eq $prevFile ? $offset - $index{$prevQName}[1] : 1e7;
    }
    ($prevQName, $prevFile, $prevOffset) = ($qName, $file, $offset);
}
$index{$prevQName}[2] = 1e7;
close $indexH;

# retrieve a single read
my ($openFile, $readH) = ("");
sub getIndexedRead {
    my ($qName) = @_; 
    my $fileName = "$index{$qName}[0].fastq.gz";

    # my $tmpFile = "$ENV{TMP_DIR_WRK}/$fileName";
    my $tmpFile = "/scratch/wilsonte_root/wilsonte0/wilsonte/7754-SA_P2_SOLO/pilot/HCT_untargeted/DEBUG.fastq";

    # print "$ENV{INPUT_DIR}/$fileName\n";
    # print "$tmpFile\n";

    unless($openFile eq $tmpFile){ # can generally expect callers to process reads in order, so keep handles open
        $readH and close $readH;
        # -e $tmpFile or system("zcat $ENV{INPUT_DIR}/$fileName > $tmpFile");
        open $readH, "<", $tmpFile or die "could not open: $tmpFile\n";
        $openFile = $tmpFile;
    }

    my $buffer = "";
    seek $readH, $index{$qName}[1], 0;
    read $readH, $buffer, $index{$qName}[2];
    return split("\n", $buffer);
}

1;
