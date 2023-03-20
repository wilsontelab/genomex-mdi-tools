use strict;
use warnings;
use File::Basename;

# write an index to allow rapid retrieval of long read sequences 
# from multiple source fastq files based on QNAME

my $indexFile = "$ENV{DATA_FILE_PREFIX}.align.fastq_index.txt.gz";
open my $indexH, "|-", "gzip -c | slurp -s 100M -o $indexFile" or die "could not open: $indexFile\n";
print STDERR "aligning:\n";

while (my $fastqFile = glob($ENV{FASTQ_FILES})) {

    my $fastqName = basename($fastqFile);
    print STDERR "  $fastqName\n";    
    $fastqName =~ s/.fastq.gz//;

    open my $inH, "-|", "slurp -s 100M gunzip -c $fastqFile" or die "could not open: $fastqFile\n";
    my $offset = 0;
    while(my $line1 = <$inH>){
        my $line2 = <$inH>;
        my $line3 = <$inH>;
        my $line4 = <$inH>;
        my $qName = (split(/\s+/, $line1))[0];
        my $read = $line1.$line2.$line3.$line4;
        print $indexH join("\t", substr($qName, 1), $fastqName, $offset), "\n";
        print $read;
        $offset += length($read);
    }
    close $inH;
}

close $indexH;

1;
