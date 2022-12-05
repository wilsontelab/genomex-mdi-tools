use strict;
use warnings;

# adjust the names on samtools faidx output to create read pairs for bwa
# output FASTA 

# initialize reporting
our $script = "parse_faidx";
our $error  = "$script error";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);

# environment variables
fillEnvVar(\our $READ_LENGTH, 'READ_LENGTH');
my $readLenNewline = $READ_LENGTH + 1;

# run the read pairs in fasta format
my $i = 1;
while(my $name1 = <STDIN>){
    my $seq1 = uc(<STDIN>);
    my $name2 = <STDIN>;
    my $seq2 = uc(<STDIN>);
    length($seq1) != $readLenNewline and die "$error: wrong seq length: $seq1"; # check for proper read parsing
    length($seq2) != $readLenNewline and die "$error: wrong seq length: $seq2";
    if($seq1 =~ m/[^N]/ or $seq2 =~ m/[^N]/){ # skip pairs with all N bases, can never align

        my $qual1 = $seq1;
        $qual1 =~ tr/ACGTN/FFFF!/;
        print "\@$i\n$seq1+\n$qual1";

        my $qual2 = $seq2;
        $qual2 =~ tr/ACGTN/FFFF!/;
        print "\@$i\n$seq2+\n$qual2";
        
        $i++;
    }
}
