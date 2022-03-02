use strict;
use warnings;

#----------------------------------------------------------
# umi loading and parsing
#----------------------------------------------------------

# working variables
our (%expectedUmis, %allowedUmis);
our ($umiLen, $nExpectedUmi, $nAllowedUmi) = (0) x 10;

# load UMIs with known values and fixed-widths
sub loadFixedUmis {
    my ($UMI_FILE, $UMI_SKIP_BASES) = @_;
    $UMI_FILE or return;
    $UMI_SKIP_BASES or $UMI_SKIP_BASES = 1; # default to a single-base 3' A extension typical of TruSeq

    # read fixed UMIs from their definition file
    open my $inH, "<", $UMI_FILE or throwError("could not open $UMI_FILE: $!");
    my $discardHeader = <$inH>; # definition file must have a header line
    while (my $line = <$inH>){
        $nExpectedUmi++;    
        chomp $line;
        $line =~ s/\r//g;
        my @line = split("\t", $line);
        my $umi = $line[0]; # the known, i.e., expected UMI value, must be in first column
        $expectedUmis{$umi}++;
        my @umi = split("", $umi);
        $umiLen = @umi;

        # allow one base mismatch from known UMIs
        foreach my $i(0..$#umi){ 
            my @allowed = @umi;
            foreach my $base(qw(A C G T)){
                $allowed[$i] = $base;
                $allowedUmis{join("", @allowed)} = $nExpectedUmi; # key=UMI as observed, value=inferred source UMI as 1-referenced index
                $nAllowedUmi++;
            }
        }
    }
    close $inH;

    # report fixed UMI metadata
    printCount($umiLen,       'umiLen',       'fixed UMI length');
    printCount($nExpectedUmi, 'nExpectedUmi', 'expected UMI values');
    printCount($nAllowedUmi,  'nAllowedUmi',  'allowed UMI values (up to 1 base mismatch from expected)');  

    # return regex to parse incoming reads
    qr/^(.{$umiLen}).{$UMI_SKIP_BASES}(.+)/;
}

1;
