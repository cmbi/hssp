#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

# Collect options and validate them

my $maxhits = 5000;
my $help = 0;

GetOptions("max-hits=i" => \$maxhits, "help|?" => \$help)
	or pod2usage(2);
pod2usage(1) if $help;

my $file = shift;
pod2usage(1) unless $file and -f $file;

die "Max hits exceeds maximum value of 9999\n"
	unless $maxhits < 10000;

# Main program, read file and collect data

my (%hits, @hits, @residues, @profile, $nchains, $kchains);

open(my $h, "<$file") or die "Could not open file $file for reading\n";
while (my $line = <$h>)
{
	chomp($line);
	die "Not a valid input file\n" unless ($line eq '# STOCKHOLM 1.0');
	
	&readChain($h);
	
#	print Dumper(\%hits);
}
close($h);



print<<EOF;
HSSP       HOMOLOGY DERIVED SECONDARY STRUCTURE OF PROTEINS , VERSION 2.0 2011
PDBID      $id
DATE       file generated on $date
SEQBASE    UniProt
THRESHOLD  according to: t(L)=(290.15 * L ** -0.562) + 5
REFERENCE  Sander C., Schneider R. : Database of homology-derived protein structures. Proteins, 9:56-68 (1991).
CONTACT    Maintained at http://www.cmbi.ru.nl/ by Maarten L. Hekkelman <m.hekkelman@cmbi.ru.nl>
HEADER     $header
COMPND     $compnd
SOURCE     $source
AUTHOR     $author
SEQLENGTH  $seqlen
NCHAIN     $nchain chain(s) in $id data set
NALIGN     $nalign
NOTATION : ID: EMBL/SWISSPROT identifier of the aligned (homologous) protein
NOTATION : STRID: if the 3-D structure of the aligned protein is known, then STRID is the Protein Data Bank identifier as taken
NOTATION : from the database reference or DR-line of the EMBL/SWISSPROT entry
NOTATION : %IDE: percentage of residue identity of the alignment
NOTATION : %SIM (%WSIM):  (weighted) similarity of the alignment
NOTATION : IFIR/ILAS: first and last residue of the alignment in the test sequence
NOTATION : JFIR/JLAS: first and last residue of the alignment in the alignend protein
NOTATION : LALI: length of the alignment excluding insertions and deletions
NOTATION : NGAP: number of insertions and deletions in the alignment
NOTATION : LGAP: total length of all insertions and deletions
NOTATION : LSEQ2: length of the entire sequence of the aligned protein
NOTATION : ACCNUM: SwissProt accession number
NOTATION : PROTEIN: one-line description of aligned protein
NOTATION : SeqNo,PDBNo,AA,STRUCTURE,BP1,BP2,ACC: sequential and PDB residue numbers, amino acid (lower case = Cys), secondary
NOTATION : structure, bridge partners, solvent exposure as in DSSP (Kabsch and Sander, Biopolymers 22, 2577-2637(1983)
NOTATION : VAR: sequence variability on a scale of 0-100 as derived from the NALIGN alignments
NOTATION : pair of lower case characters (AvaK) in the alignend sequence bracket a point of insertion in this sequence
NOTATION : dots (....) in the alignend sequence indicate points of deletion in this sequence
NOTATION : SEQUENCE PROFILE: relative frequency of an amino acid type at each position. Asx and Glx are in their
NOTATION : acid/amide form in proportion to their database frequencies
NOTATION : NOCC: number of aligned sequences spanning this position (including the test sequence)
NOTATION : NDEL: number of sequences with a deletion in the test protein at this position
NOTATION : NINS: number of sequences with an insertion in the test protein at this position
NOTATION : ENTROPY: entropy measure of sequence variability at this position
NOTATION : RELENT: relative entropy, i.e.  entropy normalized to the range 0-100
NOTATION : WEIGHT: conservation weight

## PROTEINS : identifier and alignment statistics
  NR.    ID         STRID   %IDE %WSIM IFIR ILAS JFIR JLAS LALI NGAP LGAP LSEQ2 ACCNUM     PROTEIN
EOF

if (1)
{

my ($nr, $id, $strid, $ide, $sim, $ifir, $ilas, $jfir, $jlas, $lali, $ngap, $lgap, $lseq2, $acc, $desc);
#   1 : CRAM_CRAAB          0.98  1.00    2   46    2   46   45    0    0   46  P01542     Crambin;
format =
@#### : @<<<<<<<<<<<@<<<   @#.## @#.##@####@####@####@####@####@####@####@####  @<<<<<<<<<<@*
$nr, $id, $strid, $ide, $sim, $ifir, $ilas, $jfir, $jlas, $lali, $ngap, $lgap, $lseq2, $acc, $desc
.

	for (my $i = 0; $i < scalar @hits; ++$i)
	{
		$nr = $i + 1;
		$id = $hits[$i];
		my $hit = $hits{$id};
		
		($ide, $sim) = $hit->{'HSSP'} =~ m|score=(\d+\.\d+)/(\d+\.\d+)|;
		($ifir, $ilas, $jfir, $jlas) = $hit->{'HSSP'} =~ m|aligned=(\d+)-(\d+)/(\d+)-(\d+)|;
		($lali) = $hit->{'HSSP'} =~ m|length=(\d+)|;
		($ngap) = $hit->{'HSSP'} =~ m|ngaps=(\d+)|;
		($lgap) = $hit->{'HSSP'} =~ m|gaplen=(\d+)|;
		($lseq2) = $hit->{'HSSP'} =~ m|seqlen=(\d+)|;
		
		$strid = undef;
		print Dumper($hit->{'DR'}) if defined $hit->{'DR'};
		if (defined $hit->{'DR'} and substr($hit->{'DR'}, 0, 4) eq 'PDB ') {
			my @strids = split(m/, /, substr($hit->{'DR'}, 4));
			$strid = uc shift @strids;
		}
		$strid = '' unless defined $strid;
		
		$acc = $id;
		$acc =~ s|/.+$||;
		$id = $hit->{'ID'};
		$desc = $hit->{'DE'};
		
		write;
	}
}

sub readChain
{
	my ($h) = @_;
	
	my ($date, $pdbid, $header);
	$header = '';
	
	my $state = 1;
	for (;;)
	{
		my $line = <$h>;
		chomp($line);
		last if $line eq '//';
		
		if (substr($line, 0, 8) eq '#=GF CC ')
		{
			$line = substr($line, 8);
			
			if ($line =~ m/^DATE\s+(.+)/) {
				$date = $1;
			}
			elsif ($line =~ m/^PDBID\s+(.+)/) {
				$pdbid = $1;
			}
			elsif ($line =~ m/^(HEADER|COMPND|SOURCE|AUTHOR) (.+)/) {
				$header .= $line . "\n";
			}
		}
		elsif (substr($line, 0, 8) eq '#=GF RI ')
		{
			my %ri;
			
			$ri{'info'} = substr($line, 7, 6) . substr($line, 14, 33) . substr($line, 48, 5) . '  ' . substr($line, 54, 3);
			
			push @residues, \%ri;
		}
		elsif (substr($line, 0, 8) eq '#=GF PR ')
		{
			my $pr = substr($line, 8, 5) . substr($line, 14);
			push @profile, $pr;
		}
		elsif (substr($line, 0, 5) eq '#=GS ')
		{
			my ($id, $key, $value) = substr($line, 5) =~ m/(\S+)\s+([A-Z]+) (.+)/;
			next if $key eq 'CC';

			push @hits, $id unless defined $hits{$id};
			$hits{$id}{$key} = $value;
		}
	}
}

__END__

=head1 NAME

hssp3to1.pl - Convert hssp version 3 files to original version 1 format

=head1 SYNOPSIS

hssp3to1.pl [options] file

 Options:
  --help       This help message
  --max-hits   Maximum number of hits to include

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exit.

=item B<--max-hits>

Specifies the maximum number of hits to include in the final file.
This number cannot exceed 9999. Note that in version 1 files, hits
for multiple chains are combined into one set.

=back

=head1 DESCRIPTION

B<hssp3to1.pl> will read a version 3 HSSP file and transform it into
a version 1 file.

=cut
