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
my ($date, $pdbid, $header, $nchain, $kchain, @query, $seqlen);

$seqlen = $nchain = $kchain = 0;

open(my $h, "<$file") or die "Could not open file $file for reading\n";
while (my $line = <$h>)
{
	chomp($line);
	die "Not a valid input file\n" unless ($line eq '# STOCKHOLM 1.0');
	
	$nchain += 1;
	$header = '';
	
	my $chainId;
	
	my $state = 1;
	for (;;)
	{
		my $line = <$h>;
		chomp($line);
		last if $line eq '//';
		next if length($line) == 0;
		
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
				$header = length($header) > 0 ? "$header\n$1     $2" : "$1     $2";
			}
		}
		elsif (substr($line, 0, 8) eq '#=GF RI ')
		{
			my %ri;
			
			$ri{info} = substr($line, 7, 6) . substr($line, 14, 33) . substr($line, 48, 5) . '  ' . substr($line, 54, 3);
			$ri{gap} = substr($ri{info}, 12, 1) eq '!' ? 1 : 0;

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
			if ($key eq 'CC') {
				$chainId = $id;
			}
			else {
				unless (defined $hits{$id}) {
					push @hits, $id;
					my %ins = ();
					$hits{$id}{ins} = \%ins;
					$hits{$id}{jpos} = 1;
				}
				$hits{$id}{$key} = $value;
			}
		}
		elsif (substr($line, 0, 1) ne '#')
		{
			my ($id, $data) = $line =~ m/^(\S+)\s+(\S+)$/;
			
			if ($id eq $chainId)
			{
				push @query, split(m//, $data);
			}
			else
			{
				die "Unexpected line '$line' (id = $id, data = $data)\n" unless defined $hits{$id};
				push @{$hits{$id}{aln}}, split(m//, $data);
			}
		}
	}
#	
#	print Dumper(\%hits);
}
close($h);

# 'cut' out the gaps 

my $gap = 0; my $ipos = 1;
for (my $i = 0; $i < scalar(@query); ++$i)
{
	my $r = $query[$i];
	if ($r ne '.')
	{
		++$seqlen;

		if ($gap)
		{
			foreach my $hit (values %hits)
			{
				next unless defined $hit->{ins}{$gap};
				
				my %insertion = (
					ipos => $ipos,
					jpos => $hit->{jpos} - length($hit->{ins}{$gap}),
					iseq => lc(@{$hit->{aln}}[$gap - 1]) . $hit->{ins}{$gap} . lc(@{$hit->{aln}}[$i])
				);
				
				push @{$hit->{insertions}}, \%insertion;

				@{$hit->{aln}}[$gap - 1] = lc @{$hit->{aln}}[$gap - 1];
				@{$hit->{aln}}[$i] = lc @{$hit->{aln}}[$i];
			}
		}
		
		$gap = 0;
		++$ipos;
	}
	elsif ($gap == 0)
	{
		$gap = $i;
	}
	
	foreach my $hit (values %hits)
	{
		if (@{$hit->{aln}}[$i] ne '.')
		{
			++$hit->{jpos};
			if ($gap > 0)
			{
				$hit->{ins}{$gap} .= @{$hit->{aln}}[$i];
			}
		}
	}
}

# clean up hit alignments
foreach my $hit (values %hits)
{
	my $started = 0; my $last = 0;
	for (my $i = 0; $i < scalar(@{$hit->{aln}}); ++$i)
	{
		if (@{$hit->{aln}}[$i] ne '.') {
			$started = 1;
			$last = $i;
		}
		elsif (not $started) {
			@{$hit->{aln}}[$i] = ' ';
		}
	}
	
	for (my $i = $last; $i < scalar(@{$hit->{aln}}); ++$i)
	{
		@{$hit->{aln}}[$i] = ' ';
	}
}

$seqlen = sprintf("%4d", $seqlen);
$nchain = sprintf("%4d", $nchain);
my $nalign = sprintf("%4d", scalar keys %hits);

print<<EOF;
HSSP       HOMOLOGY DERIVED SECONDARY STRUCTURE OF PROTEINS , VERSION 2.0 2011
PDBID      $pdbid
DATE       file generated on $date
SEQBASE    UniProt
THRESHOLD  according to: t(L)=(290.15 * L ** -0.562) + 5
REFERENCE  Sander C., Schneider R. : Database of homology-derived protein structures. Proteins, 9:56-68 (1991).
CONTACT    Maintained at http://www.cmbi.umcn.nl/ <hssp.cmbi\@radboudumc.nl>
$header
SEQLENGTH  $seqlen
NCHAIN     $nchain chain(s) in $pdbid data set
NALIGN     $nalign
WARNING  : The information in this file is generated from another file.
WARNING  : It is possible that information in this file is truncated to
WARNING  : accomodate the limited formatting capabilities. Please consult
WARNING  : the original file to see the full information.
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
format FMT1 =
@#### : @<<<<<<<<<<<@<<<   @#.## @#.##@####@####@####@####@####@####@####@####  @<<<<<<<<<<@*
$nr, $id, $strid, $ide, $sim, $ifir, $ilas, $jfir, $jlas, $lali, $ngap, $lgap, $lseq2, $acc, $desc
.

	my $ofh = select(STDOUT);
	$~ = 'FMT1';
	select($ofh);

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

# print alignments

my $N = scalar(keys %hits);

for (my $i = 0; $i < $N; $i += 70)
{
	my $last = $i + 70;
	$last = $N unless $last <= $N;
	
	printf("## ALIGNMENTS %4d - %4d\n SeqNo  PDBNo AA STRUCTURE BP1 BP2  ACC NOCC  VAR  ", $i + 1, $last);
	for (my $a = $i / 10; $a < $i / 10 + 7; ++$a) {
		printf("....:....%d", ($a + 1) % 10);
	}
	print "\n";

	my $resnr = 0;
	foreach my $res (@residues)
	{
		my $aln = '';
		if (not $res->{gap})
		{
			$resnr += 1 while $query[$resnr] eq '.';
			
			for (my $j = $i; $j < $last; ++$j)
			{
				$aln .= @{$hits{$hits[$j]}->{'aln'}}[$resnr];
			}
			++$resnr;
		}
		
		print $res->{info}, "  $aln\n";
	}
}

# print profile

print<<EOF;
## SEQUENCE PROFILE AND ENTROPY
 SeqNo PDBNo   V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D  NOCC NDEL NINS ENTROPY RELENT WEIGHT
EOF

print join("\n", @profile), "\n";

print<<EOF;
## INSERTION LIST
 AliNo  IPOS  JPOS   Len Sequence
EOF

for (my $nr = 1; $nr <= scalar(@hits); ++$nr)
{
	my $hit = $hits{$hits[$nr - 1]};

	next unless defined $hit->{insertions};

	my ($ipos, $jpos, $len, $iseq);

format FMT2 = 
@##### @#### @#### @#### ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$nr,   $ipos,$jpos,$len, $iseq
~~                       ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         $iseq
.
	my $ofh = select(STDOUT);
	$~ = 'FMT2';
	select($ofh);

	foreach my $ins (@{$hit->{insertions}})
	{
		$ipos = $ins->{ipos};
		$jpos = $ins->{jpos};
		$iseq = $ins->{iseq};
		$len = length($iseq) - 2;

		write;
	}
}

print "//\n";

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
