#!/usr/bin/perl

use strict;
use warnings;
use English;
use Data::Dumper;

$OUTPUT_AUTOFLUSH = 1;

# some configuration globals
my $MINIMAL_SEQUENCE_LENGTH = 25;

# Fetch parameters and check for validity
my ($in, $out) = @ARGV;
die "Usage: mkhssp2 <indir> <outfile>\n"
	unless -d $in and defined $out;

# create a tmp file for the output
my $out_tmp = "$out-$$";
open(my $oh, ">$out_tmp") or die "Could not create FastA file\n";

opendir(my $dh, $in) or die "Could not read directory $in\n";
foreach my $file (sort readdir($dh)) {
	&process_dssp($file) if -f "$in/$file" and $file =~ m/\.dssp$/;
}
closedir($dh);
close($oh);

rename($out_tmp, $out) or die "Could not rename final FastA file\n";

exit;

sub process_dssp
{
	my ($file) = @_;

	# Read DSSP file
	my (%dssp, $id);
	
	open(my $ih, "<$in/$file") or die "Could not open input file $in/$file\n";
	while (my $line = <$ih>) {
		$id = lc $1 if ($line =~ m/^HEADER  .{54}(\d...)/);
		last if ($line =~ m/^  #  RESIDUE/);
	}
	die "No ID found in DSSP file\n" unless defined $id;
	
	my $break = 0;
	my $ch = '';
	
	while (my $line = <$ih>) {
		my $rch = substr($line, 11, 1);
		my $res = substr($line, 13, 1);
		$res = 'C' if $res =~ m/[a-z]/;
	
		if ($rch eq ' ' and $res eq '!') {
			$break = 1;
			next;
		}
		
		$dssp{$rch} = '' unless defined $dssp{$rch};
		$dssp{$rch} .= '-' if ($ch eq $rch and $break);
		$dssp{$rch} .= $res;
	
		$ch = $rch;
		$break = 0;
	}
	close($ih);
	
	my %seq;
	foreach my $ch (keys %dssp) {
		push @{$seq{$dssp{$ch}}{'chain'}}, $ch;
	}
	
	my $nr = 1;
	foreach my $s (keys %seq) {
		next unless length($s) >= $MINIMAL_SEQUENCE_LENGTH;
		my @chains = @{$seq{$s}->{'chain'}};
		my $fa = $s;
		$fa =~ s/.{72}/$&\n/g;
		print $oh ">$id-$nr $id-", join('-', @chains), "\n", $fa, "\n";
		++$nr;
	}
}
