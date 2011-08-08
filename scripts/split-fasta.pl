#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my ($in, $dir) = @ARGV;
die "Usage: split-fasta.pl <fasta-file> <dest-dir>\n"
	unless -f $in and -d $dir;

open(my $ih, "<$in") or die "Could not open $in: $!\n";
my ($id, $fa);

while (my $line = <$ih>)
{
	if ($line =~ m/^>(\S+)/) {
		&test_and_replace_fa($fa, $id) if defined $fa and defined $id;
		
		$id = $1;
		$fa = $line;
	}
	else {
		die "missing ID line\n" unless defined $id;
		$fa .= $line;
	}
}
&test_and_replace_fa($fa, $id) if defined $fa and defined $id;

close $ih;

# Now process the cluster file for this FastA

open($ih, "<${in}.clstr") or die "Cluster file is missing for $in\n";
my (%tbl, $main, @ids);
while (my $line = <$ih>) {
	if ($line =~ m/^>Cluster/) {
		if (defined $main) {
			foreach my $id (@ids) {
				push @{$tbl{$id}}, $main;
			}

			$main = undef;
			undef @ids;
		}
	}
	else {
		die "malformed line $line\n"
			unless $line =~ m/\d+\t\S+, >((\S{4})-\d+)\.{3} (\*|at 100.00%)/;
		$main = $1 if $3 eq '*';
		push @ids, $2;
	}
}
close($ih);

# print out dependancy file

open(my $oh, ">$dir/hssp2.depends") or die "Could not write dependancy file\n";
foreach my $id (sort keys %tbl) {
	print $oh "\$(HSSP2DIR)$id.hssp.bz2: ", join(" \\\n\t", map {"\$(HSSP2DATADIR)$_.sto.bz2"} @{$tbl{$id}}), "\n\n";
}
close($oh);

sub test_and_replace_fa
{
	my ($fa, $id) = @_;
	
	my $file = "$dir/${id}.fa";
	if (-f $file) {
		open(my $h, "<$file") or die "Could not open file $file\n";
		local($/) = undef;
		my $test = <$h>;
		close($h);
		
		return if $test eq $fa;
	}
	
	open(my $h, ">$file") or die "Could not create FastA file $file\n";
	print $h $fa;
	close($h);
}

