# Perl module voor het parsen van UniRef
#
package MRS::Script::hssp2ix;

use Data::Dumper;

our @ISA = "MRS::Script";

sub new
{
	my $invocant = shift;

	my $self = new MRS::Script(
		name		=> 'HSSP2IX',
		url			=> 'res://',
		section		=> 'other',
		charset		=> 'utf8',
		meta		=> [ 'title' ],
		indices		=> {
			'id'		=> { name => 'Identifier' },
		},
		raw_files	=> qr/dssp-nr-100.fa.clstr/,
		@_
	);
	
	return bless $self, "MRS::Script::hssp2ix";
}

sub Parse
{
	my $self = shift;

	my %mapped;
	
	open(my $h, "</data/fasta/dssp.fa") or die "Could not open dssp.fa file\n";
	while (my $line = <$h>) {
		chomp($line);
		if ($line =~ m/^>(\S+) (.+)/) {
			$mapped{$1} = $2;
		}
	}
	close($h);
	
	my ($doc);
	my $data = $self->file;

	my %data;

	while (my $line = <$data>)
	{
		next if ($line =~ m/^>Cluster (\d+)/);
		
		my ($nr, $l, $p, $t) = ($line =~ m/^(\d+)\s+(\d+)aa, >(....-\d+)... (\*|at .+)/);

		die "Missing id in map: $p\nline: $line\n\n" unless defined $mapped{$p};
		
		my @chains = sort split(m/-/, substr($mapped{$p}, 5));
		die "Missing PDB chain IDS in $line\n" unless scalar(@chains) > 0;
		
		my $pdb = substr($p, 0, 4);
		
		$data{$pdb} .= "$chains[0]=$p\n";
	}
	
	# Now create documents
	foreach my $pdb (sort keys %data) {
		$self->Store($data{$pdb});
		$self->IndexValue('id', $pdb);
		$self->FlushDocument;
	}
}

1;
