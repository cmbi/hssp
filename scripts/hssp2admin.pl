#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use Data::Dumper;
use File::stat;

my $SEQ_DB = '/data/fasta/uniprot.fa';

my $action = shift;

# ---------------------------------------------------------------------

my $MINIMAL_SEQUENCE_LENGTH = 25;

my $dbh = DBI->connect("dbi:Pg:dbname=hssp2admin", "", "")
	or die "Geen connectie databank\n";

my $fetch_seq_stmt = $dbh->prepare(qq{
	SELECT * FROM chain_seq WHERE crc = ?
});

my $fetch_seq_stmt2 = $dbh->prepare(qq{
	SELECT * FROM chain_seq WHERE id = ?
});

my $insert_seq_stmt = $dbh->prepare(qq{
	INSERT INTO chain_seq (sequence, crc) VALUES (?, ?)
});

my $update_seq_aln_stmt = $dbh->prepare(qq{
	UPDATE chain_seq SET alignment = ? WHERE id = ?
});

my $fetch_pdb_chain_stmt = $dbh->prepare(qq{
	SELECT * FROM pdb_chain WHERE code = ? and chain = ?
});

my $insert_pdb_chain_stmt = $dbh->prepare(qq{
	INSERT INTO pdb_chain (code, chain, seq, cluster) VALUES (?, ?, ?, ?)
});

my $count_pdb_id_stmt = $dbh->prepare(qq{
	SELECT COUNT(*) FROM pdb_chain WHERE code = ?
});

my $fetch_not_ready_stmt = $dbh->prepare(qq{
	SELECT DISTINCT code FROM pdb_chain WHERE seq IN (SELECT id FROM chain_seq WHERE alignment IS NULL)
});

my $fetch_pdb_chains_stmt = $dbh->prepare(qq{
	SELECT seq, MIN(chain) AS chain FROM pdb_chain WHERE code = ? GROUP BY seq
}) or die "Fout: $dbh->errstr\n";

my $update_cluster_stmt = $dbh->prepare(qq{
	UPDATE pdb_chain SET cluster = ? WHERE seq = ?
});

# ---------------------------------------------------------------------

if ($action eq 'init')
{
	&ProcessDSSP();
}
elsif ($action eq 'update')
{
	&ProcessUpdate();
}
elsif ($action eq 'aln')
{
	&ProcessAln();
}
elsif ($action eq 'cluster')
{
	&ProcessCluster();
}
elsif ($action eq 'make')
{
#	&CreateHSSPFiles();
	&CreateHSSPMakefile();
}
elsif ($action eq 'no-hits')
{
	&ProcessNoHits();
}
else
{
	# assume id is a HSSP id
	&CreateHSSPFile($action);
}
exit;

# ---------------------------------------------------------------------

sub AddSequence
{
	my ($seq) = @_;
	
	my $crc = crc64($seq);
	my $id;
	
	$fetch_seq_stmt->execute($crc) or die "fout in execute: " . $dbh->errstr . "\n";
	while (my $r = $fetch_seq_stmt->fetchrow_hashref())
	{
		if ($r->{sequence} eq $seq)
		{
			$id = $r->{id};
			last;
		}
	}
	
	unless (defined $id) {
		$insert_seq_stmt->execute($seq, $crc) or die "Fout: " . $dbh->errstr . "\n";
		$id = &AddSequence($seq);
	}
	
	return $id;
}

sub AddChain($$$)
{
	my ($pdbid, $chainID, $sequence) = @_;
	
	$fetch_pdb_chain_stmt->execute($pdbid, $chainID) or die "Fout: " . $dbh->errstr . "\n";
	
	if (not $fetch_pdb_chain_stmt->fetchrow_hashref()) {
		my $id = &AddSequence($sequence);
		die "Geen sequence ID?\n" unless defined $id;
		$insert_pdb_chain_stmt->execute($pdbid, $chainID, $id, $id) or die "Fout: " . $dbh->errstr . "\n";
	}
}

sub ReadDSSP($)
{
	my ($file) = @_;

	# Read DSSP file
	my (%dssp, $pdbid);
	
	open(my $ih, "<$file") or die "Could not open input file $file\n";
	while (my $line = <$ih>) {
		$pdbid = lc $1 if ($line =~ m/^HEADER  .{54}(\d...)/);
		last if ($line =~ m/^  #  RESIDUE/);
	}
	die "No ID found in DSSP file\n" unless defined $pdbid;
	
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
		my $seq = $dssp{$ch};
		$seq =~ s/X+$//;	# strip off trailing X residues
		$seq =~ s/-//g;		# remove gaps
		next unless length($seq) >= $MINIMAL_SEQUENCE_LENGTH;
		AddChain($pdbid, $ch, $seq);
	}
}

sub ListDSSPFiles()
{
	opendir(my $dh, "/data/dssp") or die "Could not read directory\n";
	my @files = grep { m/\.dssp$/ } readdir($dh);
	closedir($dh);
	return @files;
}

sub ProcessDSSP()
{
	my ($dir) = @_;
	
	foreach my $file (&ListDSSPFiles()) {
		print "$file\n";
		&ReadDSSP("/data/dssp/$file") if -f "/data/dssp/$file" and $file =~ m/\.dssp$/;
	}
}

# ---------------------------------------------------------------------

sub ProcessUpdate
{
#	# update the databanks using the make script
#	system("/data/bin/update-databanks pdb dssp uniprot /data/fasta/dssp-nr-100.fa");
	
#	system("rsync -a --delete ftp.wwpdb.org::ftp_data/structures/divided/pdb/ /data/pdb/divided");
	
	my %pdb;
	opendir(my $d, "/data/pdb/divided/") or die "opendir: $!\n";
	foreach my $subdir (readdir($d)) {
		next unless -d "/data/pdb/divided/$subdir" and length($subdir) == 2;
		opendir(my $sd, "/data/pdb/divided/$subdir") or die "opendir: $!\n";
		foreach my $pdb (readdir($sd)) {
			next unless $pdb =~ m/pdb(.{4})\.ent.gz/;
			my $pdbid = $1;
			$pdb{$pdbid} = { pdbgz => "/data/pdb/divided/$subdir/$pdb", pdb => "/data/pdb/flat/pdb${pdbid}.ent" };
		}
		closedir($sd);
	}
	closedir($d);

	open(my $h, "</data/status/dssp_skip.txt");
	my %skip = map { chomp; $_ => 1 } <$h>;
	close($h);
	
	foreach my $pdbid (sort keys %pdb) {
		next if $skip{$pdbid};
		
		unless (-f $pdb{$pdbid}{pdbgz}) {
			my $cmd = "gzcat $pdb{$pdbid}{pdbgz} > $pdb{$pdbid}{pdb}";
			print "$cmd\n";
			system($cmd);
		}

		unless (-f "/data/dssp/${pdbid}.dssp") {
			my $cmd = "/home/staf/maarten/projects/mas/mkdssp $pdb{$pdbid}{pdbgz} /data/dssp/${pdbid}.dssp";
			print "$cmd\n";
			system($cmd);
			if (-f "/data/dssp/${pdbid}.dssp") {
				&ReadDSSP("/data/dssp/${pdbid}.dssp");
			}
			else {
				$skip{$pdbid} = 1;
			}
		}
	}
	
	open($h, ">/data/status/dssp_skip.txt") or die "open: $!\n";
	print $h join("\n", sort keys %skip);
	close($h);
	
	# run cd-hit
	
	open(my $dh, ">/data/fasta/dssp.fa") or die "open: $!\n";
	my $stmt = $dbh->prepare(qq{SELECT * FROM chain_seq}) or die "prepare: ".$dbh->errstr."\n";
	$stmt->execute() or die "execute: ".$dbh->errstr."\n";
	while (my $r = $stmt->fetchrow_hashref()) {
		my $id = "seq-".$r->{id};
		my $seq = $r->{sequence};
		$seq =~ s/.{72}/$&\n/g;
		print $dh ">$id\n$seq\n";
	}
	close($dh);
	
#	system("/usr/local/bin/cd-hit -i /data/fasta/dssp.fa -o /data/fasta/dssp-nr-100.fa -c 1.0 -T 0 -s 0.9");
	
	# We should have a cluster file now
	
	open($h, "</data/fasta/dssp-nr-100.fa.clstr") or die "cluster? $!\n";
	my (@clustered, $representative);
	while (my $line = <$h>) {
		if (substr($line, 0, 8) eq '>Cluster') {
			foreach my $id (@clustered) {
				$update_cluster_stmt->execute($representative, $id)
					or die "fout: ".$dbh->errstr."\n";
			}
			@clustered = ();
		}
		else {
			die "$line" unless $line =~ m/\d+\s+\d+aa, >seq-(\d+)\.\.\. (\*|at 100.00%)/;
			push @clustered, $1;
			$representative = $1 if $2 eq '*';
		}
	}
}


# ---------------------------------------------------------------------
# update the chain_seq table to have 'no hits' for all the missing alignment files

sub ProcessNoHits
{
	my $stmt = $dbh->prepare(qq{SELECT * FROM chain_seq WHERE alignment != 'no hits'});
	
	my @update;
	
	$stmt->execute() or die "execute: ".$dbh->errstr."\n";
	while (my $r = $stmt->fetchrow_hashref()) {
		push @update, $r->{id} if (not -f $r->{alignment});
	}
	
	foreach my $id (@update) {
		$update_seq_aln_stmt->execute('no hits', $id) or die "Fout: " . $dbh->errstr . "\n";
	}
	
	printf("Updated %d records\n",  scalar @update);
}

# ---------------------------------------------------------------------

sub CreateFasta
{
	my ($id) = @_;
	
	$fetch_seq_stmt2->execute($id) or die "Fout: " . $dbh->errstr . "\n";
	my $r = $fetch_seq_stmt2->fetchrow_hashref() or die "Seq $id niet gevonden: " . $dbh->errstr . "\n";
	my $seq = $r->{sequence} or die "Geen sequentie voor $id?\n";
	$seq =~ s/.{72}/$&\n/g;

	my $file = "/data/hssp2/fasta/seq-$id.fa";
	if (defined $r->{alignment} and $r->{alignment} =~ m|/data/hssp2/sto/(.+)\.aln\.bz2|) {
		$file = "/data/hssp2/fasta/$1.fa";
	}

	unless (-f $file) {
		open(my $h, ">$file") or die "Fout openen file: $!\n";
		print $h ">seq-$id\n$seq\n";
		close($h);
	}
	
	return $file;
}

sub CreateAlignment($)
{
	my ($file, $seqid) = @_;
	
	unless (-f $file) {
		my $seqfile = &CreateFasta($seqid);

		my $tmpdir = "/srv/data/tmp/$seqid";
		if (not -d $tmpdir) {
			mkdir($tmpdir) or die "mkdir: $!\n";
		}
		chdir($tmpdir);
		open(my $p, "/usr/local/bin/jackhmmer --cpu 32 -N5 --noali -A ${seqid}.sto --chkhmm ${seqid} $seqfile $SEQ_DB|")
			or die "jackhmmer: $!\n";
		while (my $line = <$p>) {
		}
		close($p);
		
		my $hfile = $file;
		$hfile =~ s/\.aln\.bz2$//;
		foreach my $n (1 .. 5) {
			rename("${seqid}-${n}.hmm", "${hfile}-${n}.hmm") if -f "${seqid}-${n}.hmm";
		}
		system("/home/staf/maarten/projects/mas/sto2fa ${seqid}.sto $file");
		unless (-f $file) {
			warn "Could not create $file\n";
			
			$update_seq_aln_stmt->execute('no hits', $seqid) or die "Fout: " . $dbh->errstr . "\n";
		}
		chdir("/");
		rmdir($tmpdir);
	}
}

sub CreateHSSPFile
{
	my ($id) = @_;
	
	my $file = "/data/hssp/$id.hssp.bz2";
	my %mapped;
	
	$fetch_pdb_chains_stmt->execute($id) or die "fout in execute: " . $dbh->errstr . "\n";
	while (my $r = $fetch_pdb_chains_stmt->fetchrow_hashref())
	{
		$fetch_seq_stmt2->execute($r->{cluster});
		my $r2 = $fetch_seq_stmt2->fetchrow_hashref() or die "Geen chain?: ".$dbh->errstr."\n";
		
		my $alignment = $r2->{alignment};
		next if $alignment and $alignment eq 'no hits';
		
		unless ($alignment) {
			$alignment = sprintf("/data/hssp2/sto/%s-%s.aln.bz2",
				$r2->{crc}, $r2->{id});
			$update_seq_aln_stmt->execute($alignment, $r2->{id}) or die "Fout: " . $dbh->errstr . "\n";
		}
		
		&CreateAlignment($alignment, $r2->{id});
		$mapped{$r->{chain}} = $alignment;
	}
	
	next unless scalar(keys %mapped) > 0;
	
	my $cmd = sprintf("/home/staf/maarten/projects/mas/aln2hssp /data/pdb/flat/pdb${id}.ent %s $file",
		join(' ', map { " --chain $_=$mapped{$_}" } sort keys %mapped));
	system($cmd);

	die "Fout: $file?\n" unless -f $file;
}

sub CreateHSSPFiles()
{
	my %hssp = map {
		my $pdbid = $1 if m/(....)\.dssp$/;
		my %e = ( hssp => "/data/hssp/${pdbid}.hssp.bz2", dssp => "/data/dssp/$_" );
		$pdbid => \%e;
	} &ListDSSPFiles();
	
	my %skip;
	$fetch_not_ready_stmt->execute() or die "Fout: " . $dbh->errstr . "\n";
	while (my $r = $fetch_not_ready_stmt->fetchrow_hashref()) {
		$skip{$r->{code}} = 1;
	}
	
	foreach my $id (keys %hssp) {
		my $file = $hssp{$id}{hssp};

		if ($skip{$id}) {
			print "skip: $id\n";
			next;
		}

		next unless -f $file;			# rebuild
		next if -M $file < 14;			# only old ones

		print "$id\n";
		
		&CreateHSSPFile($id);
	}
}

# ---------------------------------------------------------------------

sub WriteAlignmentRule
{
	my ($alnfile, $seqid, $h) = @_;
	
	my $seqfile = &CreateFasta($seqid);

	my $tmpdir = "/srv/data/tmp/$seqid";
	my $hfile = $alnfile;
	$hfile =~ s/\.aln\.bz2$//;
		
print $h <<EOF;

$alnfile: $seqfile
	test -d $tmpdir || mkdir $tmpdir
	cd $tmpdir && \$(JACKHMMER) -A ${seqid}.sto --chkhmm ${seqid} $seqfile $SEQ_DB
	\$(STO2FA) $tmpdir/${seqid}.sto \$@
	test -f ${tmpdir}/${seqid}-1.hmm && mv ${tmpdir}/${seqid}-1.hmm ${hfile}-1.hmm || true
	test -f ${tmpdir}/${seqid}-2.hmm && mv ${tmpdir}/${seqid}-2.hmm ${hfile}-2.hmm || true
	test -f ${tmpdir}/${seqid}-3.hmm && mv ${tmpdir}/${seqid}-3.hmm ${hfile}-3.hmm || true
	test -f ${tmpdir}/${seqid}-4.hmm && mv ${tmpdir}/${seqid}-4.hmm ${hfile}-4.hmm || true
	test -f ${tmpdir}/${seqid}-5.hmm && mv ${tmpdir}/${seqid}-5.hmm ${hfile}-5.hmm || true
	rm -rf $tmpdir

EOF
}

sub WriteHSSPRule
{
	my ($id, $ali, $h) = @_;
	
	my $file = "/data/hssp/$id.hssp.bz2";
	my %mapped;
	
	$fetch_pdb_chains_stmt->execute($id) or die "fout in execute ($id): " . $dbh->errstr . "\n";
	while (my $r = $fetch_pdb_chains_stmt->fetchrow_hashref())
	{
		$fetch_seq_stmt2->execute($r->{seq});
		my $r2 = $fetch_seq_stmt2->fetchrow_hashref() or die "Geen chain?: ".$dbh->errstr."\n";
		
		my $alignment = $r2->{alignment};

		next if $alignment and ($alignment eq 'no hits');
		
		unless ($alignment) {
			$alignment = sprintf("/data/hssp2/sto/%s-%s.aln.bz2",
				$r2->{crc}, $r2->{id});
			$update_seq_aln_stmt->execute($alignment, $r2->{id}) or die "Fout: " . $dbh->errstr . "\n";
		}

		$mapped{$r->{chain}} = $alignment;
		
		next if defined $ali->{$alignment};
		
		&WriteAlignmentRule($alignment, $r2->{id}, $h);
		$ali->{$alignment} = 1;
	}
	
	my $result = 0;
	if (scalar(keys %mapped) > 0) {
		my $cmd = sprintf("\$(ALN2HSSP) /data/pdb/flat/pdb${id}.ent %s $file",
			join(' ', map { "--chain $_=$mapped{$_}" } sort keys %mapped));
		my $aln = join(' ', values %mapped);

print $h <<EOF;

$file: /data/pdb/flat/pdb${id}.ent $aln
	$cmd || (echo $id >> hssp2-fail.txt)

EOF
		$result = 1;
	}
	
	return $result;
}

sub CreateHSSPMakefile()
{
	my %hssp = map {
		my $pdbid = $1 if m/(....)\.dssp$/;
		my %e = ( hssp => "/data/hssp/${pdbid}.hssp.bz2", dssp => "/data/dssp/$_" );
		$pdbid => \%e;
	} &ListDSSPFiles();
	
	open(my $h, ">hssp2.make") or die "makefile: $!\n";
	
	my @hsspids;

	print $h <<EOF;
JACKHMMER = /usr/local/bin/jackhmmer --cpu 32 -N 5 --noali -o /dev/null
ALN2HSSP  = /home/staf/maarten/projects/mas/aln2hssp
STO2FA    = /home/staf/maarten/projects/mas/sto2fa

firstTarget: all

EOF
	
	my %alignments = ();		# keep track of alignment rules
	foreach my $id (sort keys %hssp) {
		push @hsspids, $id if &WriteHSSPRule($id, \%alignments, $h);
	}
	
	my $hsspids = join(' ', @hsspids);

print $h <<EOF;

HSSPIDS   = $hsspids
HSSPFILES = \$(HSSPIDS:%=/data/hssp/%.hssp.bz2)

all: \$(HSSPFILES)

EOF
	
}

# ---------------------------------------------------------------------
#
#	ProcessAlnFile assumes first sequence in FastA formatted file is the query

sub ProcessAlnFile($)
{
	my ($file) = @_;
	
	open(my $h, "bzcat $file|") or die "open?\n";
	
	my $firstLine = <$h>;
	die "Leeg?\n" unless $firstLine;
	die "Geen fasta?\n" unless substr($firstLine, 0, 1) eq '>';
	
	my $seq = '';
	while (my $line = <$h>) {
		last if substr($line, 0, 1) eq '>';
		$seq .= $line;
	}
	
	$seq =~ s/\s+|-//g;
	my $crc = crc64($seq);
	my ($id, $alignment);
	
	$fetch_seq_stmt->execute($crc) or die "fout in execute: " . $dbh->errstr . "\n";
	while (my $r = $fetch_seq_stmt->fetchrow_hashref())
	{
		if ($r->{sequence} eq $seq)
		{
			$id = $r->{id};
			$alignment = $r->{alignment};
			last;
		}
	}
	
	if ($id) {
		if (not defined $alignment or $alignment ne $file) {
			$update_seq_aln_stmt->execute($file, $id) or die "Fout in update: " . $dbh->errstr . "\n";
		}
	}
	else {
		print "stale alignment file $file\nseq: $seq\n\n";
	}
	
	close($h);
}

sub ProcessAln()
{
	opendir(my $dh, "/data/hssp2/sto") or die "geen sto dir?\n";
	my @files = grep { m/\.aln\.bz2/ } readdir($dh);
	closedir($dh);
	
	foreach my $file (@files) {
		&ProcessAlnFile("/data/hssp2/sto/$file");
	}
}

# ---------------------------------------------------------------------
#

#sub ProcessCluster
#{
#	
#	
#	open(my $p = "/usr/local/bin/cd-hit -i /data/fasta/dssp.fa 
#}

# ---------------------------------------------------------------------
# copied from the swissknife toolkit

sub crc64
{
	my $POLY64REVh = 0xd8000000; 
	my @CRCTableh = 256;
	my @CRCTablel = 256;
	my $initialized;
	


  my $sequence = shift;
  my $crcl = 0;
  my $crch = 0;
  if (!$initialized) {
    $initialized = 1;
    for (my $i=0; $i<256; $i++) {
      my $partl = $i;
      my $parth = 0;
      for (my $j=0; $j<8; $j++) {
	my $rflag = $partl & 1;
	$partl >>= 1;
	$partl |= (1 << 31) if $parth & 1;
	$parth >>= 1;
	$parth ^= $POLY64REVh if $rflag;
      }
      $CRCTableh[$i] = $parth;
      $CRCTablel[$i] = $partl;
    }
  }
  
  foreach (split '', $sequence) {
    my $shr = ($crch & 0xFF) << 24;
    my $temp1h = $crch >> 8;
    my $temp1l = ($crcl >> 8) | $shr;
    my $tableindex = ($crcl ^ (unpack "C", $_)) & 0xFF;
    $crch = $temp1h ^ $CRCTableh[$tableindex];
    $crcl = $temp1l ^ $CRCTablel[$tableindex];
  }
  return wantarray ? ($crch, $crcl) : sprintf("%08X%08X", $crch, $crcl);
}

