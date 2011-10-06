# makefile to do the work...

HOMEDIR = $(HOME)/
RESULTS = $(HOMEDIR)results/
SCRATCH = $(TMPDIR)/

#$(SCRATCH)uniref100.fasta: $(HOMEDIR)uniref100.fasta.gz
#	zcat $< > $@

$(SCRATCH)%.fa: $(HOMEDIR)dssp-fa.tbz
	cd $(SCRATCH) && tar xf $(HOMEDIR)dssp-fa.tbz ./$(@F)

$(SCRATCH)%.sto: $(SCRATCH)%.fa
	$(HOMEDIR)bin/jackhmmer -o /dev/null -N5 --noali -A $@ $< $(SCRATCH)uniref100.fasta
	rm $<

$(RESULTS)%.sto.bz2: $(SCRATCH)%.sto
	bzip2 -c $< > $@
	rm $<

$(RESULTS)%.done: $(RESULTS)%.sto.bz2
	scp $< dba@131.174.88.34:/data/hssp2/sto/
	rm $<
	touch $@
