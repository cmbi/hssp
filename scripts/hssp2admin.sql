-- SQL script for hssp2admin system

DROP TABLE pdb_chain CASCADE;
DROP TABLE chain_seq CASCADE;
DROP SEQUENCE chain_seq_id;
CREATE SEQUENCE chain_seq_id;

CREATE TABLE chain_seq
(
	id			integer CONSTRAINT chain_id_ix PRIMARY KEY DEFAULT nextval('chain_seq_id'),
	sequence	text NOT NULL,
	crc			char(16) NOT NULL,
	alignment	text
);

CREATE INDEX chain_crc_ix ON chain_seq(crc);

CREATE TABLE pdb_chain
(
	code		char(4),
	chain		char,
	seq			integer REFERENCES chain_seq(id) ON DELETE RESTRICT,
	cluster		integer REFERENCES chain_seq(id) ON DELETE RESTRICT,
	
	PRIMARY KEY (code, chain)
);

CREATE INDEX pdb_chain_seq_ix ON pdb_chain(seq);
