from Bio import SeqIO
from sys import argv
from sqlalchemy import create_engine

file = argv[1]

# for record in SeqIO.parse(open(file, 'r'), 'fasta'):
# 	print record.id.split('|')[0], record.seq


engine = create_engine('sqlite:///../results_db/results.db')
conn = engine.connect()

for record in SeqIO.parse(open(file, 'r'), 'fasta'):
	definiton = ' '.join(record.description.split('Rhodococcus')[0]
		           .rstrip('[ ').split()[2:])

	id = record.id
	peptide_seq = str(record.seq)
	species = "Rhodococcus equi 103S"
	# print record.id, record.name, record.description

	# conn.execute("""UPDATE genes SET coding_seq='%s'
 #                WHERE transcript='%s'""" % (cds, gene_id))
	conn.execute("""INSERT INTO sequences
				(id, definition, species, peptide_seq)
                VALUES ("%s", "%s", "%s", "%s")""" %
                (id, definiton, species, peptide_seq))
	print "Inserted %s" % id
