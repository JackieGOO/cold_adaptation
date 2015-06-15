"""
format blastp output when the following -outfmt option was used:
-outfmt "6 qseqid sallseqid evalue bitscore length stitle"

the final output contains the following items per line:

Query ID, Target ID, Evalue, Target Accession, Target Definition, Species
E.g.:
2529298360 646584812   3e-40   ROP_27020   hypothetical protein    Rhodococcus opacus B4
"""

from sys import argv
import csv

with open(argv[1], 'r') as csvfile:
    blast_reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
    filename = "{}_parsed.tab".format(argv[1])
    with open(filename, 'w') as csvout:
        csvwriter = csv.writer(csvout, delimiter='\t', quotechar='|')
        for line in blast_reader:
            query_GI = line[0]
            target_GI = line[1]
            evalue = line[2]
            extra = line[5].split('Rhodococcus')
            target_accession = extra[0].split()[0]
            target_definition = ' '.join(extra[0].rstrip('[ ').strip('|')
                                   .split()[1:])
            species = "Rhodococcus{}".format(extra[1].split(':')[0]
                            .split('contig')[0].split('plasmid')[0]).strip()
            csvwriter.writerow([query_GI, target_GI, evalue, target_accession,
                          target_definition, species])