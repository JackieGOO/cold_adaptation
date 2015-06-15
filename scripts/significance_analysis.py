from Bio import SeqIO
from sqlalchemy import create_engine
import numpy as np
from scipy import stats

# load all sequences for quick access
prot_seqs = SeqIO.index("../all_proteins.faa", "fasta")


# connect to the database
engine = create_engine('sqlite:///../results_db/results.db')
conn = engine.connect()

# get the list of Rhodococcus sp. JG3 genes that have BLAST hits
JG3_genes = conn.execute("""SELECT DISTINCT query_id FROM aa_analysis""")


def aliphatic_comp(gene_list, hit_list):
    """ returns a letter for the one sample t-test of aliphatic index:
    Z - all hit aliphatic indeces are the same, therefore can't runt the test
    S - no significant difference
    A - psychrophilic query is significantly different and higher than targets
    B - psychrophilic query is significantly different and lower than targets
    """
    aliphatic_index_q = gene_list[2]
    aliphatic_index_t = [i[2] for i in hit_list]

    # check if all the indeces are the same, in that case the standard
    # deviation would be 0 and would result in division by 0.
    if len(set(aliphatic_index_t)) == 1:
        return 'Z'
    with np.errstate(divide='ignore'):
        p_value = stats.ttest_1samp(aliphatic_index_t, aliphatic_index_q)[1]
    if p_value <= 0.05:
        if aliphatic_index_q > np.mean(aliphatic_index_t):
            return 'A'
        else:
            return 'B'
    else:
        return 'S'


def aa_percent_comp(gene_list, hit_list, aa):
    """ returns a letter for the one sample t-test of amino acid percent:
    Z - all hit amino acid percents are the same, can't runt the test
    S - no significant difference
    A - psychrophilic query is significantly different and higher than targets
    B - psychrophilic query is significantly different and lower than targets
    """
    aa_dict = {"A": 3, "C": 4, "D": 5, "E": 6, "F": 7, "G": 8, "H": 9,
               "I": 10, "K": 11, "L": 12, "M": 13, "N": 14, "P": 15, "Q": 16,
               "R": 17, "S": 18, "T": 19, "V": 20, "W": 21, "Y": 22}
    aa_q = gene_list[3]
    aa_t = [i[3] for i in hit_list]
    if len(set(aa_t)) == 1:
        return 'Z'
    with np.errstate(divide='ignore'):
        p_value = stats.ttest_1samp(aa_t, aa_q)[1]
    if p_value <= 0.05:
        if aa_q > np.mean(aa_t):
            return 'A'
        else:
            return 'B'
    else:
        return 'S'


def arg_lys_ratio_comp(gene_list, hit_list):
    """ returns a letter for the one sample t-test of RK ratio:
    Z - all hit RK ratios are the same, can't runt the test
    S - no significant difference
    A - psychrophilic query is significantly different and higher than targets
    B - psychrophilic query is significantly different and lower than targets
    N - RK ratio for query couldn't be calculated due to K = 0 for query
    """
    rk_q = gene_list[43]
    if rk_q != 'N/A':  # check if RK ratio was not calculated due to K = 0
        rk_q = float(rk_q)
        # remove all the RK ratios where there were no lysines (i.e. N/A)
        rk_t = [float(i[43]) for i in hit_list if i[43] != 'N/A']
        # if rest of the numbers are the same (1), or there are less than 3
        # numbers (2) or all values are N/A and an empty list is returned (0)
        # then return Z
        if len(set(rk_t)) <= 2:
            return 'Z'
        with np.errstate(divide='ignore'):
            p_value = stats.ttest_1samp(rk_t, rk_q)[1]
        if p_value <= 0.05:
            if rk_q > np.mean(rk_t):
                return 'A'
            else:
                return 'B'
        else:
            return 'S'
    else:
        return 'N'


def acidic_res_comp(gene_list, hit_list):
    ac_q = int(gene_list[44])
    ac_t = [float(i[44]) for i in hit_list]
    if len(set(ac_t)) == 1:
        return 'Z'
    with np.errstate(divide='ignore'):
        p_value = stats.ttest_1samp(ac_t, ac_q)[1]
    if p_value <= 0.05:
        if ac_q > np.mean(ac_t):
            return 'A'
        else:
            return 'B'
    else:
        return 'S'


def aromaticity_comp(gene_list, hit_list):
    ar_q = float(gene_list[45])
    ar_t = [float(i[45]) for i in hit_list]
    if len(set(ar_t)) == 1:
        return 'Z'
    with np.errstate(divide='ignore'):
        p_value = stats.ttest_1samp(ar_t, ar_q)[1]
    if p_value <= 0.05:
        if ar_q > np.mean(ar_t):
            return 'A'
        else:
            return 'B'
    else:
        return 'S'


def gravy_comp(gene_list, hit_list):
    gr_q = float(gene_list[46])
    gr_t = [float(i[46]) for i in hit_list]
    if len(set(gr_t)) == 1:
        return 'Z'
    with np.errstate(divide='ignore'):
        p_value = stats.ttest_1samp(gr_t, gr_q)[1]
    if p_value <= 0.05:
        if gr_q > np.mean(gr_t):
            return 'A'
        else:
            return 'B'
    else:
        return 'S'

for gene in JG3_genes:
    gene = gene[0]
    gene_data = conn.execute("""SELECT * FROM aa_analysis WHERE query_id={0}
                                AND target_id={0}""".format(gene)).fetchall()
    hit_data = conn.execute("""SELECT * FROM aa_analysis WHERE query_id={0}
                               AND target_id!={0}""".format(gene)).fetchall()
    for result in gene_data:
        print gene,'\t', aliphatic_comp(result, hit_data),'\t',\
        arg_lys_ratio_comp(result, hit_data),'\t',\
        acidic_res_comp(result, hit_data),'\t',\
        aromaticity_comp(result, hit_data),'\t',\
        gravy_comp(result, hit_data)