import pandas as pd
import numpy as np
import data_cleanup as dc
from Bio.Seq import Seq
import severity_score as ses
import random

dna_sequence: str = 'ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGAC'
rna_sequence = dna_sequence.replace("T", "U")
p53_codons_kotler = [rna_sequence[i:i + 3] for i in range(0, len(rna_sequence), 3)]
p53_codons_gia = p53_codons_kotler
p53_codons_gia[71] = "CGC"


def process_dna_seq(seq: str) -> list:
    rna_seq = seq.replace("T", "U")
    res: list = [rna_seq[i:i+3] for i in range(len(seq), 3)]
    return res


def generate_codon_variations(codons: list) -> pd.DataFrame:
    variations: list = []
    for codon in codons:
        variation: list = [codon]
        for i in range(3):
            bases = ['A', 'U', 'G', 'C']
            bases.remove(codon[i])
            variation.extend([codon[:i] + base + codon[i + 1:] for base in bases])
        variations.append(variation)

    variation_matrix = pd.DataFrame(variations)
    variation_matrix.columns = ['Original'] + [f'Variation {i + 1}' for i in range(9)]
    return variation_matrix


def translate_codons_df(variation_matrix: pd.DataFrame) -> pd.DataFrame:
    """takes in a variation matrix for a codon sequence and translates all codons into AA"""
    translated_df = pd.DataFrame()

    for column in variation_matrix.columns:
        codons = variation_matrix[column]
        amino_acids = [translate_codon_to_aa(codon) for codon in codons]
        translated_df[column] = amino_acids

    return translated_df


def translate_codons_to_string(codons: list) -> str:
    """Takes in a list containing all codons. Can be used to generate AA sequence from DNA sequence found
    online. Newly generated AA sequence and sequence found in "reference_file_substitutions.csv" can be
    compared by BLAST to ensure comparability"""
    amino_acids = []

    for codon in codons:
        amino_acid = translate_codon_to_aa(codon)
        amino_acids.append(amino_acid)

    return ''.join(amino_acids)


def translate_codon_to_aa(codon: str) -> str:
    """used to translate single codons into corresponding AA"""
    codon_table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    if codon in codon_table:
        return codon_table[codon]
    else:
        return 'Unknown'


def prob_aa_position(position: int, variation_matrix: pd.DataFrame, include_original:  bool = False) -> pd.Series:
    """Function returning the probability for each AMS, resulting from a single mutation in the codon as pd.Series"""

    res: pd.Series = variation_matrix.loc[position].value_counts(normalize=True)
    if include_original:
        return res
    else:
        return res.drop(variation_matrix.Original.iloc[position])


def select_smut_position(position: int, variation_matrix: pd.DataFrame, include_original: bool = False) -> pd.Series:
    """Return all aa, resulting from a single mutation. Probability is set to 1 in a Series"""

    res: pd.Series = variation_matrix.loc[position].value_counts(normalize=True)
    res.iloc[:] = 1
    if include_original:
        return res
    else:
        return res.drop(variation_matrix.Original.iloc[position])


def exchange_prob_dict(var_frame: pd.DataFrame, bias_dms: bool = True, include_original: bool = False) -> dict:

    res: dict = {}
    if bias_dms:
        for position in range(var_frame.shape[0]):
            res[position] = prob_aa_position(position, var_frame, include_original)
    else:
        for position in range(var_frame.shape[0]):
            res[position] = select_smut_position(position, var_frame, include_original)

    return res


def clean_variation_matrix(variation_matrix: pd.DataFrame, add_val=np.nan) -> pd.DataFrame:
    """  """
    variation_matrix_cleaned = variation_matrix.replace("*", add_val)

    return variation_matrix_cleaned


def prob_smut(var_mat: pd.DataFrame, dms_scores: pd.DataFrame, bias_dms: bool = True, include_original: bool = False) -> pd.DataFrame:
    """Returns df with probabilities for single mutations"""
    res = pd.DataFrame(columns=dms_scores.columns, index=dms_scores.index, data=np.zeros(dms_scores.shape))
    prob_dict: dict = exchange_prob_dict(var_mat, bias_dms, include_original)

    for position in prob_dict.keys():
        res.loc[position] = res.loc[position].add(prob_dict[position])

    return res


def dms_smut(codon_seq: list, dms_data: pd.DataFrame, bias_dms: bool = True, include_original: bool = False) -> pd.DataFrame:
    """Returns df with probability adjusted dms_scores
        :codon_seq - Sequence of Protein (DNA - Sequence)
        :dms_data - DataFrame from DMS Experiment
    """

    codon_var_raw: pd.DataFrame = translate_codons_df(generate_codon_variations(codon_seq))
    codon_var: pd.DataFrame = clean_variation_matrix(codon_var_raw)
    # all original AA's are assigned the dms score 0
    dms_scores = dc.norm(dc.df_split(dms_data).replace(np.nan, 0))
    prob_single_mut: pd.DataFrame = prob_smut(codon_var, dms_scores, bias_dms, include_original)

    if bias_dms:
        return dc.norm(dms_scores * prob_single_mut)
    else:
        return dms_scores * prob_single_mut


def select_smut(DMS_scores: pd.DataFrame, variation_matrix: pd.DataFrame) -> pd.DataFrame:
    """Takes in a df from the DMS_data dataset, that needs to be processed by dc.df_transform.T (!!!) and a df showing all AA
     variations for each position. Returns a df, that selects only the DMS_scores for AA that are generated by a single
    mutation. Variation Matrix needs to be cleaned with "clean_variation_matrix" before! """
    variation_matrix_clean = clean_variation_matrix(variation_matrix)

    single_mutations = pd.DataFrame()
    for position in range(0, DMS_scores.shape[0]):
        mut_per_pos = prob_aa_position(position, variation_matrix_clean).index
        selected_columns_df = DMS_scores.loc[position+1, mut_per_pos]
        single_mutations = pd.concat([single_mutations, selected_columns_df], ignore_index=True)
    # used pd.concat instead of DataFrame.append, because it won't be supported in future versions
    single_mutations.set_index(DMS_scores.index)
    return single_mutations


# def prob_smut(df_smut: pd.DataFrame, variation_matrix: pd.DataFrame) -> pd.DataFrame:
#     """Takes in a df containing only single mutations and a df with all AA variants of a proteins' AA sequence.
#     Multiplies severity and probability of a DMS_score to get a better feeling for the effects of the mutation"""
#
#     all_probs = pd.DataFrame()
#     for position in range(0, df_smut.shape[0]):
#         probs_per_pos = prob_aa_position(position, variation_matrix)
#         probs = df_smut.loc[position].multiply(probs_per_pos)
#         all_probs = pd.concat([all_probs, probs], axis=1)
#
#     return all_probs.T
#


def generate_codon_variations_rdm(codons: list) -> pd.DataFrame:
    variations = []
    all_codons = [''.join(codon) for codon in itertools.product('ACGT', repeat=3)]
    for codon in codons:
        variation = [codon]
        random_codons = random.sample(all_codons, 9)
        variation.extend(random_codons)
        variations.append(variation)

    variation_matrix = pd.DataFrame(variations)
    variation_matrix.columns = ['Original'] + [f'Variation {i + 1}' for i in range(9)]
    return variation_matrix




if __name__ == '__main__':

    gia_null_eto: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Giacomelli_NULL_Etoposide_2018.csv')

    pass

