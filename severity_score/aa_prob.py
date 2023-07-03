import pandas as pd
import numpy as np
import data_cleanup as dc
from Bio.Seq import Seq
import severity_score as ses

dna_sequence: str = 'ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGAC'
rna_sequence = dna_sequence.replace("T", "U")
p53_codons_kotler = [rna_sequence[i:i + 3] for i in range(0, len(rna_sequence), 3)]
p53_codons_gia = p53_codons_kotler
p53_codons_gia[71] = "CGC"


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


def prob_aa_position(position: int, variation_matrix: pd.DataFrame) -> pd.Series:
    """Function returning the probability for each AMS, resulting from a single mutation in the codon as pd.Series"""

    res: pd.Series = variation_matrix.loc[position].value_counts(normalize=True)

    return res.drop(variation_matrix.Original.iloc[position])


def clean_variation_matrix(variation_matrix: pd.DataFrame) -> pd.DataFrame:
    """  """
    variation_matrix_cleaned = variation_matrix.replace("*", np.nan)

    return variation_matrix_cleaned


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


def dms_smut(codon_seq: list, dms_data: pd.DataFrame) -> pd.DataFrame:
    """Returns df with probability adjusted dms_scores"""
    sel_mut_dms_frame: pd.DataFrame = dc.norm(dc.df_transform(dms_data).T)

    codon_var_raw: pd.DataFrame = translate_codons_df(generate_codon_variations(codon_seq))
    codon_var: pd.DataFrame = clean_variation_matrix(codon_var_raw)

    probable_mutations: pd.DataFrame = select_smut(sel_mut_dms_frame, codon_var).sort_index(axis=1)
    dms_scores = dc.norm(dc.df_split(dms_data)).reindex_like(probable_mutations)

    res = dms_scores * probable_mutations
    return res

def prob_smut(df_smut: pd.DataFrame, variation_matrix: pd.DataFrame) -> pd.DataFrame:
    """Takes in a df containing only single mutations and a df with all AA variants of a proteins' AA sequence.
    Multiplies severity and probability of a DMS_score to get a better feeling for the effects of the mutation"""

    all_probs = pd.DataFrame()
    for position in range(0, df_smut.shape[0]):
        probs_per_pos = prob_aa_position(position, variation_matrix)
        probs = df_smut.loc[position].multiply(probs_per_pos)
        all_probs = pd.concat([all_probs, probs], axis=1)

    return all_probs.T




if __name__ == '__main__':

    gia_null_eto: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Giacomelli_NULL_Etoposide_2018.csv')
    res = dms_smut(p53_codons_gia, gia_null_eto)
    pass

