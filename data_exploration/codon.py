import pandas as pd
from Bio.Seq import Seq

dna_sequence: str = 'ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGAC'
rna_sequence = dna_sequence.replace("T", "U")
p53_codons = [rna_sequence[i:i+3] for i in range(0, len(rna_sequence), 3)]


def translate_codons_df(df: pd.DataFrame) -> pd.DataFrame:
    translated_df = pd.DataFrame()

    for column in df.columns:
        codons = df[column]
        seqs = [Seq(codon) for codon in codons]
        amino_acids = [seq.translate() for seq in seqs]
        translated_df[column] = amino_acids

    return translated_df.astype(str)


def generate_codon_variations(codons: list) -> pd.DataFrame:
    variations: list = []
    for codon in codons:
        variation: list = [codon]
        for i in range(3):
            bases = ['A', 'U', 'G', 'C']
            bases.remove(codon[i])
            variation.extend([codon[:i] + base + codon[i+1:] for base in bases])
        variations.append(variation)

    df = pd.DataFrame(variations)
    df.columns = ['Original'] + [f'Variation {i+1}' for i in range(9)]
    return df


def prob_as_position(position: int, variation_matrix: pd.DataFrame) -> pd.Series:
    """Function returning the probability for each AMS, resulting from a single mutation in the codon as pd.Series"""

    res: pd.Series = variation_matrix.loc[position].value_counts(normalize=True)

    return res


def exchange_prob_dict(frame: pd.DataFrame) -> dict:
    res: dict = {}
    for position in range(0, frame.shape[0]-1):
        res[position] = prob_as_position(position, frame)

    return res


if __name__ == '__main__':
    p53_var_frame: pd.DataFrame = translate_codons_df(generate_codon_variations(p53_codons))
    temp = exchange_prob_dict(p53_var_frame)

    print(temp[0])
    pass
