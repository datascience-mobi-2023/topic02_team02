import pandas as pd
import data_cleanup as dc
from sklearn.metrics.pairwise import euclidean_distances
import scipy.stats as stats
import data_exploration as de

def mean_substitutions(frame: pd.DataFrame) -> pd.DataFrame:
    """aus einem Datensatz, den wir gegeben haben, direkt eine Distanzmatrix zu erstellen. Zeilen in dem ausgegebenen
    DataFrame sind die alten AS, Spalten die neuen"""
    mutations_df = dc.aufteilung_mut_pos(frame)
    subs_df = mutations_df.groupby(["AS_old", "AS_new"])
    mean_scores = subs_df.DMS_score.mean()
    mean_scores_df = mean_scores.reset_index()
    frame = mean_scores_df.pivot(index="AS_old", columns="AS_new", values="DMS_score")
    return frame

def mean_substitutions_inverted(frame: pd.DataFrame) -> pd.DataFrame:
    """aus einem Datensatz, den wir gegeben haben, direkt eine Distanzmatrix zu erstellen. Zeilen in dem ausgegebenen
    DataFrame sind die neuen AS, Spalten die alten"""
    mutations_dfi = dc.aufteilung_mut_pos(frame)
    subs_dfi = mutations_dfi.groupby(["AS_old", "AS_new"])
    mean_scores = subs_dfi.DMS_score.mean()
    mean_scores_dfi = mean_scores.reset_index()
    frame = mean_scores_dfi.pivot(index="AS_new", columns="AS_old", values="DMS_score")
    return frame

def aa_distance_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    """  """
    aa_nat = frame.drop(index=[12, 18])
    labels_column = 'Letter'
    aa_rmv = aa_nat.drop(['Name', 'Abbr', 'Letter', 'Molecular Formula', 'Molecular Weight', 'Residue Formula', 'pKx3'], axis=1)
    aa_zscore = aa_rmv.apply(stats.zscore)
    aa_distances = euclidean_distances(aa_zscore.values)
    frame = pd.DataFrame(aa_distances, index=aa_nat[labels_column], columns=aa_nat[labels_column])
    return frame

