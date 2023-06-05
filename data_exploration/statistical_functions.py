import pandas as pd
import data_cleanup as dc

def p_cor():
    # TODO implement correlation
    pass

def s_cor():
    # TODO implement correlation
    pass
def distance_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    """aus einem Datensatz, den wir gegeben haben, direkt eine Distanzmatrix zu erstellen. Zeilen in dem ausgegebenen
    DataFrame sind die alten AS, Spalten die neuen"""
    mutations_df = dc.aufteilung_mut_pos(frame)
    subs_df = mutations_df.groupby(["AS_old", "AS_new"])
    mean_scores = subs_df.DMS_score.mean()
    mean_scores_df = mean_scores.reset_index()
    frame = mean_scores_df.pivot(index="AS_old", columns="AS_new", values="DMS_score")
    return frame

def distance_matrix_inverted(frame: pd.DataFrame) -> pd.DataFrame:
    """aus einem Datensatz, den wir gegeben haben, direkt eine Distanzmatrix zu erstellen. Zeilen in dem ausgegebenen
    DataFrame sind die alten AS, Spalten die neuen"""
    mutations_dfi = dc.aufteilung_mut_pos(frame)
    subs_dfi = mutations_dfi.groupby(["AS_old", "AS_new"])
    mean_scores = subs_dfi.DMS_score.mean()
    mean_scores_dfi = mean_scores.reset_index()
    frame = mean_scores_dfi.pivot(index="AS_new", columns="AS_old", values="DMS_score")
    return frame

