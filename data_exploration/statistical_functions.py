import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances
import scipy.stats as stats
import data_cleanup as dc
import data_exploration as de
from sklearn.metrics import silhouette_score
from sklearn.cluster import AgglomerativeClustering


def mean_substitutions(frame: pd.DataFrame) -> pd.DataFrame:
    """aus einem Datensatz, den wir gegeben haben, direkt eine Matrix mit Mittelwerten des Austausches zu erstellen. Zeilen in dem ausgegebenen
    DataFrame sind die alten AS, Spalten die neuen"""
    # Check for dependencies
    if "AS_new" not in frame.columns:
        raise ValueError(f"Die Funktion 'aufteilung_mut_pos()' muss vorher ausgeführt werden.")

    subs_df = frame.groupby(["AS_old", "AS_new"])
    mean_scores = subs_df.DMS_score.mean()
    mean_scores_df = mean_scores.reset_index()
    mean_scores_df = mean_scores_df.pivot(index="AS_old", columns="AS_new", values="DMS_score")
    return mean_scores_df


def mean_substitutions_inverted(frame: pd.DataFrame) -> pd.DataFrame:
    """aus einem Datensatz, den wir gegeben haben, direkt eine Distanzmatrix zu erstellen. Zeilen in dem ausgegebenen
    DataFrame sind die neuen AS, Spalten die alten"""

    if "AS_new" and "AS_old" not in frame.columns:
        raise ValueError(f"Die Funktion 'aufteilung_mut_pos()' muss vorher ausgeführt werden.")

    subs_dfi = frame.groupby(["AS_old", "AS_new"])
    mean_scores = subs_dfi.DMS_score.mean()
    mean_scores_dfi = mean_scores.reset_index()
    mean_scores_dfi = mean_scores_dfi.pivot(index="AS_new", columns="AS_old", values="DMS_score")
    return mean_scores_dfi


def aa_distance_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    """  """
    aa_nat = frame.drop(index=[12, 18])
    labels_column = 'Letter'
    aa_rmv = aa_nat.drop(['Name', 'Abbr', 'Letter', 'Molecular Formula', 'Molecular Weight', 'Residue Formula', 'pKx3'], axis=1)
    aa_zscore = dc.min_max_norm(aa_rmv.apply(stats.zscore))
    aa_distances = euclidean_distances(aa_zscore.values)
    frame = pd.DataFrame(aa_distances, index=aa_nat[labels_column], columns=aa_nat[labels_column])
    return frame


def dms_distance_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    """calculates the distances of the AA to each other based on the DMS-Scores when interchanged with another AA"""
    frame_prep = dc.rmv_na(de.mean_substitutions(frame))
    dms_distances = euclidean_distances(frame_prep.values)
    frame_prep = pd.DataFrame(dms_distances, index=frame_prep.index, columns=frame_prep.index)
    return frame_prep


def plot_dendrogram(model, **kwargs):
    from scipy.cluster.hierarchy import dendrogram
    counts = pd.Series(model.children_[:, 1])
    linkage_matrix = pd.DataFrame(model.children_, columns=['cluster_1', 'cluster_2'])
    linkage_matrix['distance'] = model.distances_
    linkage_matrix['new_count'] = counts
    dendrogram(linkage_matrix.to_numpy(), **kwargs)


def determine_clusters_silhouette(dist_matrix, min_clusters=2, max_clusters=10):
    best_score = -1
    best_clusters = 0

    for num_clusters in range(min_clusters, max_clusters + 1):
        hc = AgglomerativeClustering(n_clusters=num_clusters)
        cluster_labels = hc.fit_predict(dist_matrix)

        score = silhouette_score(dist_matrix, cluster_labels)

        if score > best_score:
            best_score = score
            best_clusters = num_clusters

    return best_clusters
