import pandas as pd
import data_cleanup as dc
import data_exploration as de
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics import silhouette_score
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.pairwise import euclidean_distances
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.decomposition import PCA
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist


def mean_substitutions(frame: pd.DataFrame) -> pd.DataFrame:
    """calculates the mean substitution values for each substitution directly from a DMS_data set."""
    # Check for dependencies
    if "AS_new" not in frame.columns:
        raise ValueError(f"Die Funktion 'aufteilung_mut_pos()' muss vorher ausgeführt werden.")

    subs_df = frame.groupby(["AS_old", "AS_new"])
    mean_scores = subs_df.DMS_score.mean()
    mean_scores_df = mean_scores.reset_index()
    mean_scores_df = mean_scores_df.pivot(index="AS_old", columns="AS_new", values="DMS_score")
    return mean_scores_df


def direct_mean_subs(fpath: str):
    """ does the same as mean_substitutions above, but with a more direct approach. Input is the path of a DMS_data
    dataset, which is cleaned. Output is cleaned better than the function above """
    df = pd.read_csv(fpath)
    mutations_df = dc.aufteilung_mut_pos(df)
    mutations_df_norm: pd.DataFrame = dc.norm(mutations_df)
    subs_df = mutations_df_norm.groupby(["AS_old", "AS_new"])
    mean_scores = subs_df.DMS_score.mean()
    mean_scores_df = mean_scores.reset_index()
    mean_subs = mean_scores_df.pivot(index="AS_old", columns="AS_new", values="DMS_score")
    return dc.rmv_na(mean_subs)


def mean_substitutions_inverted(frame: pd.DataFrame) -> pd.DataFrame:
    """calculate the inverted mean_substitutions matrix. Outdated and unnecessary, as instead mean_substitutions().T can
    be used."""

    if "AS_new" and "AS_old" not in frame.columns:
        raise ValueError(f"Die Funktion 'aufteilung_mut_pos()' muss vorher ausgeführt werden.")

    subs_dfi = frame.groupby(["AS_old", "AS_new"])
    mean_scores = subs_dfi.DMS_score.mean()
    mean_scores_dfi = mean_scores.reset_index()
    mean_scores_dfi = mean_scores_dfi.pivot(index="AS_new", columns="AS_old", values="DMS_score")
    return mean_scores_dfi


def aa_distance_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    """calculate the non-symmetric distance matrix out of the Dataset for the chemical properties of AAs"""
    aa_nat = frame.drop(index=[12, 18])
    labels_column = 'Letter'
    aa_rmv = aa_nat.drop(['Name', 'Abbr', 'Letter', 'Molecular Formula', 'Molecular Weight', 'Residue Formula', 'pKx3'],
                         axis=1)
    normed_df = dc.norm(aa_rmv)
    aa_distances = euclidean_distances(normed_df.values)
    frame = pd.DataFrame(aa_distances, index=aa_nat[labels_column], columns=aa_nat[labels_column])
    return frame


def dms_distance_matrix_wt(frame: pd.DataFrame) -> pd.DataFrame:
    """calculates the distances of the wild-type AAs to each other based on the DMS-Scores when interchanged with
    another AA"""
    frame_prep = dc.rmv_na(de.mean_substitutions(frame))
    dms_distances = euclidean_distances(frame_prep.values)
    frame_prep = pd.DataFrame(dms_distances, index=frame_prep.index, columns=frame_prep.index)
    return frame_prep


def dms_distance_matrix_mutated(frame: pd.DataFrame) -> pd.DataFrame:
    """calculates the distances of the mutated AA to each other based on the DMS-Scores when interchanged with
    another AA"""
    frame_prep = dc.rmv_na(de.mean_substitutions(frame).T)
    dms_distances = euclidean_distances(frame_prep.values)
    frame_prep = pd.DataFrame(dms_distances, index=frame_prep.index, columns=frame_prep.index)
    return frame_prep


def plot_hier_clust(distance_matrix: pd.DataFrame, title: str):
    """Takes in a non-symmetric distance matrix and condenses it. Hierarchical clustering is performed and the
    dendrogram is plotted."""
    condensed_dist = pdist(distance_matrix)
    link = hierarchy.linkage(condensed_dist, method='ward')
    plt.figure(figsize=(10, 5))
    hierarchy.dendrogram(link, labels=distance_matrix.index)
    plt.title(f'Dendrogram of hierarchical clustering based on {title}')
    plt.xlabel('Amino acids')
    plt.ylabel('Distance')
    return plt.show()


def determine_clusters_silhouette(dms_data: pd.DataFrame, min_clusters=2, max_clusters=10) -> int:
    """determine optimal number of clusters using hierarchical clustering for DMS_data, that was processed like this:
    dc.rmv_na(dc.df_transform(DMS_data)). Not perfectly optimized, implemented to proof a concept"""
    best_score = -1
    best_clusters = 0

    pca = PCA(n_components=2)
    pca_data = pca.fit_transform(dms_data)

    for num_clusters in range(min_clusters, max_clusters + 1):
        hc = AgglomerativeClustering(n_clusters=num_clusters)
        cluster_labels = hc.fit_predict(pca_data)

        score = silhouette_score(dms_data, cluster_labels)

        if score > best_score:
            best_score = score
            best_clusters = num_clusters

    return best_clusters


def pca_hierarchical_plot(dist_matrix: pd.DataFrame, optimal_num_cluster: int, title: str, show_var=False):
    """used to perform a pca on the provided distance matrix. Then, a hierarchical clustering on the pca_results is
    done. Output is a plot in which each datapoint is one AA, shown in the color of the cluster it belongs to. Showing
    explained variance by first two PCs is optional."""

    pca = PCA(n_components=2)
    pca_data = pca.fit_transform(dist_matrix)

    # Access explained variance ratio
    explained_variance = pca.explained_variance_ratio_

    hc = linkage(pca_data, method='ward')
    num_clusters = optimal_num_cluster
    cluster_labels = fcluster(hc, num_clusters, criterion='maxclust')

    pca_df = pd.DataFrame(data=pca_data, columns=['PC1', 'PC2'])
    pca_df['Cluster'] = cluster_labels.astype(str)
    pca_df['Label'] = dist_matrix.index

    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Cluster', palette='Set1').set(title=f'PCA of {title}')

    for i, point in pca_df.iterrows():
        plt.annotate(point['Label'], (point['PC1'], point['PC2']), textcoords="offset points", xytext=(0, -10),
                     ha='center')

    plt.legend(title='Cluster', loc='lower right')

    if show_var is True:
        print(f"Explained variance for PCA of {title}:")
        for i, ratio in enumerate(explained_variance):
            print(f"PC{i + 1}: {ratio:.2f}")

    return plt.show()
