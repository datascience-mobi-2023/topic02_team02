import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import data_cleanup as dc
import numpy as np


def hmap(frame: pd.DataFrame) -> None:
    hmap_frame: pd.DataFrame = frame.pivot(index='AS_new', columns=['position_mut', 'AS_old'], values='DMS_score')
    plt.figure(figsize=(50, 8))
    sns.set(font_scale=2)
    sns.heatmap(hmap_frame, cmap='seismic')
    plt.title('DMS Scores for Mutations')
    plt.xlabel('Position and Wildtype AA')
    plt.ylabel('Mutated AA')
    plt.show()
    return None


def mult_hmap(data1: pd.DataFrame, data2: pd.DataFrame, data3: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(30, 20), sharex=True)
    sns.set(font_scale=2)
    frames = [data1, data2, data3]
    for frame in frames:
        dc.aufteilung_mut_pos(frame)
    names = ["Gia_null_eto", "Gia_wt_nut", "Gia_null_nut"]
    for i, ax in enumerate(axes):
        hmap_frame: pd.DataFrame = frames[i].pivot(index='AS_new', columns=['position_mut', 'AS_old'],
                                                   values='DMS_score')
        sns.heatmap(hmap_frame, cmap='seismic', ax=ax)
        ax.set_title('Custom Title {}'.format(i+1))
        if i < 3:
            ax.set(ylabel='Mutated AA')
        ax.set_title('')
        if i < 2:
            ax.set(xlabel='')
        ax.text(0.05, 0.95, names[i], transform=ax.transAxes, fontsize=24,
                verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    plt.tight_layout()
    plt.xlabel('Position and Wildtype AA')
    plt.show()


def linegraph(dataframes):
    plt.figure(figsize=(35, 7))
    for df in dataframes:
        x = df.columns.get_level_values('position_mut').astype(int)
        y = df.iloc[0].values.astype(float)
        label = df.name
        plt.plot(x, y, label=label, marker='o')
    plt.xlabel('position')
    plt.ylabel('DMS_score')
    plt.title('Mean DMS_scores of different datasets throughout each position of tp53')
    plt.xticks(np.arange(0, 393+1, 10))
    plt.legend()
    plt.show()


def multiple_linegraph(*dataframes):
    plt.figure(figsize=(50, 8))

    for df in dataframes:
        df_transposed = df.transpose()
        df_transposed.plot.line()

    plt.xlabel('position_mut')
    plt.ylabel('DMS_score')
    plt.legend(dataframes)
    return plt.show()


def hmap_mean_variance(df: pd.DataFrame) -> None:
    plt.figure(figsize=(50, 10))
    sns.heatmap(df, cmap='coolwarm', annot=True, fmt=".2f", linewidths=0.5)
    plt.title('Mean substitution Heatmap of Mutated AA vs Wildtype AA')
    plt.xlabel('Mutated AA')
    plt.ylabel('Wildtype AA')
    plt.show()
    return None
