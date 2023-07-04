import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import data_cleanup as dc



def load_data_frame(file_path: list) -> pd.DataFrame:
    'Load list of .csv files into one data Frame'
    # list comprehensions are better practice
    container: list = [pd.read_csv(filename, index_col=None, header=0) for filename in file_path]

    # transformation into data frame (dict)
    frame = pd.concat(container, axis=0, ignore_index=True)
    print("The frame has been successfully loaded")
    return frame


def hmap(frame: pd.DataFrame) -> None:
    hmap_frame: pd.DataFrame = frame.pivot(index='AS_new', columns=['position_mut', 'AS_old'], values='DMS_score')
    plt.figure(figsize=(50, 8))
    sns.set(font_scale=2)
    sns.heatmap(hmap_frame, cmap='seismic')
    plt.title('DMS Scores for Mutations')
    plt.show()
    return None


def mult_hmap(Daten1: pd.DataFrame, Daten2: pd.DataFrame, Daten3: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(40, 29), sharex=True)
    sns.set(font_scale=2)
    frames = [Daten1, Daten2, Daten3]
    for frame in frames:
        frame = dc.aufteilung_mut_pos(frame)
    names = ["Gia_null_eto", "Gia_wt_nut", "Gia_null_nut"]
    for i, ax in enumerate(axes):
        hmap_frame: pd.DataFrame = frames[i].pivot(index='AS_new', columns=['position_mut', 'AS_old'], values='DMS_score')
        sns.heatmap(hmap_frame, cmap='seismic', ax=ax)
        ax.set_title('')
        if i < 2:
            ax.set(xlabel='')
        ax.text(0.05, 0.95, names[i], transform=ax.transAxes, fontsize=24,
                verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    plt.tight_layout()
    plt.show()

def df_mean(df: pd.DataFrame) -> pd.DataFrame:
    df_trafo: pd.DataFrame = dc.df_transform(df)
    df_trafo_narmv: pd.DataFrame = dc.rmv_na(df_trafo)
    df_trafo_narmv_mean = pd.DataFrame(columns=df_trafo_narmv.columns)

    for column in df_trafo_narmv.columns:
        if column != 'position_mut' and column != 'AS_old':
            column_mean = df_trafo_narmv[column].mean()
            df_trafo_narmv_mean.loc[0, column] = column_mean

    df_trafo_narmv_mean = df_trafo_narmv_mean.reset_index(drop=True)

    return df_trafo_narmv_mean

def df_mean(df_trafo: pd.DataFrame) -> pd.DataFrame:
    df_trafo_narmv: pd.DataFrame = dc.rmv_na(df_trafo)
    df_trafo_narmv_mean = pd.DataFrame(columns=df_trafo_narmv.columns)

    for column in df_trafo_narmv.columns:
        if column != 'position_mut' and column != 'AS_old':
            column_mean = df_trafo_narmv[column].mean()
            df_trafo_narmv_mean.loc[0, column] = column_mean

    df_trafo_narmv_mean = df_trafo_narmv_mean.reset_index(drop=True)

    return df_trafo_narmv_mean



def multiple_linegraph(*dataframes):

    plt.figure(figsize=(50, 8))

    for df in dataframes:
        df_transposed = df.transpose()
        df_transposed.plot.line()

    plt.xlabel('position_mut')
    plt.ylabel('DMS_score')
    plt.legend(dataframes)
    plt.show()

def calculate_average_dms_score_new(*args):
    results = {}

    for arg in args:
        df_name = arg[0]
        df = arg[1]
        grouped = df.groupby('AS_new')
        sums = grouped['DMS_score'].sum()
        counts = grouped['DMS_score'].count()
        averages = sums / counts
        results[df_name] = averages

    result_df = pd.DataFrame(results)
    return result_df

def calculate_average_dms_score_old(*args):
    results = {}

    for arg in args:
        df_name = arg[0]
        df = arg[1]
        grouped = df.groupby('AS_old')
        sums = grouped['DMS_score'].sum()
        counts = grouped['DMS_score'].count()
        averages = sums / counts
        results[df_name] = averages

    result_df = pd.DataFrame(results)
    return result_df

def hmap_mean_variance (df = pd.DataFrame) -> None:
    plt.figure(figsize=(50, 10))
    sns.heatmap(df, cmap='coolwarm', annot=True, fmt=".2f", linewidths=0.5)
    plt.title('Heatmap of AS_new vs AS_old')
    plt.xlabel('AS_new')
    plt.ylabel('AS_old')
    plt.show()
    return None