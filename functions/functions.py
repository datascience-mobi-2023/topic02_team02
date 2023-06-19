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


def multiple_hmap(*frames: pd.DataFrame) -> None:
    hmap_frames = []
    for frame in frames:
        hmap_frame = frame.pivot(index='AS_new', columns=['position_mut', 'AS_old'], values='DMS_score')
        hmap_frames.append(hmap_frame)

    plt.figure(figsize=(50, 8 * len(frames)))
    sns.set(font_scale=2)

    for i, hmap_frame in enumerate(hmap_frames):
        plt.subplot(len(frames), 1, i + 1)
        sns.heatmap(hmap_frame, cmap='seismic')
        plt.title(f'DMS Scores for Mutations - Dataset {i + 1}')

    plt.tight_layout()
    plt.show()
    return None
