import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


### this file is for general functions, which are not related to a specific part of the project.

def fib(x: int) -> int:
    # primitive implementation der Fibonacci Folge -> einfache Rekursion
    if x == 1:
        return 0
    elif x <= 1:
        return 1
    else:
        return x + fib(x-1)

def nums(n):
    yield n
    yield from nums(n+1)

s = nums(2)
def sieb(s):
    n = next(s)
    yield n
    yield from sieb(i for i in s if i%n!=0)

# generiert die ersten x Primzahlen durch das Sieb des Eratostenes.
def primzahlen(x: int) -> list:
    res = []
    prime = sieb(nums(2))
    for i in range(x):
        res.append(next(prime))
    return res

def load_data_frame(file_path: list) -> pd.DataFrame:
    'Load list of .csv files into one data Frame'
    # list comprehensions are better practice
    container: list = [pd.read_csv(filename, index_col=None, header=0) for filename in file_path]

    # transformation into data frame (dict)
    frame = pd.concat(container, axis=0, ignore_index=True)
    print("The frame has been successfully loaded")
    return frame

def hmap(frame: pd.DataFrame):
    hmap_frame: pd.DataFrame = frame.pivot(index='AS_new', columns=['position_mut', 'AS_old'], values='DMS_score')
    plt.figure(figsize=(50, 8))
    sns.set(font_scale=2)
    sns.heatmap(hmap_frame, cmap='seismic')
    plt.title('DMS Scores for Mutations')
    plt.show()
    pass