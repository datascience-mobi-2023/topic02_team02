## all added functions are to be declared in the init file, to ensure effortles usage
import pandas as pd
def normalisierung(frame: pd.DataFrame, upper_border: float = 1.0, lower_border: float = -1.0) -> pd.DataFrame:
    """Der gegebene Datensatz wird normalisiert auf die Grenzen upper und lower, die per default auf 1, -1 stehen"""
    NormalisierungsDatensatz: pd.DataFrame = frame.copy()
    val_min: float = NormalisierungsDatensatz.DMS_score.min()
    val_max: float = NormalisierungsDatensatz.DMS_score.max()
    counter: int = 1

    while counter <= NormalisierungsDatensatz.shape[0] - 1:
        NormalisierungsDatensatz.iloc[counter, 2] = (NormalisierungsDatensatz.iloc[counter, 2] - val_min)/(val_max - val_min) * (upper_border - lower_border) -1
        counter += 1
    return NormalisierungsDatensatz