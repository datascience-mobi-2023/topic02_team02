## all added functions are to be declared in the init file, to ensure effortless usage
import pandas as pd

def min_max_norm(frame: pd.DataFrame, upper_border: float = 1.0, lower_border: float = -1.0) -> pd.DataFrame:
    """Der gegebene Datensatz wird normalisiert auf die Grenzen upper und lower, die per default auf 1, -1 stehen"""
    NormalisierungsDatensatz: pd.DataFrame = frame.copy()
    val_min: float = NormalisierungsDatensatz.DMS_score.min()
    val_max: float = NormalisierungsDatensatz.DMS_score.max()
    counter: int = 1

    while counter <= NormalisierungsDatensatz.shape[0] - 1:
        NormalisierungsDatensatz.iloc[counter, 2] = (NormalisierungsDatensatz.iloc[counter, 2] - val_min)/(val_max - val_min) * (upper_border - lower_border) -1
        counter += 1
    return NormalisierungsDatensatz

def set_transform_norm(frame: pd.DataFrame) -> pd.DataFrame:
    """Set transformation operation"""
    mean_val: float = frame.DMS_score.mean()
    var_val: float = frame.DMS_score.std()
    frame.DMS_score = (frame.DMS_score - mean_val) / var_val
    return frame

def aufteilung_mut_pos(frame: pd.DataFrame) -> pd.DataFrame:
    """Mutation, Position der Mutation und neue AS seden in drei separate Spalten des Frames aufgeteilt"""
    frame['position_mut'] = frame.mutant.str.slice(start=1, stop=-1).astype(int)
    frame['AS_old'] = frame.mutant.str.get(0)
    frame['AS_new'] = frame.mutant.str.get(-1)
    return frame

def df_transform(frame: pd.DataFrame) -> pd.DataFrame:
    """df wird in Format WT_Sequenz-Mutierte_AS umgewandelt"""
    frame['position_mut'] = frame.mutant.str.slice(start=1, stop=-1).astype(int)
    frame['AS_old'] = frame.mutant.str.get(0)
    frame['AS_new'] = frame.mutant.str.get(-1)
    pd.DataFrame = frame.pivot(index='AS_new', columns=['position_mut', 'AS_old'], values='DMS_score')
    return frame