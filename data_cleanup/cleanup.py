import pandas as pd
import data_cleanup as dc
import scipy.stats as stats


def isfloat(series: pd.Series):
    for pos in range(series.shape[0]):
        try:
            float(series.iloc[pos])
            return True
        except ValueError:
            return False


def min_max_norm(norm_df: pd.DataFrame, upper_border: float = 1.0, lower_border: float = -1.0) -> pd.DataFrame:
    """Der gegebene Datensatz wird normalisiert auf die Grenzen upper und lower, die per default auf 1, -1 stehen"""

    for col in norm_df.columns:
        if isfloat(norm_df[col]):
            val_min: float = norm_df[col].min()
            val_max: float = norm_df[col].max()
            norm_df[col] = (norm_df[col] - val_min) / (val_max - val_min) * (upper_border - lower_border) - 1
        else:
            # print(f'{col} contains values, which arent floats')
            continue

    return norm_df


def z_transform(frame: pd.DataFrame) -> pd.DataFrame:
    """Set transformation operation"""
    for col in frame:
        if isfloat(frame[col]):
            mean_val: float = frame[col].mean()
            var_val: float = frame[col].std()
            frame[col] = (frame[col] - mean_val) / var_val
        else:
            continue
    return frame


def norm(frame: pd.DataFrame) -> pd.DataFrame:
    """Combining z-transform and min_max_norm for easy usage"""
    res: pd.DataFrame = min_max_norm(z_transform(frame))

    return res


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
    pivoted_frame = frame.pivot(index='AS_new', columns=['position_mut', 'AS_old'], values='DMS_score')
    return pivoted_frame


def df_split(frame: pd.DataFrame) -> pd.DataFrame:
    frame['position_mut'] = frame.mutant.str.slice(start=1, stop=-1).astype(int)-1
    frame['aa_new'] = frame.mutant.str.get(-1)
    res = frame.pivot(index='position_mut', values='DMS_score', columns='aa_new')
    return res


def rmv_na(frame: pd.DataFrame) -> pd.DataFrame:
    """Ändert die Werte eines df von NaN zu 0, wenn mutierte AS der WT AS entspricht"""
    for col in frame.columns:
        for row in frame.index:
            col_str = str(col)
            row_str = str(row)
            letters_row = ''.join(filter(str.isalpha, row_str))
            letters_col = ''.join(filter(str.isalpha, col_str))

            if letters_row == letters_col and pd.isna(frame.loc[row, col]):
                frame.loc[row, col] = 0

    return frame


def min_max_val(frame: pd.DataFrame) -> dict:
    """returns lowes and highest value of frame with position"""
    m_val: pd.Series = frame['DMS_score'].idxmin()
    h_val: pd.Series = frame['DMS_score'].idxmax()

    return {"min": m_val, "max": h_val}


def high_val(df: pd.DataFrame, num_high: int, ) -> pd.DataFrame:
    """Zeigt die x höchsten Werte eines Datensatzes an"""
    gia_null_eto_dft = df_transform(df)
    gia_null_eto_dft_narm = rmv_na(gia_null_eto_dft)
    sum_df = pd.DataFrame(gia_null_eto_dft_narm.sum(), columns=['Sum'])
    highest_values = sum_df.nlargest(num_high, "Sum")
    # print(f"\nHighest {num_high} values in Sum :")
    # print(highest_values)

    return highest_values


def low_val(df: pd.DataFrame, num_low: int) -> pd.DataFrame:
    """Zeigt die x niedrigsten Werte eines Datensatzes an"""
    gia_null_eto_dft = df_transform(df)
    gia_null_eto_dft_narm = rmv_na(gia_null_eto_dft)
    sum_df = pd.DataFrame(gia_null_eto_dft_narm.sum(), columns=['Sum'])
    lowest_values = sum_df.nsmallest(num_low, "Sum")
    # print(f"\nLowest {num_low} values in Sum:")
    # print(lowest_values)

    return lowest_values


def clean_aa(frame: pd.DataFrame) -> pd.DataFrame:
    """calculate the cleaned Dataset for the chemical properties of AAs"""
    aa_nat = frame.drop(index=[12, 18])
    labels_column = 'Letter'
    aa_rmv = aa_nat.drop(['Name', 'Abbr', 'Letter', 'Molecular Formula', 'Molecular Weight', 'Residue Formula', 'pKx3'],
                         axis=1)
    aa_zscore = dc.min_max_norm(aa_rmv.apply(stats.zscore))
    frame = pd.DataFrame(aa_zscore.values, index=aa_nat[labels_column], columns=aa_rmv.columns)
    return frame


def load_data_frame(file_path: list) -> pd.DataFrame:
    """Load list of .csv files into one data Frame"""
    # list comprehensions are better practice
    container: list = [pd.read_csv(filename, index_col=None, header=0) for filename in file_path]

    # transformation into data frame (dict)
    frame = pd.concat(container, axis=0, ignore_index=True)
    print("The frame has been successfully loaded")
    return frame


if __name__ == "__main__":
    df: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Giacomelli_WT_Nutlin_2018.csv')
    print(min_max_norm(df))
    # functions.hmap(aufteilung_mut_pos(min_max_norm(df)))
