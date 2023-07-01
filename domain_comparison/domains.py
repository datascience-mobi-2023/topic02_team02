## all added functions are to be declared in the init file, to ensure effortless usage
import pandas as pd
import pip
from scipy.stats import shapiro
import data_cleanup as dc

def slice_domain(DMS_data: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    """Takes in a Data set in the form of the DataFrames in the DMS_data folder.
    Also takes in positions to be sliced as start and end position."""
    mutations_df_copy = DMS_data.copy()
    mutations_df = dc.aufteilung_mut_pos(mutations_df_copy)
    sliced_domains = mutations_df.loc[(mutations_df['position_mut'] >= start) & (mutations_df['position_mut'] <= end)].copy()

    return sliced_domains


def test_normality(data, alpha=0.05):
    """
    Performs a Shapiro-Wilk test for normality on a list of values.
    """
    statistic, p_value = shapiro(data)
    is_normal = p_value > alpha
    results = {
        'statistic': statistic,
        'p-value': p_value,
        'is_normal': is_normal
    }
    return results
