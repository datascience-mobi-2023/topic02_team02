## all added functions are to be declared in the init file, to ensure effortless usage
import pandas as pd
import data_cleanup as dc


def slice_domain(DMS_data: pd.DataFrame, start: str, end: str) -> pd.DataFrame:
    """ Takes in a Data set in the form of the DataFrames in the DMS_data folder. Also takes in Positions to bel sliced
    as start and end position. """
    mutations_df = dc.aufteilung_mut_pos(DMS_data)
    sliced_domains = mutations_df.loc[(mutations_df['position_mut'] >= start) & (mutations_df["position_mut"] <= end)]
    return sliced_domains
