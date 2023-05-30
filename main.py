import pandas as pd

import data_cleanup as dc
import functions as func

if __name__ == '__main__':
    fpath = ['./DMS_data/P53_HUMAN_Giacomelli_NULL_Etoposide_2018.csv',
             './DMS_data/P53_HUMAN_Giacomelli_WT_Nutlin_2018.csv',
             './DMS_data/P53_HUMAN_Giacomelli_NULL_Nutlin_2018.csv']

    for data_set in fpath:
        df = pd.read_csv(data_set)
        datensatz_normalisiert = dc.min_max_norm(dc.set_transform_norm(df))
        dat_norm_aufgeteilt: pd.DataFrame = dc.aufteilung_mut_pos(datensatz_normalisiert)
        func.hmap(dat_norm_aufgeteilt)
