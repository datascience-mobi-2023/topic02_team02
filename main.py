import pandas as pd

import data_cleanup as dc
import functions as func


if __name__ == '__main__':
    fpath = './DMS_data/P53_HUMAN_Giacomelli_NULL_Etoposide_2018.csv'
    df = pd.read_csv(fpath)
    datensatz_normalisiert = dc.set_transform_norm(df)
    dat_norm_aufgeteilt: pd.DataFrame = dc.aufteilung_mut_pos(datensatz_normalisiert)
    func.hmap(dat_norm_aufgeteilt)