import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import data_cleanup as dc
import Documentation as doc

gia_null_eto: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Giacomelli_NULL_Etoposide_2018.csv')
gia_null_nut: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Giacomelli_NULL_Nutlin_2018.csv')
gia_wt_nut: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Giacomelli_WT_Nutlin_2018.csv')
kot_hum: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Kotler_2018.csv')

anzahl: int = 10
def LLV():
    GNELLV: pd.DataFrame = dc.low_val(gia_null_eto, anzahl);
    GNESLV = dc.low_val(gia_null_eto, anzahl);
    GNNLLV = dc.low_val(gia_null_nut, anzahl);
    GNNSLV = dc.low_val(gia_null_nut, anzahl);
    GWNLLV = dc.low_val(gia_wt_nut, anzahl);
    GWNSLV = dc.low_val(gia_wt_nut, anzahl);
    KLLV = dc.low_val(kot_hum, anzahl);
    KSLV = dc.low_val(kot_hum, anzahl);
    return None