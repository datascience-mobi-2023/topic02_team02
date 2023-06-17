
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import data_cleanup as dc

###################
#Darios Daten zur Visualisierung:
gia_null_eto: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Giacomelli_NULL_Etoposide_2018.csv')
gia_null_nut: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Giacomelli_NULL_Nutlin_2018.csv')
gia_wt_nut: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Giacomelli_WT_Nutlin_2018.csv')
kot_hum: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Kotler_2018.csv')



gia_null_eto_auf: pd.DataFrame = dc.aufteilung_mut_pos(gia_null_eto)
gia_null_nut_auf: pd.DataFrame = dc.aufteilung_mut_pos(gia_null_nut)
gia_wt_nut_auf: pd.DataFrame = dc.aufteilung_mut_pos(gia_wt_nut)
kot_hum_auf: pd.DataFrame = dc.aufteilung_mut_pos(kot_hum)


lowest_vals = pd.DataFrame(columns=['Name of the Dataset', 'Location of lowest DMS_score sum', "Sum", "Original AA"])

lowest_vals['Name of the Dataset'] = ["Giacomelli Null Etoposide", "Giacomelli NULL Nutlin", "Giacomelli WT Nutlin", "Kotler"]
lowest_vals['Location of lowest DMS_score sum'] = [280, 205, 245, 245]
lowest_vals['Sum'] = [-33.450339, -27.798457, -41.124490, -4.352254]
lowest_vals['Original AA'] = ['R', 'Y', 'G', "G"]

lowest_vals_gesammelt = pd.DataFrame(columns=['Giacomelli NULL Etoposide location','Giacomelli NULL Etoposide sums', 'Giacomelli NULL Nutlin location', 'Giacomelli NULL Nutlin sums', "Giacomelli WT Nutlin location", "Giacomelli WT Nutlin sums",  "Kotler location", "Kotler sums"])

anzahl: int = 10

GNELLV = dc.low_val(gia_null_eto, anzahl)
GNESLV = dc.low_val(gia_null_eto, anzahl)
GNNLLV = dc.low_val(gia_null_nut, anzahl)
GNNSLV = dc.low_val(gia_null_nut, anzahl)
GWNLLV = dc.low_val(gia_wt_nut, anzahl)
GWNSLV = dc.low_val(gia_wt_nut, anzahl)
KLLV = dc.low_val(kot_hum, anzahl)
KSLV = dc.low_val(kot_hum, anzahl)

lowest_vals_gesammelt = pd.DataFrame()

lowest_vals_gesammelt["Giacomelli NULL Etoposide location"] = GNELLV.index.get_level_values(0).values
lowest_vals_gesammelt["Giacomelli NULL Etoposide sums"] = GNELLV['Sum'].values
lowest_vals_gesammelt["Giacomelli NULL Nutlin location"] = GNNLLV.index.get_level_values(0).values
lowest_vals_gesammelt["Giacomelli NULL Nutlin sums"] = GNNLLV['Sum'].values
lowest_vals_gesammelt["Giacomelli WT Nutlin location"] = GWNLLV.index.get_level_values(0).values
lowest_vals_gesammelt["Giacomelli WT Nutlin sums"] = GWNLLV['Sum'].values
lowest_vals_gesammelt["Kotler location"] = KLLV.index.get_level_values(0).values
lowest_vals_gesammelt["Kotler sums"] = KSLV['Sum'].values

###########################