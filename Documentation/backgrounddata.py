import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import data_cleanup as dc
import functions as fun

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

gia_null_eto_df: pd.DataFrame = dc.df_transform(gia_null_eto)
gia_null_nut_df: pd.DataFrame = dc.df_transform(gia_null_nut)
gia_wt_nut_df: pd.DataFrame = dc.df_transform(gia_wt_nut)
kot_hum_df: pd.DataFrame = dc.df_transform(kot_hum)


#lowest_vals = pd.DataFrame(columns=['Name of the Dataset', 'Location of lowest DMS_score sum', "Sum", "Original AA"])

#lowest_vals['Name of the Dataset'] = ["Giacomelli Null Etoposide", "Giacomelli NULL Nutlin", "Giacomelli WT Nutlin", "Kotler"]
#lowest_vals['Location of lowest DMS_score sum'] = [280, 205, 245, 245]
#lowest_vals['Sum'] = [-33.450339, -27.798457, -41.124490, -4.352254]
#lowest_vals['Original AA'] = ['R', 'Y', 'G', "G"]

#lowest_vals_gesammelt = pd.DataFrame(columns=['Giacomelli NULL Etoposide location','Giacomelli NULL Etoposide sums', 'Giacomelli NULL Nutlin location', 'Giacomelli NULL Nutlin sums', "Giacomelli WT Nutlin location", "Giacomelli WT Nutlin sums",  "Kotler location", "Kotler sums"])

#anzahl: int = 10

#GNELLV = dc.low_val(gia_null_eto, anzahl)
#GNESLV = dc.low_val(gia_null_eto, anzahl)
#GNNLLV = dc.low_val(gia_null_nut, anzahl)
#GNNSLV = dc.low_val(gia_null_nut, anzahl)
#GWNLLV = dc.low_val(gia_wt_nut, anzahl)
#GWNSLV = dc.low_val(gia_wt_nut, anzahl)
#KLLV = dc.low_val(kot_hum, anzahl)
#KSLV = dc.low_val(kot_hum, anzahl)

#lowest_vals_gesammelt = pd.DataFrame()

#lowest_vals_gesammelt["Giacomelli NULL Etoposide location"] = GNELLV.index.get_level_values(0).values
#lowest_vals_gesammelt["Giacomelli NULL Etoposide sums"] = GNELLV['Sum'].values
#lowest_vals_gesammelt["Giacomelli NULL Nutlin location"] = GNNLLV.index.get_level_values(0).values
#lowest_vals_gesammelt["Giacomelli NULL Nutlin sums"] = GNNLLV['Sum'].values
#lowest_vals_gesammelt["Giacomelli WT Nutlin location"] = GWNLLV.index.get_level_values(0).values
#lowest_vals_gesammelt["Giacomelli WT Nutlin sums"] = GWNLLV['Sum'].values
#lowest_vals_gesammelt["Kotler location"] = KLLV.index.get_level_values(0).values
#lowest_vals_gesammelt["Kotler sums"] = KSLV['Sum'].values



gia_null_eto_mean = fun.df_mean(gia_null_eto_df)
gia_null_nut_mean = fun.df_mean(gia_null_nut_df)
gia_wt_nut_mean = fun.df_mean(gia_wt_nut_df)
kot_hum_mean = fun.df_mean(kot_hum_df)

gia_null_eto_mean.name = 'gia_null_eto_mean'
gia_null_nut_mean.name = 'gia_null_nut_mean'
gia_wt_nut_mean.name = 'gia_wt_nut_mean'
kot_hum_mean.name = 'kot_hum_mean'
# Erwartet eine Liste an dataframes ("dataframes=[gia_null_eto_mean, gia_null_nut_mean, gia_wt_nut_mean, kot_hum_mean]")
#def plot_multiple_datasets(dataframes):
    # figsize kann verändert werden je nach Präferenz
    #plt.figure(figsize=(40, 6))
    # for-loop der jeden Eintrag in der Liste "Dataframes" durchgeht
    #for df in dataframes:
        # Werte für x-Achse aus 'position_mut' als int extrahiert
        #x = df.columns.get_level_values('position_mut').astype(int)
        # Werte für y-Achse aus erster Zeile von df als float extrahiert
        #y = df.iloc[0].values.astype(float)
        #label = df.name
        #plt.plot(x, y, label=label, marker='o')
    #plt.xlabel('position')
    #plt.ylabel('DMS_score')
    #plt.title('Mean DMS_scores of different datasets throughout each position of tp53')
    #WARUM ZUR VERDAMMTEN HÖLLE MACHT ER NUR 102 BIS 292 ICH VERSTEHE ES NICHT BITTE BEENDE ES
    #plt.xticks(x[::10])
    #plt.legend()
    #plt.show()
    #plt.plot(x[::10], y[::10], label=label, marker='o')












###########################