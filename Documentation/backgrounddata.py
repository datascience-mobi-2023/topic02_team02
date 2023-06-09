#IMPORTS
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import data_cleanup as dc
import functions as fun
import resources as rs
import severity_score as ss

###################
#DATA FOR VISUALIZATION

gia_null_eto: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Giacomelli_NULL_Etoposide_2018.csv')
gia_null_nut: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Giacomelli_NULL_Nutlin_2018.csv')
gia_wt_nut: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Giacomelli_WT_Nutlin_2018.csv')
kot_hum: pd.DataFrame = pd.read_csv('../DMS_data/P53_HUMAN_Kotler_2018.csv')

gia_null_eto_norm: pd.DataFrame = dc.norm(gia_null_eto)
gia_null_nut_norm: pd.DataFrame = dc.norm(gia_null_nut)
gia_wt_nut_norm: pd.DataFrame = dc.norm(gia_wt_nut)
kot_hum_norm: pd.DataFrame = dc.norm(kot_hum)

gia_null_eto_norm_amp: pd.DataFrame = dc.aufteilung_mut_pos(gia_null_eto_norm)
gia_null_nut_norm_amp: pd.DataFrame = dc.aufteilung_mut_pos(gia_null_nut_norm)
gia_wt_nut_norm_amp: pd.DataFrame = dc.aufteilung_mut_pos(gia_wt_nut_norm)
kot_hum_norm_amp: pd.DataFrame = dc.aufteilung_mut_pos(kot_hum_norm)

#gia_null_eto_auf: pd.DataFrame = dc.aufteilung_mut_pos(gia_null_eto)
#gia_null_nut_auf: pd.DataFrame = dc.aufteilung_mut_pos(gia_null_nut)
#gia_wt_nut_auf: pd.DataFrame = dc.aufteilung_mut_pos(gia_wt_nut)
#kot_hum_auf: pd.DataFrame = dc.aufteilung_mut_pos(kot_hum)

gia_null_eto_df: pd.DataFrame = dc.df_transform(gia_null_eto)
gia_null_nut_df: pd.DataFrame = dc.df_transform(gia_null_nut)
gia_wt_nut_df: pd.DataFrame = dc.df_transform(gia_wt_nut)
kot_hum_df: pd.DataFrame = dc.df_transform(kot_hum)

#kot_hum_auf_norm: pd.DataFrame = dc.norm(kot_hum_auf)

lowest_vals = pd.DataFrame(columns=['Name of the Dataset', 'Location of lowest DMS_score sum', "Sum", "Original AA"])

lowest_vals['Name of the Dataset'] = ["Giacomelli Null Etoposide", "Giacomelli NULL Nutlin", "Giacomelli WT Nutlin", "Kotler"]
lowest_vals['Location of lowest DMS_score sum'] = [280, 205, 245, 245]
lowest_vals['Sum'] = [-6.190289, -13.762829, -15.419176, -6.568038]
lowest_vals['Original AA'] = ['R', 'Y', 'G', "G"]

lowest_vals_gesammelt = pd.DataFrame(columns=['Giacomelli NULL Etoposide location','Giacomelli NULL Etoposide sums', 'Giacomelli NULL Nutlin location', 'Giacomelli NULL Nutlin sums', "Giacomelli WT Nutlin location", "Giacomelli WT Nutlin sums",  "Kotler location", "Kotler sums"])

anzahl: int = 10
lowest_vals_gesammelt = pd.DataFrame()

#####

GNELV= dc.low_val(gia_null_eto_norm, 5)
GNEHV= dc.high_val(gia_null_eto_norm, 5)



####
gia_null_eto_mean = fun.df_mean(gia_null_eto_df)
gia_null_nut_mean = fun.df_mean(gia_null_nut_df)
gia_wt_nut_mean = fun.df_mean(gia_wt_nut_df)
kot_hum_mean = fun.df_mean(kot_hum_df)

gia_null_eto_mean.name = 'gia_null_eto_mean'
gia_null_nut_mean.name = 'gia_null_nut_mean'
gia_wt_nut_mean.name = 'gia_wt_nut_mean'
kot_hum_mean.name = 'kot_hum_mean'
###########################

# Liniendiagrammcode:
#%%
# RNA SEQUENCE AND SLICING
#%%
mutated_p53 = ss.generate_codon_variations(ss.p53_codons_gia)
aa = ss.translate_codons_df(mutated_p53)
# PART 1 NUR SINGLE MUTATIONS
#%%
# 1.2.) INVERSE DF_TRAFO : NOTE : NO LONGER NEEDED, CHANGED FUNCTION SELECT_SMUT
#%%
gia_null_eto_dfi: pd.DataFrame = dc.df_transform(gia_null_eto).T
gia_null_nut_dfi: pd.DataFrame = dc.df_transform(gia_null_nut).T
gia_wt_nut_dfi: pd.DataFrame = dc.df_transform(gia_wt_nut).T
kot_hum_dfi: pd.DataFrame = dc.df_transform(kot_hum).T
#%%
# 1.3.) CREATE VARIATION_MATRIX AND CLEAN IT
#%%

variation_matrix_gia = ss.translate_codons_df(ss.generate_codon_variations(ss.p53_codons_gia))
cleaned_vm_gia = ss.clean_variation_matrix(variation_matrix_gia)
cleaned_vm_kot = cleaned_vm_gia.iloc[101:292].copy()

#%%
# 1.4.) HIER WERDEN BIS AUF SINGLE MUTATIONS ALLE ENTFERNT
#%%
gia_null_eto_dfi_cvm_aa: pd.DataFrame = ss.select_smut(gia_null_eto_dfi, cleaned_vm_gia)
gia_null_nut_cvm_aa: pd.DataFrame = ss.select_smut(gia_null_nut_dfi, cleaned_vm_gia)
gia_wt_nut_cvm_aa: pd.DataFrame = ss.select_smut(gia_wt_nut_dfi, cleaned_vm_gia)

#%%
# PART 2 NICHT NUR SINGLE MUTATIONS
#%%
# 2.1.) Z-TRAFO
#%%
gia_null_eto_z: pd.DataFrame = dc.z_transform(gia_null_eto)
gia_null_nut_z: pd.DataFrame = dc.z_transform(gia_null_nut)
gia_wt_nut_z: pd.DataFrame = dc.z_transform(gia_wt_nut)
kot_hum_z: pd.DataFrame = dc.z_transform(kot_hum)
#%%
# 2.2.) MIN MAX NORMIERUNG
#%%
gia_null_eto_z_mmn: pd.DataFrame = dc.min_max_norm(gia_null_eto_z)
gia_null_nut_z_mmn: pd.DataFrame = dc.min_max_norm(gia_null_nut_z)
gia_wt_nut_z_mmn: pd.DataFrame = dc.min_max_norm(gia_wt_nut_z)
kot_hum_z_mmn: pd.DataFrame = dc.min_max_norm(kot_hum_z)
#%%
# 2.3.) DF TRAFO
#%%
gia_null_eto_z_mmn_df: pd.DataFrame = dc.df_transform(gia_null_eto_z_mmn)
gia_null_nut_z_mmn_df: pd.DataFrame = dc.df_transform(gia_null_nut_z_mmn)
gia_wt_nut_z_mmn_df: pd.DataFrame = dc.df_transform(gia_wt_nut_z_mmn)
kot_hum_z_mmn_df: pd.DataFrame = dc.df_transform(kot_hum_z_mmn)
#%%
# 2.4.) MEAN AN JEDER STELLE
#%%
gia_null_eto_z_mmn_df_mean: pd.DataFrame = fun.df_mean(gia_null_eto_z_mmn_df)
gia_null_nut_z_mmn_df_mean: pd.DataFrame = fun.df_mean(gia_null_nut_z_mmn_df)
gia_wt_nut_z_mmn_df_mean: pd.DataFrame = fun.df_mean(gia_wt_nut_z_mmn_df)
kot_hum_z_mmn_df_mean: pd.DataFrame = fun.df_mean(kot_hum_z_mmn_df)
#%%
# 2.5.) LINIENGRAPH
#%%
kot_hum_z_mmn_df_mean.name = 'kot_hum_z_mmn_df_mean'
gia_null_nut_z_mmn_df_mean.name = 'gia_null_nut_z_mmn_df_mean'
gia_wt_nut_z_mmn_df_mean.name = 'gia_wt_nut_z_mmn_df_mean'
gia_null_eto_z_mmn_df_mean.name = 'gia_null_eto_z_mmn_df_mean'

def liniengraph(dataframes):
    plt.figure(figsize=(35, 7))
    for df in dataframes:
        x = df.columns.get_level_values('position_mut').astype(int)
        y = df.iloc[0].values.astype(float)
        label = df.name
        plt.plot(x, y, label=label, marker='o')
    plt.xlabel('position')
    plt.ylabel('DMS_score')
    plt.title('Mean DMS_scores of different datasets throughout each position of tp53')
    plt.xticks(np.arange(0, 393+1, 10))
    plt.legend()
    plt.show()

# mean distance
#%% md
## GIA NULL ETO
#%%
fpathGNE = '../DMS_data/P53_HUMAN_Giacomelli_NULL_Etoposide_2018.csv'
dfGNE = pd.read_csv(fpathGNE)
mutations_dfGNE = dc.aufteilung_mut_pos(dfGNE)
mutations_dfGNE_norm: pd.DataFrame = dc.norm(mutations_dfGNE)
subs_df = mutations_dfGNE_norm.groupby(["AS_old", "AS_new"])
mean_scoresGNE = subs_df.DMS_score.mean()
mean_scores_dfGNE = mean_scoresGNE.reset_index()
mean_substitutionsGNE = mean_scores_dfGNE.pivot(index="AS_old", columns="AS_new", values="DMS_score")
dc.rmv_na(mean_substitutionsGNE)
#%% md
## GIA NULL NUT
#%%
fpathGNN = '../DMS_data/P53_HUMAN_Giacomelli_NULL_Nutlin_2018.csv'
dfGNN = pd.read_csv(fpathGNN)
mutations_dfGNN = dc.aufteilung_mut_pos(dfGNN)
mutations_dfGNN_norm: pd.DataFrame = dc.norm(mutations_dfGNN)
subs_df = mutations_dfGNN_norm.groupby(["AS_old", "AS_new"])
mean_scoresGNN = subs_df.DMS_score.mean()
mean_scores_dfGNN = mean_scoresGNN.reset_index()
mean_substitutionsGNN = mean_scores_dfGNN.pivot(index="AS_old", columns="AS_new", values="DMS_score")
dc.rmv_na(mean_substitutionsGNN)
#%% md
## GIA WT NUT
#%%
fpathGWN = '../DMS_data/P53_HUMAN_Giacomelli_WT_Nutlin_2018.csv'
dfGWN = pd.read_csv(fpathGWN)
mutations_dfGWN = dc.aufteilung_mut_pos(dfGWN)
mutations_dfGWN_norm: pd.DataFrame = dc.norm(mutations_dfGWN)
subs_df = mutations_dfGWN_norm.groupby(["AS_old", "AS_new"])
mean_scoresGWN = subs_df.DMS_score.mean()
mean_scores_dfGWN = mean_scoresGWN.reset_index()
mean_substitutionsGWN = mean_scores_dfGWN.pivot(index="AS_old", columns="AS_new", values="DMS_score")
dc.rmv_na(mean_substitutionsGWN)
#%% md
## KOT HUM
#%%
fpathKH = '../DMS_data/P53_HUMAN_Giacomelli_WT_Nutlin_2018.csv'
dfKH = pd.read_csv(fpathKH)
mutations_dfKH = dc.aufteilung_mut_pos(dfKH)
mutations_dfKH_norm: pd.DataFrame = dc.norm(mutations_dfKH)
subs_df = mutations_dfKH_norm.groupby(["AS_old", "AS_new"])
mean_scoresKH = subs_df.DMS_score.mean()
mean_scores_dfKH = mean_scoresKH.reset_index()
mean_substitutionsKH = mean_scores_dfKH.pivot(index="AS_old", columns="AS_new", values="DMS_score")
dc.rmv_na(mean_substitutionsKH)
###############
# Domain comparison diagramms