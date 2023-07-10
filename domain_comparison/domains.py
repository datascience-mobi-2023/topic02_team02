# all added functions are to be declared in the init file, to ensure effortless usage
import pandas as pd
from scipy.stats import shapiro
import data_cleanup as dc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import gaussian_kde

domains_regulatory: dict = {'t1_domain': (1, 43),
                            't2_domain': (44, 63),
                            'pr_domain': (64, 92),
                            'dna_b_domain': (102, 292),
                            'tetra_domain': (320, 355),
                            'reg_domain': (356, 393)
                            }


def slice_domain(dms_data: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    """Takes in a Data set in the form of the DataFrames in the DMS_data folder.
    Also takes in positions to be sliced as start and end position."""
    mutations_df_copy = dms_data.copy()
    mutations_df = dc.aufteilung_mut_pos(mutations_df_copy)
    sliced_domains = mutations_df.loc[
        (mutations_df['position_mut'] >= start) & (mutations_df['position_mut'] <= end)].copy()

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


def adjust_domain(regulatory: dict, name: str, frame: pd.DataFrame) -> dict:
    domain = dc.rmv_na(dc.df_transform(doc.slice_domain(frame, start=regulatory[name][0], end=regulatory[name][1])))
    domain_list = doc.slice_domain(frame, start=regulatory[name][0], end=regulatory[name][1])
    mean = domain.mean().rename('mean')
    res: dict = {'domain': domain, 'domain_list': domain_list, 'mean': mean}
    return res


def distr_and_hmap(domain_reg: dict, domain_name: str) -> None:
    dms_scores = domain_reg['domain_list']['DMS_score']
    plt.hist(dms_scores, bins=50)
    plt.xlabel('DMS Score')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of DMS Scores in the {domain_name} Domain')
    plt.show()

    print(f'Mean: {dms_scores.mean()}')
    print(f'Median: {dms_scores.median()}')
    sns.heatmap(domain_reg['domain'])

    return None


def comp_dms_distr(regulatory: dict, frame) -> None:
    regions = [adjust_domain(regulatory, elem, frame)['domain_list']['DMS_score'] for elem in regulatory.keys()]
    dms_scores = np.concatenate(regions)

    fig, ax = plt.subplots(figsize=(10, 6))

    datasets = ['Transactivation domain 1', 'Transactivation domain 2', 'Proline rich region',
                'DNA binding domain', 'Tetramerization domain', 'Regulatory domain']
    for data, label in zip(regions, datasets):
        kde = gaussian_kde(data)
        x_vals = np.linspace(np.min(dms_scores), np.max(dms_scores), 1000)
        y_vals = kde(x_vals)

        ax.plot(x_vals, y_vals, linewidth=2, label=label)

    ax.set_xlabel('Value')
    ax.set_ylabel('Density')
    ax.set_title('Different Domains with only Single Mutation DMS Scores')
    ax.legend()

    plt.show()
    return None
