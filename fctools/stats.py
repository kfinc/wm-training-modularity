import numpy as np
from scipy import stats
import pandas as pd

def m_wei(x, w):
    """Weighted Mean"""
    return np.sum(x * w) / np.sum(w)

def cov_wei(x, y, w):
    """Weighted Covariance"""
    return np.sum(w * (x - m_wei(x, w)) * (y - m_wei(y, w))) / np.sum(w)

def corr_wei(x, y, w):
    """Weighted Correlation"""
    return cov_wei(x, y, w) / np.sqrt(cov_wei(x, x, w) * cov_wei(y, y, w))

def bootstrap_replicate_1d(data, func):
    """Calculate bootstrap replicate"""
    return func(np.random.choice(data, size=len(data)))

def ttest_rel_cond(group, variable, data):
    """ Function which calculates paired t-test comparing variable between two task conditions"""
    sess = ['ses-1', 'ses-2', 'ses-3', 'ses-4']
    conds = ['1-back', '2-back']
    table = np.zeros((4,2))

    for i, ses in enumerate(sess):
        group_1 = data[(data.Session == ses) & (data.Condition == '1-back') & (data.Group == group)]
        group_2 = data[(data.Session == ses) & (data.Condition == '2-back') & (data.Group == group)]

        t, pval = stats.ttest_rel(group_1[variable].values, group_2[variable].values,)
        table[i, 0] = t
        table[i, 1] = pval

    return pd.DataFrame(table.tolist(), columns = ['statistic', 'pval'])


def ttest_rel_sess(group, variable, data):
    """ Function which calculates paired t-test comparing variable between each session"""
    sess = ['ses-1', 'ses-2', 'ses-3', 'ses-4']
    table = np.zeros((4,4,2))
    
    for i, ses1 in enumerate(sess):
        for j, ses2 in enumerate(sess):
            ses_1 = data[(data.Session == ses1) & (data.Group == group)]
            ses_2 = data[(data.Session == ses2) & (data.Group == group)]

            t, pval = stats.ttest_rel(ses_1[variable].values, ses_2[variable].values,)
            table[i, j, 0] = t
            table[i, j, 1] = pval

    return table

def ttest_ind_groups(variable, data):
    """ Function which calculates two sample t-test comparing variables between groups"""
    sess = ['ses-1', 'ses-2', 'ses-3', 'ses-4']

    table = np.zeros((4,2))

    for i, ses in enumerate(sess):
        group_1 = data[(data.Session == ses) & (data.Group == 'Experimental')]
        group_2 = data[(data.Session == ses) & (data.Group == 'Control')]

        t, pval = stats.ttest_ind(group_1[variable].values, group_2[variable].values,)
        table[i, 0] = t
        table[i, 1] = pval
        
    return  pd.DataFrame(table.tolist(), columns = ['statistic', 'pval'])

