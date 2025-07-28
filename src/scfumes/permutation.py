import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import zscore

def avg_multiply_by_group(
    subgroup,
    group_names,
    cell_labels,
    exp_mat,
    exp_mat_columns,
    cutoff_exp,
    met_target,
):
    """
    Average expression per group, z‑score rows, multiply by metabolite targets to obtain a Commu_Score matrix.
    """
    avg_exp_all, prop_exp_all = [], []
    for x in group_names:      
        cells_list = set(subgroup.obs[cell_labels == x].index) #index for a specific cell type; subgroup: cell x gene
        cell_loc = np.where(pd.Series(exp_mat_columns).isin(cells_list))[0]
        cell_met = exp_mat[cell_loc,:]
        avg_cell_mat = sp.csc_matrix(np.mean(cell_met, axis=0).reshape(-1, 1)) # calculate mean for each column, then change the row into column, type: sparse matrix
        avg_exp_all.append(avg_cell_mat)
        mask = (cell_met > cutoff_exp).astype(np.float32)
        
        # Calculate the proportion of values in each row that exceed the cutoff
        proportions_cell = sp.csr_matrix((mask.sum(axis=0) / mask.shape[0]).reshape(-1,1)) # calculate proportion for each column, then change the row into column, type: sparse matrix
        prop_exp_all.append(proportions_cell)

    avg_exp = sp.hstack(avg_exp_all) # sparse matrix
    prop_exp = sp.hstack(prop_exp_all)

    #pd.DataFrame(avg_exp.toarray()).to_csv("avg_exp.tsv",sep="\t",index=False)
    z_scale_avg_exp = sp.csr_matrix(np.nan_to_num(zscore(avg_exp.toarray(), axis=1, ddof=1)))
    #pd.DataFrame(z_scale_avg_exp.toarray()).to_csv("z_avg_exp.tsv",sep="\t",index=False)

    met_data_mat = sp.csr_matrix(np.array(met_target['pValue_merge'].to_list())).T

    commu_score = sp.csr_matrix.multiply(z_scale_avg_exp,met_data_mat)
    print("commu_score: %s" % commu_score)

    return commu_score,prop_exp

def multiply_for_permutation(    
    subgroup,
    group_names,
    cell_labels,
    exp_mat,
    exp_mat_columns,
    met_target
):
    """
    Codes only for permutation testing.
    """
    avg_exp_all = []
    for x in group_names:
        cells_list = set(subgroup.obs[cell_labels == x].index) #index for a specific cell type
        cell_loc = np.where(pd.Series(exp_mat_columns).isin(cells_list))[0]
        cell_met = exp_mat[cell_loc,:]
        avg_cell_mat = sp.csr_matrix(np.mean(cell_met, axis=0).reshape(-1, 1)) # calculate mean for each column, then change the row into column      
        avg_exp_all.append(avg_cell_mat)

    avg_exp = sp.hstack(avg_exp_all)
    z_scale_avg_exp = sp.csr_matrix(np.asarray(np.nan_to_num(zscore(avg_exp.toarray(), axis=1, ddof=1))))  # z-score for each row (for cell types)
    met_data_mat = sp.csr_matrix(np.array(met_target['pValue_merge'].to_list())).T
    shuffled_score = sp.csr_matrix.multiply(z_scale_avg_exp,met_data_mat)
    print("Shuffled_score: %s" % shuffled_score)

    return shuffled_score

def shuffle_cells(exp_mat,cell_labels, seed=42):
    rng = np.random.default_rng(seed)
    indices = np.arange(len(cell_labels))
    return exp_mat[rng.permutation(indices), :]

def permutation_test(
    subgroup,
    cell_labels,
    group_names,
    observed_score,
    exp_mat,
    exp_mat_columns,
    met_target,
    n_perm=1000,
):
    comparison_count = np.zeros_like(observed_score.toarray(), dtype=float)
    
    for i in range(n_perm):

        shuffled_mat = shuffle_cells(exp_mat, cell_labels)
        shuffled_score = multiply_for_permutation(subgroup,group_names, cell_labels,shuffled_mat, exp_mat_columns, met_target)
        comparison_count += (shuffled_score > observed_score)
        p_values = comparison_count / float(n_perm)

        p_vals = sp.csr_matrix(p_values)

    return p_vals
