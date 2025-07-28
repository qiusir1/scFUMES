import numpy as np
import pandas as pd
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests
from datetime import datetime

from .permutation import (
    avg_multiply_by_group,
    permutation_test,
)

def met_tar_net(
    subgroup,
    met_target,
    cutoff_exp,
    output_name,
    out_dir,
    n_perm=1000,
    s_m_all_name=None,
    cluster_label="cell_type",
):
    ncell = len(subgroup.obs_names) # cell, row
    ngene = len(subgroup.var_names) # gene, column
    print('We get expression data with {n1} genes and {n2} cells.'.format(n1 = ngene, n2 = ncell))
    
    exp_mat = sp.csr_matrix(subgroup.X)
    exp_mat_columns = subgroup.obs_names
    cell_labels = subgroup.obs[cluster_label] # Cell type labels are defaulted in the adata.obs under the column name "cell_type"
    group_names = cell_labels.unique().tolist()
    

    commu_score, prop_exp = avg_multiply_by_group(
        subgroup, group_names, cell_labels, exp_mat,
        exp_mat_columns, cutoff_exp, met_target
    )

    print(f"[{datetime.now():%F %T}] Start permutations…")
    p_vals = permutation_test(
        subgroup, cell_labels, group_names, commu_score,
        exp_mat, exp_mat_columns, met_target, n_perm=n_perm
    )
    
    p_values_flat = p_vals.toarray().flatten()
    _, p_values_corrected_flat, _, _ = multipletests(p_values_flat, method='fdr_bh')
    
    fdr_vals = p_values_corrected_flat.reshape(p_vals.shape)
    
    print(f"[{datetime.now():%F %T}] Permutations Completed.")
    
    target_exp_prop = pd.DataFrame(prop_exp.toarray(), columns = group_names)
    target_exp_prop.index = s_m_all_name
    target_exp_prop.to_csv(out_dir / f"{output_name}_Prop.tsv", sep='\t')
    
    target_exp_prop_unstacked = (target_exp_prop).unstack().reset_index()
    target_exp_prop_unstacked.columns = ['Cell type', 'metabolite-target pairs', 'Target_Prop']
    target_exp_prop_unstacked.to_csv(out_dir / f"{output_name}_Prop_unstacked.tsv", sep='\t')
    
    metcellScore_data = pd.DataFrame(commu_score.toarray(), columns = group_names)
    metcellScore_data.index = s_m_all_name
    metcellScore_data.to_csv(out_dir / f"{output_name}_Commu_score.tsv", sep='\t')
    metcellScore_data_unstack = metcellScore_data.unstack().reset_index()
    metcellScore_data_unstack.columns = ['Cell type', 'metabolite-target pairs', 'Commu_score']
    metcellScore_data_unstack.to_csv(out_dir / f"{output_name}_Commu_score_unstacked.tsv", sep = '\t')
    
    p_values = pd.DataFrame(p_vals.toarray(),columns = group_names)
    p_values.index = s_m_all_name
    #p_values.columns = ['Cell type', 'metabolite-target pairs', 'p value']
    p_values.to_csv(out_dir / f"{output_name}_p_value_perm.tsv", sep = '\t')
    p_values_unstacked = (p_values).unstack().reset_index()
    p_values_unstacked.columns = ['Cell type', 'metabolite-target pairs', 'p value']
    p_values_unstacked.to_csv(out_dir / f"{output_name}_p_value_perm_unstacked.tsv", sep = '\t')

    fdr_vals_data = pd.DataFrame(fdr_vals, columns = group_names, index= s_m_all_name)
    fdr_vals_data.to_csv(out_dir / f"{output_name}_fdr_perm.tsv", sep="\t")
    fdr_vals_data_unstacked = (fdr_vals_data).unstack().reset_index()
    fdr_vals_data_unstacked.columns = ['Cell type', 'metabolite-target pairs', 'fdr_bh']
    fdr_vals_data_unstacked.to_csv(out_dir / f"{output_name}_fdr_perm_unstacked.tsv", sep = '\t')


    p_values_unstacked = p_values_unstacked[["p value"]]
    fdr_vals_data_unstacked = fdr_vals_data_unstacked[["fdr_bh"]]
    target_exp_prop_unstacked = target_exp_prop_unstacked[["Target_Prop"]]

    con1 = pd.concat([metcellScore_data_unstack, p_values_unstacked], axis=1)
    con2 = pd.concat([con1, fdr_vals_data_unstacked], axis=1)
    con3 = pd.concat([con2, target_exp_prop_unstacked], axis=1)
    con3["fdr_bh"] = con3["fdr_bh"].replace({0.000:2.2e-16})
    con3["neglogfdr_bh"] = -np.log10(con3["fdr_bh"])  

    con3.to_csv(out_dir / f"{output_name}_fdr_perm_merge.tsv",sep="\t")

    print(f"[{datetime.now():%F %T}] Done – results in {out_dir}")
