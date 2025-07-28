from pathlib import Path
import scanpy as sc
import pandas as pd

from .preprocess import prepare_adata
from .met_target import met_tar_net
from .cutoff import auto_cutoff


class ScFUMESCalculator:
    """
    Integrating preprocessing, network building, and permutation  testing.
    """

    def __init__(self, adata, cluster_label="cell_type"):
        self.adata = adata
        self.cluster_label = cluster_label

    @classmethod
    def from_files(cls, h5ad_path, cluster_label="cell_type"):
        adata = sc.read_h5ad(h5ad_path)
        prepare_adata(adata, cluster_label)
        return cls(adata, cluster_label)

    def process(
        self,
        met_target_file="MetTarget.tsv",
        comparison="Diagnosis",
        group="AD",
        output_dir="results",
        q=0.25,
        n_perm=1000,
    ):
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        met_target = pd.read_csv(met_target_file, sep="\t")
        met_target = met_target[met_target['Gene_name'].isin(self.adata.var_names) & met_target['INCHIKEY']] # Gene names should be usingthe same IDs.
        met_target = pd.DataFrame(met_target)
        print(met_target)
        
        s_m_all = met_target['Gene_name'] + ' ~ ' + met_target['INCHIKEY']
        s_m_all = s_m_all.tolist()
        s_m_all_name = met_target['Gene_name'] + ' ~ ' + met_target['NAME']
        s_m_all_name = s_m_all_name.tolist()
        target = [s.split(' ~ ')[0] for s in s_m_all] # target (gene) list
        target_locs = [self.adata.var_names.tolist().index(t) for t in target]
        target_mat = self.adata[:, target_locs] 

        cutoff_exp = auto_cutoff(target_mat, q=q)

        subgroup = target_mat[target_mat.obs[comparison].isin([group])]
        met_tar_net(
            subgroup=subgroup,
            met_target=met_target,
            cutoff_exp=cutoff_exp,
            output_name= comparison + "_" + group,
            out_dir=out_dir,
            n_perm=n_perm,
            s_m_all_name = s_m_all_name,
            cluster_label=self.cluster_label,
        )
