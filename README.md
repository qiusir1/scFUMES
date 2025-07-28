# scFUMES

Single cell FUnctional MEtabolite-Sensor (scFUMES), a tool to analyze metabolite‑sensor networks in single‑cell data in the context of diseases. scFUMES will utilize single cell/nuclei data and buit-in metabolite-sensor interaction to prioritize cell type-specific metabolite-sensor network. Single-cell/nuclei data will be preprocessed, but quality controls were hihgly recommanded before running the program. Reference: Cell type-specific master metabolic regulators of Alzheimer’s disease. Yunguang Qiu, Yuan Hou, Liam Wetzel, Jessica Z.K. Caldwell, Xiongwei Zhu, Andrew A. Pieper, Tian Liu, Feixiong Cheng. bioRxiv 2025.07.11.664443; doi: https://doi.org/10.1101/2025.07.11.664443.

## Installation (conda)

```bash
git clone
cd scFUMES
conda env create -f environment.yml
conda activate scFUMES
pip install .
```

## Example

```bash
scFUMES \
  --dataset tests/adata.h5ad \
  --met_target data/MetTarget.tsv \
  --cluster_label "cell_type" \
  --comparison "disease" \
  --group "normal" \
  --outdir results \
  --n_perm 10
```