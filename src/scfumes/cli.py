import argparse
from pathlib import Path
from .calculator import ScFUMESCalculator

def build_parser() -> argparse.ArgumentParser:
    scfumes = argparse.ArgumentParser(
        prog="scFUMES",
        description="Single cell FUnctional MEtabolite-Sensor (scFUMES), \
        a tool to analyze metabolite‑sensor networks in single‑cell data in the context of diseases.",
    )
    scfumes.add_argument("--dataset", required=True, help=".h5ad file")
    scfumes.add_argument("--met_target", required=True,
                   help="Metabolite‑target TSV")
    scfumes.add_argument("--cluster_label", default="cell_type",
                   help="Column in metadata that defines clusters")
    scfumes.add_argument("--comparison", default="Diagnosis",
                   help="comparison of interest, e.g., 'Diagnosis'")
    scfumes.add_argument("--group", default="AD",
                   help="Group name of interest, e.g., 'AD, 'Control'")
    scfumes.add_argument("--outdir", default="results",
                   help="Output directory")
    scfumes.add_argument("--n_perm", type=int, default=1000,
                   help="Number of permutations")
    scfumes.add_argument(
        "--regress_covariates",
        default=None,
        help="Comma-separated list of covariates in adata.obs to regress out",
    )
    return scfumes


def main(argv = None):
    args = build_parser().parse_args(argv)

    regress_covariates = None
    if args.regress_covariates:
        regress_covariates = [
            covariate.strip()
            for covariate in args.regress_covariates.split(",")
            if covariate.strip()
        ]

    calc = ScFUMESCalculator.from_files(
        h5ad_path=args.dataset,
        cluster_label=args.cluster_label,
    )
    calc.process(
        met_target_file=args.met_target,
        comparison=args.comparison,
        group=args.group,
        output_dir=Path(args.outdir),
        n_perm=args.n_perm,
        regress_covariates=regress_covariates,
    )
