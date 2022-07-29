#!/home/grahman/miniconda3/envs/qadabra-analyses/bin/python
#SBATCH --chdir=/home/grahman/projects/qadabra-analyses/nearing
#SBATCH --output=/home/grahman/projects/qadabra-analyses/nearing/slurm_out/%x.slurm.out
#SBATCH --partition=short
#SBATCH --cpus-per-task=38
#SBATCH --mem-per-cpu=8G
#SBATCH --time=6:00:00

import logging
import os
import time
import yaml

import biom
from jinja2 import Template
from joblib import Parallel, delayed
import numpy as np
import pandas as pd

NEARING_DATASET_FPATHS = "datasets/sorted_input.tsv"
NEARING_DATASET_NAMES = "datasets/dataset_name_mapping.csv"
PROCESSED_DATASET_DIR = "all_results/"


def process_single_dataset(
    logger: logging.Logger,
    dataset: str,
    dataset_files: list,
):
    this_outdir = f"{PROCESSED_DATASET_DIR}/{dataset}"
    os.makedirs(this_outdir, exist_ok=True)
    logger.info(f"Dataset: {dataset}")

    this_dataset_table_file, this_dataset_metadata_file = [
        x for x in dataset_files
        if dataset in x
        and "No_filt_Results" not in x
        and "rare" not in x
    ]

    new_table_fpath = f"{this_outdir}/table.biom"
    new_metadata_fpath = f"{this_outdir}/metadata.tsv"

    logger.info(f"Original table file: {this_dataset_table_file}")
    logger.info(f"Metadata file: {this_dataset_metadata_file}")

    # Convert TSV -> BIOM
    try:
        tbl_tsv = pd.read_table(this_dataset_table_file, sep="\t", index_col=0)
    except IndexError:
        # At least one of the datasets has a # comment row
        tbl_tsv = pd.read_table(this_dataset_table_file, sep="\t",
                                index_col=0, skiprows=[0])
        logger.info("Skipping first row")

    # Remove taxonomy column from tables that have it
    if "taxonomy" in tbl_tsv.columns:
        logger.info("Dropping taxonomy column")
        tbl_tsv = tbl_tsv.drop(columns=["taxonomy"])

    # Append S_ to sample IDs to coerce to string
    tbl_tsv.columns = [f"S_{x}" for x in tbl_tsv.columns]

    # Load metadata
    this_dataset_metadata = pd.read_table(
        this_dataset_metadata_file,
        sep="\t",
        index_col=0
    ).squeeze()
    this_dataset_metadata.index = [
        f"S_{x}" for x in this_dataset_metadata.index
    ]
    grp_name = this_dataset_metadata.name

    # Remove completely discriminatory taxa
    logger.info("Filtering completely discriminatory taxa...")
    joint_df = tbl_tsv.T.join(this_dataset_metadata)
    gb = joint_df.groupby(grp_name).sum()
    feat_presence = gb.apply(lambda x: x.all())
    feats_to_keep = feat_presence[feat_presence].index
    tbl_tsv = tbl_tsv.loc[feats_to_keep]

    # Convert to BIOM
    tbl_biom = biom.Table(tbl_tsv.values, sample_ids=tbl_tsv.columns,
                          observation_ids=tbl_tsv.index)

    # 10% prevalence filter
    logger.info("Filtering at 10% prevalence...")
    prev = tbl_biom.pa(inplace=False).sum(axis="observation")
    n = int(tbl_biom.shape[1] * 0.1)
    feats_to_keep = tbl_biom.ids("observation")[np.where(prev >= n)]
    tbl_biom.filter(feats_to_keep, "observation")
    tbl_biom.remove_empty()

    # Filter table and metadata to common samples
    samps_to_keep = list(
        set(tbl_biom.ids()).intersection(this_dataset_metadata.index)
    )
    tbl_biom.filter(samps_to_keep)
    logger.info(f"New table shape: {tbl_biom.shape}")

    with biom.util.biom_open(new_table_fpath, "w") as f:
        tbl_biom.to_hdf5(f, "converted")

    this_dataset_metadata.loc[samps_to_keep].to_csv(new_metadata_fpath,
                                                    sep="\t", index=True)


def main(logger: logging.Logger):
    logger.info("Loading dataset filepaths...")
    with open(NEARING_DATASET_FPATHS, "r") as f:
        dataset_files = f.read().splitlines()

    logger.info("Loading dataset names...")
    dataset_names = pd.read_table(NEARING_DATASET_NAMES, index_col=0, sep=",")
    dataset_names = dataset_names.index.tolist()

    dataset_files = [
        os.path.join("datasets", x[1:])  # Remove slash as first character
        for x in dataset_files
    ]

    os.makedirs(f"{PROCESSED_DATASET_DIR}", exist_ok=True)

    logger.info("Processing dataset files...")
    parallel_args = {
        "backend": "multiprocessing",
        "verbose": 100
    }
    Parallel(n_jobs=38, **parallel_args)(
        delayed(process_single_dataset)(
            logger,
            dataset,
            dataset_files,
        )
        for dataset in dataset_names
    )


if __name__ == "__main__":
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    sh = logging.StreamHandler()
    formatter = logging.Formatter(
        "[%(asctime)s - %(levelname)s] :: %(message)s",
    )
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    start_time = time.time()
    main(logger)
    end_time = time.time()
    logger.info(f"Total script time: {end_time - start_time}")
