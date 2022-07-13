#!/home/grahman/miniconda3/envs/qadabra-analyses/bin/python
#SBATCH --chdir=/home/grahman/projects/qadabra-analyses
#SBATCH --output=/home/grahman/projects/qadabra-analyses/%x.slurm.out
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --time=6:00:00

import logging
import os
import time

import biom
from jinja2 import Template
import numpy as np
import pandas as pd

LOGFILE = "create_workflow.py.log"
NEARING_DIR = "data/nearing_datasets/"
NEARING_DATASET_FPATHS = os.path.join(NEARING_DIR, "sorted_input.tsv")
NEARING_DATASET_NAMES = os.path.join(NEARING_DIR, "dataset_name_mapping.csv")
SNKFILE_TEMPLATE = "templates/Snakefile.jinja2"
CFG_TEMPLATE = "templates/config.yaml.jinja2"


def main(logger: logging.Logger):
    logger.info("Loading dataset filepaths...")
    with open(NEARING_DATASET_FPATHS, "r") as f:
        dataset_files = f.read().splitlines()

    logger.info("Loading dataset names...")
    dataset_names = pd.read_table(NEARING_DATASET_NAMES, index_col=0, sep=",")
    dataset_names = dataset_names.index.tolist()

    dataset_files = [
        os.path.join(NEARING_DIR, x[1:])  # Remove slash as first character
        for x in dataset_files
    ]

    os.makedirs("all_results", exist_ok=True)
    all_dataset_dict = dict()
    logger.info("Processing dataset files...")
    for dataset in dataset_names:
        os.makedirs(f"all_results/{dataset}", exist_ok=True)
        logger.info("=========================================================")
        logger.info(f"Dataset: {dataset}")

        this_dataset_table_file, this_dataset_metadata_file = [
            x for x in dataset_files
            if dataset in x
            and "No_filt_Results" not in x
            and "rare" not in x
        ]

        new_table_fpath = f"all_results/{dataset}/table.biom"
        new_metadata_fpath = f"all_results/{dataset}/metadata.tsv"

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
        tbl_biom = biom.Table(tbl_tsv.values, sample_ids=tbl_tsv.columns,
                              observation_ids=tbl_tsv.index)
        this_dataset_metadata = pd.read_table(
            this_dataset_metadata_file,
            sep="\t",
            index_col=0
        ).squeeze()
        this_dataset_metadata.index = [
            f"S_{x}" for x in this_dataset_metadata.index
        ]

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
        logger.info(f"Saved converted table to {new_table_fpath}!")

        this_dataset_metadata = this_dataset_metadata.loc[samps_to_keep]
        logger.info(f"Metadata shape: {this_dataset_metadata.shape}")

        this_dataset_metadata.to_csv(new_metadata_fpath, sep="\t", index=True)
        logger.info(f"Saved filtered table to {new_metadata_fpath}")

        covariate_name = this_dataset_metadata.name
        logger.info(f"Covariate name: {covariate_name}")

        group_counts = this_dataset_metadata.value_counts()
        logger.info(f"Group counts: {group_counts.to_dict()}")

        max_grp = group_counts.idxmax()
        min_grp = group_counts.idxmin()
        logger.info(f"Using group with more samples as reference: {max_grp}")

        all_dataset_dict[dataset] = {
            "table_file": new_table_fpath,
            "metadata_file": new_metadata_fpath,
            "covariate_name": covariate_name,
            "target_name": min_grp,
            "reference_name": max_grp
        }

    logger.info("=========================================================")

    with open(SNKFILE_TEMPLATE, "r") as f:
        snakefile_text = (
            Template(f.read())
            .render({"all_dataset_dict": all_dataset_dict})
        )

    logger.info("Writing Snakefile...")
    os.makedirs("workflow", exist_ok=True)
    with open("workflow/Snakefile", "w") as f:
        f.write(snakefile_text)

    with open(CFG_TEMPLATE, "r") as f:
        cfg_text = (
            Template(f.read())
            .render({"all_dataset_dict": all_dataset_dict})
        )

    os.makedirs("config", exist_ok=True)
    logger.info("Writing config file...")
    with open("config/config.yaml", "w") as f:
        f.write(cfg_text)


if __name__ == "__main__":
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(LOGFILE, mode="w+")
    formatter = logging.Formatter(
        "[%(asctime)s - %(levelname)s] :: %(message)s",
    )
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    start_time = time.time()
    main(logger)
    end_time = time.time()
    logger.info(f"Total script time: {end_time - start_time}")
