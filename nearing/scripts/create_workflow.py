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

NEARING_DATASET_NAMES = "datasets/dataset_name_mapping.csv"
PROCESSED_DATASET_DIR = "all_results"


with open("templates/main_snakefile.jinja2", "r") as f:
    MAIN_SNKFILE_TEMPLATE = Template(f.read())

with open("templates/dataset_snakefile.jinja2", "r") as f:
    DATASET_SNKFILE_TEMPLATE = Template(f.read())

with open("templates/dataset_config.yaml.jinja2", "r") as f:
    DATASET_CFG_TEMPLATE = Template(f.read())

with open("config/main_config.yaml", "r") as f:
    ALL_TOOLS = yaml.safe_load(f)["tools"]


def process_single_dataset(
    logger: logging.Logger,
    dataset_name: str,
):
    tbl_fpath = f"{PROCESSED_DATASET_DIR}/{dataset_name}/table.biom"
    metadata_fpath = f"{PROCESSED_DATASET_DIR}/{dataset_name}/metadata.tsv"
    metadata = pd.read_table(metadata_fpath, sep="\t", index_col=0)

    covariate_name = metadata.squeeze().name
    logger.info(f"Covariate name: {covariate_name}")

    group_counts = metadata[covariate_name].value_counts()
    logger.info(f"Group counts: {group_counts.to_dict()}")

    max_grp = group_counts.idxmax()
    min_grp = group_counts.idxmin()
    logger.info(f"Using group with more samples as reference: {max_grp}")

    songbird_formula = f"C({covariate_name}, Treatment('{max_grp}'))"

    logger.info("Writing dataset specific Snakefile...")

    dataset_snkfile_text = (
        DATASET_SNKFILE_TEMPLATE
        .render({"dataset_name": dataset_name, "all_tools": ALL_TOOLS})
    )

    with open(f"{PROCESSED_DATASET_DIR}/{dataset_name}/Snakefile", "w") as f:
        f.write(dataset_snkfile_text)

    d = {
        "table_file": tbl_fpath,
        "metadata_file": metadata_fpath,
        "covariate_name": covariate_name,
        "target_name": min_grp,
        "reference_name": max_grp,
        "songbird_formula": songbird_formula
    }
    return d


def main(logger: logging.Logger):
    logger.info("Loading dataset names...")
    dataset_names = pd.read_table(NEARING_DATASET_NAMES, index_col=0, sep=",")
    dataset_names = dataset_names.index.tolist()

    all_dataset_dict = dict()

    logger.info("Processing dataset files...")
    parallel_args = {
        "backend": "multiprocessing",
        "verbose": 100
    }
    dict_list = Parallel(n_jobs=38, **parallel_args)(
        delayed(process_single_dataset)(
            logger,
            dataset,
        )
        for dataset in dataset_names
    )

    for dataset, d in zip(dataset_names, dict_list):
        all_dataset_dict[dataset] = d

    main_snakefile_text = (
        MAIN_SNKFILE_TEMPLATE
        .render({"all_dataset_dict": all_dataset_dict})
    )

    logger.info("Writing main Snakefile...")
    with open("Snakefile", "w") as f:
        f.write(main_snakefile_text)

    cfg_text = (
        DATASET_CFG_TEMPLATE
        .render({"all_dataset_dict": all_dataset_dict})
    )

    os.makedirs("config", exist_ok=True)
    logger.info("Writing dataset config file...")
    with open("config/dataset_config.yaml", "w") as f:
        f.write(cfg_text)


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
