import logging
import os
import glob

from jinja2 import Template
import pandas as pd

LOGFILE = f"{__file__}.log"
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

    all_dataset_dict = dict()
    logger.info("Processing dataset files...")
    for dataset in dataset_names:
        logger.info("=========================================================")
        logger.info(f"Dataset: {dataset}")

        this_dataset_table_file, this_dataset_metadata_file = [
            x for x in dataset_files
            if dataset in x
            and "No_filt_Results" not in x
            and "rare" not in x
        ]
        logger.info(f"Table file: {this_dataset_table_file}")
        logger.info(f"Metadata file: {this_dataset_metadata_file}")

        this_dataset_metadata = pd.read_table(
            this_dataset_metadata_file,
            sep="\t",
            index_col=0
        ).squeeze()
        logger.info(f"Metadata shape: {this_dataset_metadata.shape}")

        covariate_name = this_dataset_metadata.name
        logger.info(f"Covariate name: {covariate_name}")

        group_counts = this_dataset_metadata.value_counts()
        logger.info(f"Group counts: {group_counts.to_dict()}")

        max_grp = group_counts.idxmax()
        min_grp = group_counts.idxmin()
        logger.info(f"Using group with more samples as reference: {max_grp}")

        all_dataset_dict[dataset] = {
            "table_file": this_dataset_table_file,
            "metadata_file": this_dataset_metadata_file,
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
    fh = logging.FileHandler(LOGFILE)
    formatter = logging.Formatter(
        "[%(asctime)s - %(levelname)s] :: %(message)s",
    )
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    main(logger)
