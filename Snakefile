import yaml

from snakemake.utils import min_version
min_version("6.0")


configfile: "config/main_config.yaml"


with open("config/dataset_config.yaml", "r") as f:
    dataset_cfg = yaml.safe_load(f)
all_datasets = list(dataset_cfg.keys())


def create_cfg(dataset_name):
    cfg = config.copy()
    cfg.update(dataset_cfg[dataset_name])
    return cfg

target_files = expand(
    "all_results/{dataset}/figures/{viz}",
    dataset=all_datasets,
    viz=["pca.svg", "spearman_heatmap.svg"]
)
target_files.extend(expand(
    "all_results/{dataset}/figures/{curve}/{curve}.pctile_{pctile}.svg",
    dataset=all_datasets,
    curve=["pr", "roc"],
    pctile=config["log_ratio_feat_pcts"]
))

rule all:
    input:
        target_files,


# Separate Snakefile for each dataset
include: "all_results/Chemerin/Snakefile"
include: "all_results/Blueberry/Snakefile"
include: "all_results/MALL/Snakefile"
include: "all_results/Exercise/Snakefile"
include: "all_results/BISCUIT/Snakefile"
include: "all_results/Office/Snakefile"
include: "all_results/art_scher/Snakefile"
include: "all_results/asd_son/Snakefile"
include: "all_results/cdi_schubert/Snakefile"
include: "all_results/cdi_vincent/Snakefile"
include: "all_results/crc_baxter/Snakefile"
include: "all_results/crc_zeller/Snakefile"
include: "all_results/edd_singh/Snakefile"
include: "all_results/hiv_dinh/Snakefile"
include: "all_results/hiv_lozupone/Snakefile"
include: "all_results/hiv_noguerajulian/Snakefile"
include: "all_results/ibd_gevers/Snakefile"
include: "all_results/ibd_papa/Snakefile"
include: "all_results/ob_goodrich/Snakefile"
include: "all_results/ob_ross/Snakefile"
include: "all_results/ob_turnbaugh/Snakefile"
include: "all_results/ob_zhu/Snakefile"
include: "all_results/par_scheperjans/Snakefile"
include: "all_results/t1d_alkanani/Snakefile"
include: "all_results/t1d_mejialeon/Snakefile"
include: "all_results/GWMC_ASIA_NA/Snakefile"
include: "all_results/GWMC_HOT_COLD/Snakefile"
include: "all_results/sw_sed_detender/Snakefile"
include: "all_results/sw_plastic_frere/Snakefile"
include: "all_results/sed_plastic_hoellein/Snakefile"
include: "all_results/wood_plastic_kesy/Snakefile"
include: "all_results/seston_plastic_mccormick/Snakefile"
include: "all_results/glass_plastic_oberbeckmann/Snakefile"
include: "all_results/sed_plastic_rosato/Snakefile"
include: "all_results/ArcticTransects/Snakefile"
include: "all_results/ArcticFireSoils/Snakefile"
include: "all_results/ArcticFreshwaters/Snakefile"
include: "all_results/Ji_WTP_DS/Snakefile"