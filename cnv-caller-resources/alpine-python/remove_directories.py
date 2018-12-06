import glob
import os
import shutil
from argparse import ArgumentParser


def remove_empty_folder(folder):
    try:
        os.rmdir(folder)
    except OSError as e:
        if "Directory not empty" not in str(e):
            raise e


if __name__ == "__main__":
    parser = ArgumentParser(description="Delete unused runs")
    parser.add_argument("cohort")
    parser.add_argument("run_type")
    parser.add_argument("gene")
    args = parser.parse_args()

    for gene_folder in glob.glob(f"/mnt/output/{args.cohort}/*/{args.run_type}/{args.gene}/"):
        shutil.rmtree(gene_folder)
    for run_type_folder in glob.glob(f"/mnt/output/{args.cohort}/*/{args.run_type}"):
        remove_empty_folder(run_type_folder)
    for date_folder in glob.glob(f"/mnt/output/{args.cohort}/*"):
        remove_empty_folder(date_folder)
