import glob
import subprocess
import csv
import os

def get_cnv_patissier_dir():
    """Returns the base directory of the project from the scripts folder"""
    return os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

cnv_pat_dir = get_cnv_patissier_dir()

class SymlinkBams:
    def __init__(self, cohort):
        self.cohort = cohort
        self.gene_lists = glob.glob(
            f"{cnv_pat_dir}/input/{self.cohort}/sample-lists/*"
        )

    def symlink_gene_bams(self, gene_list, cnv_status):
        with open(gene_list) as handle:
            samples = csv.DictReader(handle, delimiter="\t")
            for sample in samples:
                if sample["result_type"] == cnv_status:
                    base_output_dir = (f"{cnv_pat_dir}/input/{self.cohort}/bam/"
                                       f"{sample['target_gene']}")
                    output_dir = f"{base_output_dir}/{cnv_status}"
                    if not os.path.exists(output_dir):
                        os.makedirs(output_dir)
                    subprocess.run(
                        ["ln", "-s", f"{sample['sample_name']}", output_dir], check=True
                    )
                    subprocess.run(
                        ["ln", "-s", f"{sample['sample_name']}.bai", output_dir], check=True
                    )

    def main(self):
        for gene_list in self.gene_lists:
            self.symlink_gene_bams(gene_list, "normal")
            self.symlink_gene_bams(gene_list, "normal-panel")
            self.symlink_gene_bams(gene_list, "positive")


if __name__ == "__main__":
    symlinker = SymlinkBams("ICR")
    symlinker.main()
