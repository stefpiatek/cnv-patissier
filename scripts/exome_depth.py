"""

"""

import glob
import subprocess
import os
import csv

import toml

from . import utils, base_classes

cnv_pat_dir = utils.get_cnv_patissier_dir()


class ExomeDepthBase(base_classes.BaseCNVTool):
    def __init__(self, cohort, gene, start_time, normal_panel):
        super().__init__(cohort, gene, start_time, normal_panel)

        self.settings = {**self.settings, "docker_image": "stefpiatek/exomedepth:1.1.10", "min_mapq": 20}

    def run_command(self, args):
        """Create dir for output and runs a GATK command in docker"""
        print(f"*** Running  ExomeDepth: {self.run_type} \n output: {args[-1]} ***")

        self.run_docker_subprocess(["Rscript", f"/mnt/cnv-caller-resources/exome-depth/{self.run_type}.R", *args])
        print(f"*** Completed  ExomeDepth: {self.run_type} {args[-1]} ***")


class ExomeDepthCohort(ExomeDepthBase):
    def __init__(self, cohort, gene, start_time):
        super().__init__(cohort, gene, start_time, normal_panel=True)
        self.run_type = "exome-depth_cohort"
        self.output_base, self.docker_output_base = self.base_output_dirs()

    def run_workflow(self):
        # write bam locations to file to be read by R script
        try:
            os.makedirs(self.output_base)
        except FileExistsError:
            pass

        with open(f"{self.output_base}/bam_table.txt", "w") as handle:
            writer = csv.DictWriter(handle, fieldnames=["path"], delimiter="\t")
            writer.writeheader()
            for bam in self.settings["bams"]:
                writer.writerow({"path": bam})
        bam_table = f"{self.docker_output_base}/bam_table.txt"

        self.run_command(
            [
                f"--bam-table={bam_table}",
                f"--capture-bed={self.settings['capture_path']}",
                f"--ref-fasta={self.settings['ref_fasta']}",
                f"--min-mapq={self.settings['min_mapq']}",
                "--out-path",
                f"{self.docker_output_base}/cohort.Rdata",
            ]
        )


class ExomeDepthCase(ExomeDepthBase):
    def __init__(self, cohort, gene, start_time):
        super().__init__(cohort, gene, start_time, normal_panel=False)
        self.run_type = "exome-depth_case"
        self.output_base, self.docker_output_base = self.base_output_dirs()

        normal_panel_start = self.get_normal_panel_time()
        self.normal_path_base = (
            f"/mnt/output/{self.cohort}/{normal_panel_start}/{self.run_type.replace('case', 'cohort')}/{self.gene}"
        )

        self.settings = {**self.settings, "normal_panel_start_time": normal_panel_start}

    def run_workflow(self):
        # write bam locations to file to be read by R script
        try:
            os.makedirs(self.output_base)
        except FileExistsError:
            pass

        for bam in self.settings["bams"]:
            sample_name = bam.replace(".bam", "").replace(self.sample_suffix, "").split("/")[-1]
            self.run_command(
                [
                    f"--bam={bam}",
                    f"--sample-name={sample_name}",
                    f"--cohort-rdata={self.normal_path_base}/cohort.Rdata",
                    f"--min-mapq={self.settings['min_mapq']}",
                    "--out-base",
                    f"{self.docker_output_base}/{sample_name}",
                ]
            )
