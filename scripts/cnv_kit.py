"""
Default settings following suggested germline alterations from 
https://cnvkit.readthedocs.io/en/stable/germline.html


"""
import subprocess
import os

import toml

from . import utils, base_classes

cnv_pat_dir = utils.get_cnv_patissier_dir()


class CNVKit(base_classes.BaseCNVTool):
    def __init__(self, cohort, gene, start_time):
        super().__init__(cohort, gene, start_time, normal_panel=True)
        self.run_type = "cnvkit"

        self.output_base, self.docker_output_base = self.base_output_dirs()

        sample_ids, bams = utils.SampleUtils.select_samples(self.gene_list, normal_panel=False)
        self.bam_mount = utils.SampleUtils.get_mount_point(bams)
        docker_bams = [f"/mnt/bam-input/{bam.split(self.bam_mount)[-1]}" for bam in bams]

        self.settings = {**self.settings, "docker_image": "etal/cnvkit:0.9.5", "unknown_bams": docker_bams}
        self.settings["normal_bams"] = self.settings.pop("bams")

    def run_cnvkit_command(self, args, stdout=None):
        """Create dir for output and runs a CNV-tool command in docker"""
        try:
            os.makedirs(f"{self.output_base}/")
        except FileExistsError:
            print(f"*** Output folder already exists ***")

        print(f"*** Running  CNVkit: {args[0]} \n output: {args[-1]} ***")
        self.run_docker_subprocess(["cnvkit.py", *args], stdout=stdout)
        print(f"*** Completed  CNVkit: {args[0]} {args[-1]} ***")

    def run_workflow(self):
        self.run_cnvkit_command(
            [
                "batch",
                *self.settings["unknown_bams"],
                "--normal",
                *self.settings["normal_bams"],
                "--targets",
                self.settings["capture_path"],
                "--fasta",
                self.settings["ref_fasta"],
                "--output-reference",
                f"{self.docker_output_base}/reference.cnn",
                "--output-dir",
                f"{self.docker_output_base}/batch-results/",
                "-p",
                f"{self.settings['max_cpu']}",
            ]
        )

        for bam_path in self.settings["unknown_bams"]:
            sample = bam_path.split("/")[-1].replace(".bam", "")
            gene_metric_out = f"{self.output_base}/{sample.replace(self.sample_suffix, '')}.txt"

            # segmetrics to get confidence intervals for filtering in call
            self.run_cnvkit_command(
                [
                    "segmetrics",
                    "-s",
                    f"{self.docker_output_base}/batch-results/{sample}.cns",
                    "--ci",
                    f"{self.docker_output_base}/batch-results/{sample}.cnr",
                    "-o",
                    f"{self.docker_output_base}/batch-results/{sample}_ci.cns",
                ]
            )

            self.run_cnvkit_command(
                [
                    "call",
                    "-m",
                    "clonal",
                    f"{self.docker_output_base}/batch-results/{sample}_ci.cns",
                    "--filter",
                    "ci",
                    "-o",
                    f"{self.docker_output_base}/batch-results/{sample}_calls.cns",
                ]
            )
