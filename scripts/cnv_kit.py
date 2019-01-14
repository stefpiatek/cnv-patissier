"""
Default settings following suggested germline alterations from 
https://cnvkit.readthedocs.io/en/stable/germline.html


"""
import csv
import subprocess
import os

import toml

from . import utils, base_classes


class CNVKit(base_classes.BaseCNVTool):
    def __init__(self, capture, gene, start_time):
        super().__init__(capture, gene, start_time, normal_panel=True)
        self.run_type = "cnvkit"
        self.extra_db_fields = ["probes", "cn", "log2", "depth", "weight", ""]

        self.output_base, self.docker_output_base = self.base_output_dirs()

        sample_ids, bams = utils.SampleUtils.select_samples(self.sample_sheet, normal_panel=False)
        self.bam_mount = utils.SampleUtils.get_mount_point(bams)
        docker_bams = [f"/mnt/bam-input/{bam.split(self.bam_mount)[-1]}" for bam in bams]

        self.settings = {**self.settings, "docker_image": "etal/cnvkit:0.9.5", "unknown_bams": docker_bams}
        self.settings["normal_bams"] = self.settings.pop("bams")

    def parse_output_file(self, file_path, sample_id):
        cnvs = []
        with open(file_path, "r") as handle:
            output = csv.DictReader(handle, delimiter="\t")
            for row in output:
                if row:
                    cnv = dict(row)
                    cnv["chrom"] = cnv.pop("chromosome")
                    cnv.pop("gene")
                    cnv["sample_id"] = sample_id
                    cnv["cn"] = int(cnv["cn"])
                    if cnv["cn"] < 2:
                        cnv["alt"] = "DEL"
                    elif cnv["cn"] > 2:
                        cnv["alt"] = "DUP"
                    else:
                        # skip as has a copy number of 2
                        continue
                    cnvs.append(cnv)
        return cnvs

    def run_cnvkit_command(self, args, stdout=None):
        """Create dir for output and runs a CNV-tool command in docker"""
        try:
            os.makedirs(f"{self.output_base}/")
        except FileExistsError:
            base_classes.logger.info(f"Output folder already exists")

        base_classes.logger.info(f"Running  CNVkit: {args[0]} \n output: {args[-1]}")
        self.run_docker_subprocess(["cnvkit.py", *args], stdout=stdout)
        base_classes.logger.info(f"Completed  CNVkit: {args[0]} {args[-1]}")

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
                f"{self.max_cpu}",
            ]
        )

        for bam in self.settings["unknown_bams"]:
            sample = bam.split("/")[-1].replace(".bam", "")
            sample_name = self.bam_to_sample[bam]
            gene_metric_out = f"{self.output_base}/{sample_name}.txt"

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
