"""
Exome Depth run with defaults from vignette https://cran.r-project.org/web/packages/ExomeDepth/vignettes/ExomeDepth-vignette.pdf

- Uses GC content normalisation
"""

import csv
import glob
import subprocess
import os

import toml

from . import utils, base_classes


class ExomeDepthBase(base_classes.BaseCNVTool):
    def __init__(self, capture, gene, start_time, normal_panel):
        super().__init__(capture, gene, start_time, normal_panel)

        self.settings = {**self.settings, "docker_image": "stefpiatek/exomedepth:1.1.10", "min_mapq": 20}

    def run_command(self, args):
        """Create dir for output and runs a GATK command in docker"""
        base_classes.logger.info(f"Running  ExomeDepth: {self.run_type} \n output: {args[-1]}")

        self.run_docker_subprocess(["Rscript", f"/mnt/cnv-caller-resources/exome-depth/{self.run_type}.R", *args])
        base_classes.logger.info(f"Completed  ExomeDepth: {self.run_type} {args[-1]}")


class ExomeDepthCohort(ExomeDepthBase):
    def __init__(self, capture, gene, start_time):
        super().__init__(capture, gene, start_time, normal_panel=True)
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
    def __init__(self, capture, gene, start_time):
        super().__init__(capture, gene, start_time, normal_panel=False)
        self.extra_db_fields = [
            "nexons",
            "id",
            "BF",
            "reads.expected",
            "reads.observed",
            "reads.ratio",
            "start.p",
            "end.p",
        ]
        self.run_type = "exome-depth_case"
        self.output_base, self.docker_output_base = self.base_output_dirs()

        normal_panel_start = self.get_normal_panel_time()
        self.normal_path_base = (
            f"/mnt/output/{self.capture}/{normal_panel_start}/{self.run_type.replace('case', 'cohort')}/{self.gene}"
        )

        self.settings = {**self.settings, "normal_panel_start_time": normal_panel_start}

    def parse_output_file(self, file_path, sample_id):
        cnvs = []
        with open(file_path, "r") as handle:
            output = csv.DictReader(handle, delimiter="\t")
            for row in output:
                if row:
                    cnv = dict(row)
                    cnv["chrom"] = f"{self.settings['chromosome_prefix']}{cnv.pop('chromosome')}"
                    cnv["sample_id"] = sample_id
                    call = cnv.pop("type")
                    if call == "deletion":
                        cnv["alt"] = "DEL"
                    elif call == "duplication":
                        cnv["alt"] = "DUP"
                    else:
                        raise Exception(f"non-deletion or duplication called in Exome depth:\n {cnv}")
                    cnvs.append(cnv)
        return cnvs

    def run_workflow(self):
        # write bam locations to file to be read by R script
        try:
            os.makedirs(self.output_base)
        except FileExistsError:
            pass

        for bam in self.settings["bams"]:
            sample_name = self.bam_to_sample[bam]
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
