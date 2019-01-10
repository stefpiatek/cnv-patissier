"""
CopywriteR set up using suggested settings for targeted exome 
https://bioconductor.org/packages/release/bioc/vignettes/CopywriteR/inst/doc/CopywriteR.pdf
    - 50 kb window

Note: only has paired sample mode, so normal samples are just randomly paired with unknowns
"""

import csv
import subprocess
import os
import pathlib

from . import utils, base_classes


class Copywriter(base_classes.BaseCNVTool):
    def __init__(self, capture, gene, start_time, normal_panel=True):
        super().__init__(capture, gene, start_time, normal_panel)
        self.run_type = "copywriter"

        self.output_base, self.docker_output_base = self.base_output_dirs()

        sample_ids, bams = utils.SampleUtils.select_samples(self.gene_list, normal_panel=False)
        self.bam_mount = utils.SampleUtils.get_mount_point(bams)
        docker_bams = [f"/mnt/bam-input/{bam.split(self.bam_mount)[-1]}" for bam in bams]

        self.settings = {**self.settings, "docker_image": "stefpiatek/copywriter:2.2.0", "unknown_bams": docker_bams}
        self.settings["normal_bams"] = self.settings.pop("bams")

    def parse_output_file(self, file_path, sample_id):
        cnvs = []
        with open(file_path, "r") as handle:
            output = csv.DictReader(handle, delimiter="\t")
            sample_to_bam = {sample: bam for (bam, sample) in self.bam_to_sample.items()}
            bam_name = pathlib.Path(sample_to_bam[sample_id]).name
            bamfile_to_sample = {pathlib.Path(bam_path).name: sample for bam_path, sample in self.bam_to_sample.items()}
            for row in output:
                if row['unknown'] == bam_name:
                    cnv = dict(row)
                    cnv["chrom"] = f"{self.settings['chromosome_prefix']}{cnv['chrom']}"
                    cnv["sample_id"] = sample_id
                    cnv["control_id"] = bamfile_to_sample[cnv.pop('control')]
                    cnv["seg.mean"] = float(cnv["seg.mean"])
                    # TODO: set this threshold using ROC curve
                    if cnv["seg.mean"] <= -1.3:
                        cnv["alt"] = "DEL"
                    elif cnv["seg.mean"] >= 1.3:
                        cnv["alt"] = "DUP"
                    else:
                        continue
                    for field in ["num.mark", "unknown"]:
                        cnv.pop(field)
                    cnvs.append(cnv)
        return cnvs

    def run_command(self, args):
        self.run_docker_subprocess(["Rscript", f"/mnt/cnv-caller-resources/copywriter/copywriter_runner.R", *args])

    def run_workflow(self):
        # write bam locations to file to be read by R script
        # only paired settings so just pair up unknowns with controls at random
        # if it looks good, could use exomedepth choice of normal to select the appropriate control

        # assume 3 gb per worker so memory doesn't run out
        max_workers = min([int(self.settings["max_cpu"]), int(self.settings["max_mem"]) // 3])
        total_batches = len(self.settings["unknown_bams"]) // 30
        if len(self.settings["unknown_bams"]) % 30:
            total_batches += 1
        for batch_number in range(total_batches):
            batch_output = f"{self.output_base}/batch_{batch_number}"
            docker_batch_output = f"{self.docker_output_base}/batch_{batch_number}"
            try:
                os.makedirs(batch_output)
            except FileExistsError:
                pass

            with open(f"{batch_output}/all_samples.txt", "w") as handle:
                writer = csv.DictWriter(handle, fieldnames=["samples", "controls"], delimiter="\t")
                writer.writeheader()
                for normal_index in range(30):
                    unknown_index = normal_index + 30 * batch_number
                    normal_bam = self.settings["normal_bams"][normal_index]
                    try:
                        unknown_bam = self.settings["unknown_bams"][unknown_index]
                    except IndexError:
                        # will get index error for the last batch if not divisible by 30. Skip these
                        continue
                    writer.writerow({"samples": normal_bam, "controls": normal_bam})
                    writer.writerow({"samples": unknown_bam, "controls": normal_bam})

            base_classes.logger.info(f"Running CopywriteR on batch {batch_number}")
            self.run_command(
                [
                    f"--max-cpu={max_workers}",
                    f"--output-path={docker_batch_output}",
                    f"--capture-regions={self.settings['capture_path']}",
                ]
            )
