"""
DECoN using default settings, no exclusion of failed samples
"""

import csv
import pathlib

from . import utils, base_classes


class DECoN(base_classes.BaseCNVTool):
    """
    DECoN class, main in BaseCNVTool will cause self.run_workflow() to be triggered.
    """
    def __init__(self, capture, gene, start_time, normal_panel=True):
        self.run_type = "decon"
        super().__init__(capture, gene, start_time, normal_panel=normal_panel)
        self.extra_db_fields = [
            "CNV.ID",
            "Correlation",
            "N.comp",
            "Start.b",
            "End.b",
            "N.exons",
            "Genomic.ID",
            "BF",
            "Reads.expected",
            "Reads.observed",
            "Reads.ratio",
            "Gene",
        ]
        self.settings = {**self.settings, "docker_image": "stefpiatek/decon:1.0.2"}

    def parse_output_file(self, file_path, sample_id):
        sample_to_bam = {sample: bam for (bam, sample) in self.bam_to_sample.items()}
        bam_name = pathlib.Path(sample_to_bam[sample_id]).name.replace(".bam", "")
        bamfile_to_sample = {
            pathlib.Path(bam_path).name.replace(".bam", ""): sample for bam_path, sample in self.bam_to_sample.items()
        }
        cnvs = []
        with open(file_path, "r") as handle:
            output = csv.DictReader(handle, delimiter="\t")
            for row in output:
                if row["Sample"] == bam_name:
                    cnv = dict(row)
                    cnv["chrom"] = f"{self.settings['chromosome_prefix']}{cnv.pop('Chromosome').lstrip('chr')}"
                    cnv["start"] = cnv.pop("Start")
                    cnv["end"] = cnv.pop("End")
                    cnv["sample_id"] = sample_id
                    cnv["sample_id"] = bamfile_to_sample[cnv.pop("Sample")]
                    call = cnv.pop("CNV.type")
                    if call == "deletion":
                        cnv["alt"] = "DEL"
                    elif call == "duplication":
                        cnv["alt"] = "DUP"

                    cnvs.append(cnv)
        return cnvs

    def run_command(self, args):
        base_classes.logger.info(f"Running  {self.run_type}: {args[0]} \n output: {args[-1]}")
        self.run_docker_subprocess(["Rscript", *args])

        base_classes.logger.info(f"Completed  {self.run_type}: {args[0]} {args[-1]}")

    def run_workflow(self):
        """
        Will run entire workflow and return the final output data paths, and the sample names analysed
        :return: (output_paths, sample_names)

        """
        pathlib.Path(self.output_base).mkdir(parents=True, exist_ok=True)

        with open(f"{self.output_base}/bams.txt", "w") as handle:
            for bam in self.settings["normal_bams"] + self.settings["unknown_bams"]:
                handle.write(f"{bam}\n")

        read_in = f"{self.docker_output_base}/ReadInBams"
        self.run_command(
            [
                "ReadInBams.R",
                "--bams",
                f"{self.docker_output_base}/bams.txt",
                "--bed",
                self.settings["capture_path"],
                "--fasta",
                self.settings["ref_fasta"],
                "--out",
                read_in,
            ]
        )

        cnv_calls = f"{self.docker_output_base}/cnv_calls"
        self.run_command(["makeCNVcalls.R", "--Rdata", f"{read_in}.RData", "--plot", "None", "--out", cnv_calls])

        sample_names = [f"{self.bam_to_sample[unknown_bam]}" for unknown_bam in self.settings["unknown_bams"]]
        output_paths = [f"{self.output_base}/cnv_calls_all.txt" for unknown_bam in self.settings["unknown_bams"]]

        return output_paths, sample_names
