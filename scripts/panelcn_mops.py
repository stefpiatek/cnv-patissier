"""
panelcn.MOPS set up using basic instructions in manual https://bioconductor.org/packages/release/bioc/vignettes/panelcn.mops/inst/doc/panelcn.mops.pdf

Used splitting of BED file with defaults (100bp bins with 50 bp overlap)
"""

import csv
import pathlib

from . import base_classes


class panelcnMOPS(base_classes.BaseCNVTool):
    def __init__(self, capture, gene, start_time, normal_panel=True):
        self.run_type = "panelcn_mops"
        super().__init__(capture, gene, start_time, normal_panel=normal_panel)
        self.extra_db_fields = ["gene", "exon", "rc", "medrc", "rc.norm", "medrc.norm", "lowqual", "cn"]
        self.settings = {**self.settings, "docker_image": "stefpiatek/panelcn_mops:1.4.0"}

    def parse_output_file(self, file_path, sample_id):
        cnvs = []
        with open(file_path, "r") as handle:
            output = csv.DictReader(handle, delimiter="\t")
            for row in output:
                if row["Sample"] == sample_id:
                    cnv = {key.lower(): value for key, value in row.items()}
                    cnv["chrom"] = f"{self.settings['chromosome_prefix']}{cnv.pop('chr')}"
                    cnv["sample_id"] = cnv.pop("sample")
                    copy_number = int(cnv["cn"].lstrip("CN"))
                    if copy_number < 2:
                        cnv["alt"] = "DEL"
                    elif copy_number > 2:
                        cnv["alt"] = "DUP"
                    else:
                        raise Exception(f"row doesn't have a copy number change {cnv}")
                    cnvs.append(cnv)

        return cnvs

    def run_command(self, args):
        self.run_docker_subprocess(["Rscript", f"/mnt/cnv-caller-resources/panelcn_mops/panelcn_mops_runner.R", *args])

    def run_workflow(self):
        bed_path = self.settings["capture_path"].replace("/mnt", base_classes.cnv_pat_dir)

        pathlib.Path(self.output_base).mkdir(parents=True, exist_ok=True)

        with open(bed_path, "r") as input_bed:
            with open(f"{self.output_base}/capture.bed", "w") as output_bed:
                for line in input_bed:
                    chrom, start, end, gene = line.split()
                    output_bed.write(f"{chrom}\t{start}\t{end}\t{gene}.{chrom}.{start}.{end}\n")

        with open(f"{self.output_base}/samples.tsv", "w") as handle:
            handle.write(f"bam_path\tsample_name\tsample_type\n")
            for bam in self.settings["unknown_bams"] + self.settings["normal_bams"]:
                sample = self.bam_to_sample[bam]
                if bam in self.settings["unknown_bams"]:
                    sample_type = "unknown"
                else:
                    sample_type = "normal_panel"

                handle.write(f"{bam}\t{sample}\t{sample_type}\n")

        self.run_command([f"--output-path={self.docker_output_base}", f"--gene={self.gene}"])

        sample_names = [f"{self.bam_to_sample[unknown_bam]}" for unknown_bam in self.settings["unknown_bams"]]
        output_paths = [f"{self.output_base}/calls.tsv" for sample_name in sample_names]

        return output_paths, sample_names
