"""
CODEX2 using settins from targeted NGS demo script
https://github.com/yuchaojiang/CODEX2/blob/master/targeted_sequencing/CODEX2_targeted_demo.R

Used 10 K iterations as suggested by their overall readme as the maximum required

"""
import csv
import pathlib

from scripts import base_classes


class CODEX2(base_classes.BaseCNVTool):
    def __init__(self, capture, gene, start_time, normal_panel=True):
        self.run_type = "codex2"
        super().__init__(capture, gene, start_time, normal_panel=normal_panel)
        self.extra_db_fields = ["gene", "length_kb", "length_exon", "raw_cov", "norm_cov", "copy_no", "lratio", "mBIC"]
        self.settings = {**self.settings, "docker_image": "stefpiatek/codex2:26e796c"}

    def parse_output_file(self, file_path, sample_id):
        cnvs = []
        with open(file_path, "r") as handle:
            output = csv.DictReader(handle, delimiter="\t")
            for row in output:
                if row["sample_name"] == sample_id:
                    cnv = dict(row)
                    cnv["chrom"] = f"{self.settings['chromosome_prefix']}{cnv.pop('chr').replace('chr', '')}"
                    cnv["start"] = cnv.pop("st_bp")
                    cnv["end"] = cnv.pop("ed_bp")
                    cnv["sample_id"] = cnv.pop("sample_name")
                    cnv["alt"] = cnv.pop("cnv").upper()

                    cnvs.append(cnv)
        return cnvs

    def run_command(self, args):
        """Runs CODEX2 script in docker"""
        self.run_docker_subprocess(["Rscript", "/mnt/cnv-caller-resources/codex2/run_codex.R", *args])

    def run_workflow(self):
        pathlib.Path(self.output_base).mkdir(parents=True, exist_ok=True)

        with open(f"{self.output_base}/samples.tsv", "w") as handle:
            handle.write(f"bam_path\tsample_name\tsample_type\n")
            for bam in self.settings["unknown_bams"] + self.settings["normal_bams"]:
                sample = self.bam_to_sample[bam]
                if bam in self.settings["unknown_bams"]:
                    sample_type = "unknown"
                else:
                    sample_type = "normal_panel"

                handle.write(f"{bam}\t{sample}\t{sample_type}\n")

        self.run_command(
            [
                f"--output-path={self.docker_output_base}",
                f"--capture-bed={self.settings['capture_path']}",
                f"--chrom-prefix={self.settings['chromosome_prefix']}",
            ]
        )


        sample_names = [f"{self.bam_to_sample[unknown_bam]}" for unknown_bam in self.settings["unknown_bams"]]
        output_paths = [f"{self.output_base}/calls.txt" for sample_name in sample_names]

        return output_paths, sample_names