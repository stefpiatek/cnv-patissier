"""
SavvyCNV set up using default settings from the readme 
https://github.com/rdemolgen/SavvySuite


"""
import csv
import pathlib

from scripts import base_classes


class SavvyCNV(base_classes.BaseCNVTool):
    def __init__(self, capture, gene, start_time, normal_panel=True):
        self.run_type = "savvycnv"
        super().__init__(capture, gene, start_time, normal_panel=normal_panel)
        self.extra_db_fields = [
            "supporting_chunks",
            "width_of_cnv",
            "phred_score",
            "normalised_phred",
            "relative_dosage",
            "coverage_filename",
        ]
        self.settings = {**self.settings, "docker_image": "stefpiatek/savvycnv:f996a83"}

    def parse_output_file(self, file_path, sample_id):
        call_fields = ["chrom", "start", "end", "call", *self.extra_db_fields]
        coverage_filename = f"{self.docker_output_base}CoverageBinner/{sample_id}.coverageBinner"
        cnvs = []

        with open(file_path, "r") as handle:
            output = csv.DictReader(handle, delimiter="\t", fieldnames=call_fields)
            for row in output:
                if row["coverage_filename"] == coverage_filename:
                    cnv = dict(row)
                    cnv["chrom"] = f"{self.settings['chromosome_prefix']}{cnv['chrom'].lstrip('chr')}"
                    call = cnv.pop("call")
                    if call == "Deletion":
                        cnv["alt"] = "DEL"
                    elif call == "Duplication":
                        cnv["alt"] = "DUP"
                    else:
                        raise Exception(f"No call has been made for {coverage_filename}")
                    cnvs.append(cnv)
        return cnvs

    def run_command(self, args, stdout):
        """Runs SavvyCNV command in docker"""

        base_classes.logger.info(f"Running {args[0]} step of SavvyCNV for {args[-1]} ")
        self.run_docker_subprocess(
            ["java", f"-Xmx{self.max_mem}g", f"-XX:ConcGCThreads={self.max_cpu}", *args], stdout=stdout
        )
        base_classes.logger.info(f"Completed {args[0]} step of SavvyCNV for {args[-1]}")

    def run_workflow(self):
        coverage_dir = pathlib.Path(f"{self.output_base}/CoverageBinner")
        coverage_dir.mkdir(parents=True, exist_ok=True)

        docker_coverage_files = []
        for bam in self.settings["unknown_bams"] + self.settings["normal_bams"]:
            sample_name = self.bam_to_sample[bam]
            docker_coverage_files.append(f"{self.docker_output_base}/CoverageBinner/{sample_name}.coverageBinner")
            with open(coverage_dir / f"{sample_name}.coverageBinner", "w") as handle:
                self.run_command(["CoverageBinner", bam], stdout=handle)

        with open(f"{self.output_base}/cnv_calls.txt", "w") as handle:
            self.run_command(["SavvyCNV", "-d", "400", "-trans", "0.01", *docker_coverage_files], stdout=handle)

        sample_names = [f"{self.bam_to_sample[unknown_bam]}" for unknown_bam in self.settings["unknown_bams"]]
        output_paths = [f"{self.output_base}/cnv_calls.txt" for sample_name in sample_names]

        return output_paths, sample_names
