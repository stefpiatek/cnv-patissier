"""
Base class for running of a CNV tool - each tool inherits from this class
"""

import csv
import glob
import os
import subprocess

import toml

from . import utils

cnv_pat_dir = utils.get_cnv_patissier_dir()
genome_fasta_path = "/var/reference_sequences/MiSeq/genome.fa"


class BaseCNVTool:
    def __init__(self, cohort, gene, start_time, normal_panel=True):
        self.start_time = start_time
        self.cohort = cohort
        self.gene = gene
        self.sample_suffix = ".sorted"
        self.gene_list = f"{cnv_pat_dir}/input/{cohort}/sample-sheets/{gene}_samples.txt"

        sample_ids, bams = utils.SampleUtils.select_samples(self.gene_list, normal_panel=normal_panel)

        self.bam_mount = utils.SampleUtils.get_mount_point(bams)

        docker_bams = [f"/mnt/bam-input/{bam.split(self.bam_mount)[-1]}" for bam in bams]

        max_cpu = "30"
        max_mem = "50"

        with open(self.gene_list) as handle:
            sample_sheet = csv.DictReader(handle, dialect="excel", delimiter="\t")
            capture_file = set(row["capture_filename"] for row in sample_sheet)
            assert len(capture_file) == 1, (
                "Single capture file should be used for all samples within a sample sheet"
                f"Gene {gene} has multiple captures for its samples, please fix this and run again "
            )

        self.settings = {
            "bams": docker_bams,
            "ref_fasta": f"/mnt/ref_genome/{genome_fasta_path.split('/')[-1]}",
            "intervals": f"/mnt/input/{cohort}/bed/{gene}.bed",
            "max_cpu": max_cpu,
            "max_mem": max_mem,
            "docker_image": None,
            "chromosome_prefix": "chr",
            "cohort": self.cohort,
            "gene": self.gene,
            "start_time": start_time,
            "capture_path": f"/mnt/input/{cohort}/bed/{capture_file.pop()}",
        }

    def base_output_dirs(self):
        """Returns base directory for output: (system_base, docker_base)"""
        output_base = f"{cnv_pat_dir}/output/{self.cohort}/{self.start_time}/{self.run_type}/{self.gene}"
        docker_output_base = output_base.replace(cnv_pat_dir, "/mnt")

        return (output_base, docker_output_base)

    def delete_unused_runs(self):
        print(f"*** Removing any old or unsuccessful runs for {self.cohort}, {self.run_type} ,{self.gene}***")
        subprocess.run(
            [
                "docker",
                "run",
                "--rm",
                "-v",
                f"{cnv_pat_dir}/output/:/mnt/output/:rw",
                "-v",
                f"{cnv_pat_dir}/cnv-caller-resources/:/mnt/cnv-caller-resources/:ro",
                "frolvlad/alpine-python3",
                "python3.6",
                "/mnt/cnv-caller-resources/alpine-python/remove_directories.py",
                self.cohort,
                self.run_type,
                self.gene,
            ]
        )

    def get_normal_panel_time(self):
        normal_path = f"{cnv_pat_dir}/successful-run-settings/{self.cohort}/{self.run_type.replace('case', 'cohort')}/{self.gene}.toml"
        with open(normal_path) as handle:
            normal_config = toml.load(handle)
        return f"{normal_config['start_time']}"

    @staticmethod
    def parse_vcf_4_2(vcf_path):
        """Parses VCF v4.2, if positive cnv, returns dicts of information within a list"""
        cnvs = []
        output_path = vcf_path.split("output/")[-1]
        cohort, time_start, cnv_caller, gene, *args = output_path.split("/")
        sample_id = args[-1].replace(".vcf", "").replace("_segments", "").replace("_intervals", "")
        with open(vcf_path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue

                fields = ["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "data"]
                row = {field: data for (field, data) in zip(fields, line.split())}

                call_data = {key: value for (key, value) in zip(row["format"].split(":"), row["data"].split(":"))}
                if call_data["GT"] != "0":
                    row["start"] = row.pop("pos")
                    if call_data["GT"] == "1":
                        row["alt"] = "DEL"
                    elif call_data["GT"] == "2":
                        row["alt"] = "DUP"
                    cnv = {
                        **row,
                        "end": row["info"].replace("END=", ""),
                        "cnv_caller": cnv_caller,
                        "gene": gene,
                        "sample_id": sample_id,
                        "cohort": cohort,
                    }

                    cnvs.append(cnv)
        return cnvs

    def run_docker_subprocess(self, args, stdout=None):
        """Run docker subprocess as non root user, mounting input and reference genome dir"""
        ref_genome_dir = os.path.dirname(genome_fasta_path)
        subprocess.run(
            [
                "docker",
                "run",
                "--rm",
                "-v",
                f"{ref_genome_dir}:/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/:ro",
                "-v",
                f"{cnv_pat_dir}/input:/mnt/input/:ro",
                "-v",
                f"{self.bam_mount}:/mnt/bam-input/:ro",
                "-v",
                f"{cnv_pat_dir}/cnv-caller-resources/:/mnt/cnv-caller-resources/:ro",
                "-v",
                f"{cnv_pat_dir}/output/:/mnt/output/:rw",
                self.settings["docker_image"],
                *args,
            ],
            check=True,
            stdout=stdout,
        )

    def run_required(self, previous_run_settings_path):
        """Returns True if workflow hasn't been run before or settings have changed since then"""
        if not os.path.exists(previous_run_settings_path):
            print(f"*** Run required for {self.cohort} {self.run_type} {self.gene} ***")
            self.delete_unused_runs()
            return True
        with open(previous_run_settings_path) as handle:
            previous_settings = toml.load(handle)
            previous_settings.pop("start_time")
            current_settings = dict(self.settings)
            current_settings.pop("start_time")
            if current_settings != previous_settings:
                print(f"*** Run required for {self.cohort} {self.run_type} {self.gene} ***")
                self.delete_unused_runs()
                return True
        print(f"*** Re-run not required for {self.cohort} {self.run_type} {self.gene} ***")

    def run_workflow(self):
        """Placeholder for individual tool running"""
        pass

    def main(self):
        previous_run_settings_path = (
            f"{cnv_pat_dir}/successful-run-settings/{self.cohort}/{self.run_type}/{self.gene}.toml"
        )
        if self.run_required(previous_run_settings_path):
            self.run_workflow()
            self.write_settings_toml()
        # TODO: read vcfs into database

    def write_settings_toml(self):
        """Write case toml data for successful run"""
        output_dir = f"{cnv_pat_dir}/successful-run-settings/{self.cohort}/{self.run_type}/"
        try:
            os.makedirs(output_dir)
        except FileExistsError:
            print(f"*** run config directory exists for {self.cohort}/{self.run_type} ***")
        output_path = f"{output_dir}/{self.gene}.toml"
        with open(output_path, "w") as out_file:
            toml.dump(self.settings, out_file)
