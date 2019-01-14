"""
Base class for running of a CNV tool - each tool inherits from this class
"""

import csv
import glob
import os
import subprocess
import sys

import toml
from loguru import logger

from scripts import utils, models
from settings import cnv_pat_settings
from scripts.db_session import DbSession

cnv_pat_dir = utils.get_cnv_patissier_dir()
# set up logger, replacing built in stderr and adding cnv-patissier log file
logger.remove(0)
logger.add(
    sys.stderr,
    colorize=True,
    format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> <level>{level}</level> {message}",
    level="DEBUG",
)
docker_level = logger.level("DOCKER", no=15, color="<yellow>")
logger.add(f"{cnv_pat_dir}/logs/cnv-patissier.log")
logger.add(f"{cnv_pat_dir}/logs/error.log", level="ERROR", mode="w")


class BaseCNVTool:
    def __init__(self, capture, gene, start_time, normal_panel=True):
        self.session = DbSession.factory()
        self.start_time = start_time
        self.capture = capture
        self.extra_db_fields = []
        self.gene = gene
        self.max_cpu = cnv_pat_settings["max_cpu"]
        self.max_mem = cnv_pat_settings["max_mem"]
        self.sample_sheet = f"{cnv_pat_dir}/input/{capture}/sample-sheets/{gene}_samples.txt"

        sample_ids, bams = utils.SampleUtils.select_samples(self.sample_sheet, normal_panel=normal_panel)
        bam_to_sample = utils.SampleUtils.get_bam_to_id(self.sample_sheet)

        self.bam_mount = utils.SampleUtils.get_mount_point(bams)

        docker_bams = [f"/mnt/bam-input/{bam.split(self.bam_mount)[-1]}" for bam in bams]
        self.bam_to_sample = {
            f"/mnt/bam-input/{bam.split(self.bam_mount)[-1]}": sample_id for (bam, sample_id) in bam_to_sample.items()
        }
        self.sample_to_bam = {sample_id: bam for (bam, sample_id) in self.bam_to_sample.items()}

        with open(self.gene_list) as handle:
            sample_sheet = csv.DictReader(handle, dialect="excel", delimiter="\t")
            capture_file = set(f"{row['capture_name']}.bed" for row in sample_sheet)
            assert len(capture_file) == 1, (
                "Single capture file should be used for all samples within a sample sheet"
                f"Gene {gene} has multiple captures for its samples, please fix this and run again "
            )

        self.settings = {
            "bams": docker_bams,
            "ref_fasta": f"/mnt/ref_genome/{cnv_pat_settings['genome_fasta_path'].split('/')[-1]}",
            "intervals": f"/mnt/input/{capture}/bed/{gene}.bed",
            "docker_image": None,
            "chromosome_prefix": "chr",
            "capture": self.capture,
            "gene": self.gene,
            "start_time": start_time,
            "capture_path": f"/mnt/input/{capture}/bed/{capture_file.pop()}",
        }

    def base_output_dirs(self):
        """Returns base directory for output: (system_base, docker_base)"""
        output_base = f"{cnv_pat_dir}/output/{self.capture}/{self.start_time}/{self.run_type}/{self.gene}"
        docker_output_base = output_base.replace(cnv_pat_dir, "/mnt")

        return (output_base, docker_output_base)

    def delete_unused_runs(self):
        logger.info(f"Removing any old or unsuccessful runs for {self.capture}, {self.run_type}, {self.gene}")
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
                self.capture,
                self.run_type,
                self.gene,
            ]
        )

    def filter_capture(self):
        capture_path = f"{cnv_pat_dir}/{self.settings['capture_path'][5:]}"
        gene_bed = []
        with open(capture_path, "r") as capture:
            for line in capture:
                name = line.split()[3]
                if name == self.gene:
                    gene_bed.append(line)
        return gene_bed

    def filter_cnvs(self, cnvs, gene_bed):
        filtered_cnvs = []
        for cnv in cnvs:
            for line in gene_bed:
                chrom, start, end, *_ = line.split()
                start, end = int(start), int(end)
                cnv_start, cnv_end = int(cnv["start"]), int(cnv["end"])
                if cnv["chrom"] == chrom:
                    within_region = (start <= cnv_start <= end) or (start <= cnv_end <= end)
                    spanning_region = (cnv_start <= start) and (cnv_end >= end)
                    if within_region or spanning_region:
                        cnv["json_data"] = {field: cnv.pop(field) for field in self.extra_db_fields if field}
                        filtered_cnvs.append(cnv)
                        break
        return filtered_cnvs

    def get_bam_header(self, sample_id):
        docker_bam = self.sample_to_bam[sample_id]
        samtools = self.run_docker_subprocess(
            ["samtools", "view", "-H", docker_bam], docker_image="stefpiatek/samtools:1.9", stdout=subprocess.PIPE
        )
        header = str(samtools.stdout, "utf-8")
        return header

    def get_normal_panel_time(self):
        normal_path = (
            f"{cnv_pat_dir}/successful-run-settings/{self.capture}/"
            f"{self.run_type.replace('case', 'cohort')}/{self.gene}.toml"
        )
        with open(normal_path) as handle:
            normal_config = toml.load(handle)
        return f"{normal_config['start_time']}"

    def parse_output_file(file_path, sample_id=None):
        """Dummy method to be overwritten by each CNV-caller class"""
        pass

    @staticmethod
    def parse_vcf_4_2(input_vcf, sample_id):
        """Parses VCF v4.2, if positive cnv, returns dicts of information within a list"""
        cnvs = []
        for line in input_vcf:
            if line.startswith("#") or not line:
                continue
            fields = ["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "data"]
            row = {field: data for (field, data) in zip(fields, line.split())}
            format_data = {key: value for (key, value) in zip(row.pop("format").split(":"), row.pop("data").split(":"))}
            info_data = {}
            for info_field in row.pop("info").split(";"):
                try:
                    field, value = info_field.split("=")
                except ValueError as e:
                    if info_field and info_field != "IMPRECISE":
                        raise e
                    elif info_field:
                        field, value = ["calling", "IMPRECISE"]
                    else:
                        field, value = ["empty", "empty"]
                info_data[field] = value
            genotype = format_data["GT"].split("/")[0]
            if genotype != "0":
                row["start"] = row.pop("pos")
                if genotype == "1":
                    row["alt"] = "DEL"
                elif genotype == "2":
                    row["alt"] = "DUP"
                end = info_data["END"]
                cnv = {**row, "end": end, "sample_id": sample_id, "format_data": format_data, "info_data": info_data}
                cnvs.append(cnv)
        return cnvs

    @logger.catch(reraise=True)  # TODO: remove decorator once in production
    def process_caller_output(self, vcf_path, sample_id=None):
        cnvs = self.parse_output_file(vcf_path, sample_id)
        gene_bed = self.filter_capture()
        filtered_cnvs = self.filter_cnvs(cnvs, gene_bed)
        if filtered_cnvs:
            logger.debug(f"{self.run_type}\n {filtered_cnvs}\n\n")

    def run_docker_subprocess(self, args, stdout=None, docker_image=None, docker_genome="/mnt/ref_genome/"):
        """Run docker subprocess as root user, mounting input and reference genome dir"""
        ref_genome_dir = os.path.dirname(cnv_pat_settings["genome_fasta_path"])
        if not docker_image:
            docker_image = self.settings["docker_image"]
        docker_command = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{ref_genome_dir}/:{docker_genome}:ro",
            "-v",
            f"{cnv_pat_dir}/input:/mnt/input/:ro",
            "-v",
            f"{self.bam_mount}:/mnt/bam-input/:ro",
            "-v",
            f"{cnv_pat_dir}/cnv-caller-resources/:/mnt/cnv-caller-resources/:ro",
            "-v",
            f"{cnv_pat_dir}/output/:/mnt/output/:rw",
            docker_image,
            *args,
        ]
        logger.log("DOCKER", " ".join(docker_command))
        process = subprocess.run(docker_command, check=True, stdout=stdout)
        return process

    def run_required(self, previous_run_settings_path):
        """Returns True if workflow hasn't been run before or settings have changed since then"""
        if not os.path.exists(previous_run_settings_path):
            logger.info(f"Run required for {self.capture} {self.run_type} {self.gene}")
            self.delete_unused_runs()
            return True
        with open(previous_run_settings_path) as handle:
            previous_settings = toml.load(handle)
            previous_settings.pop("start_time")
            current_settings = dict(self.settings)
            current_settings.pop("start_time")
            if current_settings != previous_settings:
                logger.info(f"Run required for {self.capture} {self.run_type} {self.gene}")
                self.delete_unused_runs()
                return True
        logger.info(f"Re-run not required for {self.capture} {self.run_type} {self.gene}")

    def run_workflow(self):
        """Placeholder for individual tool running"""
        pass

    @logger.catch(reraise=True)
    def main(self):
        previous_run_settings_path = (
            f"{cnv_pat_dir}/successful-run-settings/{self.capture}/{self.run_type}/{self.gene}.toml"
        )
        if self.run_required(previous_run_settings_path):
            self.run_workflow()
            self.write_settings_toml()
        # TODO: read vcfs into database

    def write_settings_toml(self):
        """Write case toml data for successful run"""
        output_dir = f"{cnv_pat_dir}/successful-run-settings/{self.capture}/{self.run_type}/"
        try:
            os.makedirs(output_dir)
        except FileExistsError:
            logger.info(f"run config directory exists for {self.capture}/{self.run_type}")
        output_path = f"{output_dir}/{self.gene}.toml"
        with open(output_path, "w") as out_file:
            toml.dump(self.settings, out_file)
