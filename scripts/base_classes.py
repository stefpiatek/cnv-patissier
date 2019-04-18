"""
Base class for running of a CNV tool - each tool inherits from this class
"""

import csv
import datetime
import glob
import json
import os
import pathlib
import subprocess
import sys
import time

from sqlalchemy import sql
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
        self.script_dirs = [f"{cnv_pat_dir}/{folder}" for folder in ["scripts", "cnv-caller-resources"]]

        if "PYTEST_CURRENT_TEST" in os.environ.keys():
            self.settings = {"chromosome_prefix": "chr"}
        else:
            # will not be done during pytest running of tests
            self.max_cpu = cnv_pat_settings["max_cpu"]
            self.max_mem = cnv_pat_settings["max_mem"]
            self.sample_sheet = f"{cnv_pat_dir}/input/{capture}/sample-sheets/{gene}.txt"
            self.output_base, self.docker_output_base = self.base_output_dirs()

            normal_sample_ids, normal_bams = utils.SampleUtils.select_samples(self.sample_sheet, normal_panel=True)
            unknown_sample_ids, unknown_bams = utils.SampleUtils.select_samples(self.sample_sheet, normal_panel=False)
            self.bam_mount = utils.SampleUtils.get_mount_point(unknown_bams + normal_bams)
            normal_docker_bams = [f"/mnt/bam-input/{bam.split(self.bam_mount)[-1]}" for bam in normal_bams]
            unknown_docker_bams = [f"/mnt/bam-input/{bam.split(self.bam_mount)[-1]}" for bam in unknown_bams]

            bam_to_sample = utils.SampleUtils.get_bam_to_id(self.sample_sheet)
            self.bam_to_sample = {
                f"/mnt/bam-input/{bam.split(self.bam_mount)[-1]}": sample_id
                for (bam, sample_id) in bam_to_sample.items()
            }
            self.sample_to_bam = {sample_id: bam for (bam, sample_id) in self.bam_to_sample.items()}

            self.settings = {
                "normal_bams": normal_docker_bams,
                "ref_fasta": f"/mnt/ref_genome/{cnv_pat_settings['genome_fasta_path'].split('/')[-1]}",
                "genome_build_name": cnv_pat_settings["genome_build_name"],
                "intervals": f"/mnt/input/{capture}/bed/{gene}.bed",
                "docker_image": None,
                "chromosome_prefix": cnv_pat_settings["chromosome_prefix"],
                "capture": capture,
                "gene": gene,
                "start_time": start_time,
                # To add in after sure all genes run correctly
                # "sample_sheet_md5sum": self.get_md5sum(self.sample_sheet),
                "capture_path": f"/mnt/input/{capture}/bed/{capture}.bed",
                "unknown_bams": unknown_docker_bams,
            }

    def base_output_dirs(self):
        """Returns base directory for output: (system_base, docker_base)"""
        output_base = f"{cnv_pat_dir}/output/{self.capture}/{self.start_time}/{self.run_type}/{self.gene}"
        docker_output_base = output_base.replace(cnv_pat_dir, "/mnt")

        return (output_base, docker_output_base)

    def check_chrom_prefix(self, bed_file):
        """Raises exception if chromosome prefix doesn't match chromsome in bed file"""
        chromosome_names = [x for x in range(1, 23)] + ["X", "Y", "M", "mt"]
        chromosomes = [f"{self.settings['chromosome_prefix']}{chromsome}\t" for chromsome in chromosome_names]
        with open(bed_file, "r") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not any([line.startswith(chrom) for chrom in chromosomes]):
                    raise Exception(
                        "BED file contains line which has an invalid chromosome:\n"
                        f"Line number: {line_number}\n"
                        "Line: '{}'\n".format(line.replace("\t", "<tab>").rstrip()) + 
                        "Expected format: '{}'\n".format(chromosomes[0].replace("\t", "<tab>start<tab>end<tab>gene"))
                        + "Please update 'chromosome_prefix' in local settings file, or alter the BED file."
                    )

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
                        cnv["json_data"] = {field: cnv.pop(field) for field in self.extra_db_fields}
                        filtered_cnvs.append(cnv)
                        break
        return filtered_cnvs

    def get_bam_header(self, sample_id):
        docker_bam = self.sample_to_bam[sample_id]
        samtools = self.run_docker_subprocess(
            ["samtools", "view", "-H", docker_bam], docker_image="stefpiatek/samtools:1.9", stdout=subprocess.PIPE
        )
        header = str(samtools.stdout, "utf-8")

        for line in header.split("\n"):
            for field_value in line.split():
                if field_value.startswith("SM:"):
                    sample_name = field_value.strip("SM:")
        if not sample_name:
            raise Exception(f"No SM tag found in file {docker_bam}")
        assert (
            sample_name == sample_id
        ), f"sample-sheet sample_id '{sample_id}' for '{docker_bam}' does not match bam SM: '{sample_name}'"

        return header

    def get_md5sum(self, file_path):
        md5sum_proc = subprocess.run(["md5sum", file_path], check=True, stdout=subprocess.PIPE)
        md5sum, path = str(md5sum_proc.stdout, "utf-8").split()
        return md5sum, path.replace(cnv_pat_dir, "cnv-pat")

    def get_normal_panel_time(self):
        normal_path = (
            f"{cnv_pat_dir}/successful-run-settings/{self.capture}/"
            f"{self.run_type.replace('case', 'cohort')}/{self.gene}.toml"
        )
        with open(normal_path) as handle:
            normal_config = toml.load(handle)
        return f"{normal_config['start_time']}"

    def get_normal_panel_duration(self):
        normal_path = (
            f"{cnv_pat_dir}/successful-run-settings/{self.capture}/"
            f"{self.run_type.replace('case', 'cohort')}/{self.gene}.toml"
        )
        with open(normal_path) as handle:
            normal_config = toml.load(handle)
        duration = normal_config["end_datetime"] - normal_config["start_datetime"]
        return duration

    def parse_output_file(file_path, sample_id=None):
        """Dummy method to be overwritten by each CNV-caller class"""
        pass

    @staticmethod
    def parse_vcf(input_vcf, sample_id):
        """Parses VCF v4.0 - v4.2, if positive cnv, returns dicts of information within a list"""
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
                        # empty value so just skip field
                        continue
                info_data[field] = value
            if "CN" in format_data.keys():  # copy number is copy info
                copy_info = int(format_data["CN"])
            else:  # genotype gives copy number info as a string
                copy_info = format_data["GT"]
            if copy_info != "0" and copy_info != 2:
                row["start"] = row.pop("pos")
                end = info_data["END"]
                if isinstance(copy_info, int):
                    if copy_info < 2:
                        row["alt"] = "DEL"
                    elif copy_info > 2:
                        row["alt"] = "DUP"
                else:
                    if copy_info == "1":
                        row["alt"] = "DEL"
                    elif copy_info == "2":
                        row["alt"] = "DUP"
                cnv = {**row, "end": end, "sample_id": sample_id, "format_data": format_data, "info_data": info_data}
                cnvs.append(cnv)
        return cnvs

    def prerun_steps(self, sample_sheet_path, ref_genome_path):
        """
        Returns dictionary of sample_id: bam header
        Checks: 
          - filenames have no invalid characters (check_files)
          - file paths exist (check_files)
          - file paths are unique (check_files)
          - sample_ids are unique (check_unique)
          - reference genome files exist
          - SN tag in bam header (from get_bam_header)
          - sample id matches the sample ID given (from get_bam_header)

        """
        bam_headers = {}
        sample_paths = []
        sample_ids = []
        with open(sample_sheet_path, "r") as handle:
            sample_sheet = csv.DictReader(handle, dialect="excel", delimiter="\t")
            for line in sample_sheet:
                sample_paths.append(line["sample_path"])
                sample_ids.append(line["sample_id"])

        utils.SampleUtils.check_files(sample_paths)
        utils.SampleUtils.check_unique(sample_ids, "sample_id")
        for extension in ["", ".fai", ".dict"]:
            ref_genome = pathlib.Path(ref_genome_path + extension)
            assert (
                ref_genome.exists()
            ), f"{ref_genome} does not exist\nPlease edit your settings file or create the file"

        for sample_id, sample_path in zip(sample_ids, sample_paths):
            bam_header = self.get_bam_header(sample_id)
            bam_headers[sample_id] = bam_header
            #  to avoid `docker: Error response from daemon: container did not start before the specified timeout.`
            time.sleep(5)

        return bam_headers

    def process_caller_output(self, sample_path, sample_id=None):
        try:
            cnvs = self.parse_output_file(sample_path, sample_id)
        except FileNotFoundError as e:
            if self.run_type == "excavator2":
                # excavator2 sometimes doesn't produce vcf file
                cnvs = []
            else:
                raise e
        gene_bed = self.filter_capture()
        filtered_cnvs = self.filter_cnvs(cnvs, gene_bed)
        return filtered_cnvs

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
            "-v",
            f"{cnv_pat_dir}/tests/test_files/:/mnt/test_files/:rw",
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
            for field in ["start_time", "start_datetime", "end_datetime"]:
                previous_settings.pop(field)
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

    def upload_all_called_cnvs(self, output_paths, sample_ids):
        for path_and_id in zip(output_paths, sample_ids):
            cnv_calls = self.process_caller_output(path_and_id[0], path_and_id[1])
            for cnv in cnv_calls:
                self.upload_called_cnv(cnv)

    def upload_all_known_data(self):
        self.upload_cnv_caller()
        self.upload_gene()
        known_cnv_table = self.upload_samples(self.sample_sheet)
        self.upload_positive_cnvs(known_cnv_table)

    def upload_all_md5sums(self, run_id):
        for folder in self.script_dirs:
            folder_path = pathlib.Path(folder)
            files = [python for python in folder_path.glob("**/*.py")] + [R for R in folder_path.glob("**/*.R")]
            for file in files:
                md5sum, file_path = self.get_md5sum(file)
                Queries.update_or_create(
                    models.File, self.session, defaults={"run_id": run_id, "relative_path": file_path}, md5sum=md5sum
                )

    def upload_cnv_caller(self):
        Queries.get_or_create(models.Caller, self.session, defaults=dict(name=self.run_type))

    def upload_gene(self):
        bed = self.filter_capture()
        defaults = {"name": self.gene, "genome_build": self.settings["genome_build_name"], "capture": self.capture}
        upload_data = {"chrom": bed[0].split()[0], "start": int(bed[0].split()[1]), "end": int(bed[-1].split()[2])}

        Queries.update_or_create(models.Gene, self.session, defaults=defaults, **upload_data)
        self.session.commit()

    def upload_positive_cnvs(self, known_cnv_table):
        for row in known_cnv_table:
            Queries.get_or_create(models.CNV, self.session, defaults=row["cnv"])
            sample_instance = self.session.query(models.Sample).filter_by(**row["sample_defaults"]).first()
            cnv_instance = self.session.query(models.CNV).filter_by(**row["cnv"]).first()
            known_cnv = dict(cnv_id=cnv_instance.id, sample_id=sample_instance.id)
            Queries.get_or_create(models.KnownCNV, self.session, defaults=known_cnv)
            self.session.commit()

    def upload_samples(self, sample_sheet_path):
        with open(sample_sheet_path, "r") as handle:
            known_cnv_table = []
            gene_instance = (
                self.session.query(models.Gene)
                .filter_by(name=self.gene, capture=self.capture, genome_build=self.settings["genome_build_name"])
                .first()
            )

            sample_sheet = csv.DictReader(handle, dialect="excel", delimiter="\t")
            for line in sample_sheet:
                bam_header = self.bam_headers[line["sample_id"]]
                sample_defaults = {"name": line["sample_id"], "path": line["sample_path"], "gene_id": gene_instance.id}
                sample_data = {"bam_header": bam_header, "result_type": line["result_type"]}

                Queries.update_or_create(models.Sample, self.session, defaults=sample_defaults, **sample_data)

                if line["result_type"] == "positive":
                    cnv = {
                        "alt": line["cnv_call"],
                        "genome_build": self.settings["genome_build_name"],
                        "chrom": line["chromosome"],
                        "start": line["start"],
                        "end": line["end"],
                    }
                    known_cnv_table.append({"cnv": cnv, "sample_defaults": sample_defaults})
                self.session.commit()

        return known_cnv_table

    def upload_called_cnv(self, cnv_call):
        json_data = json.dumps(cnv_call.pop("json_data"))
        sample_name = cnv_call.pop("sample_id")
        cnv_call.update({"genome_build": self.settings["genome_build_name"]})

        Queries.get_or_create(models.CNV, self.session, defaults=cnv_call)
        caller_instance = self.session.query(models.Caller).filter_by(name=self.run_type).first()
        cnv_instance = self.session.query(models.CNV).filter_by(**cnv_call).first()
        gene_instance = (
            self.session.query(models.Gene)
            .filter_by(name=self.gene, capture=self.capture, genome_build=self.settings["genome_build_name"])
            .first()
        )
        sample_instance = (
            self.session.query(models.Sample).filter_by(name=sample_name, gene_id=gene_instance.id).first()
        )

        called_cnv_defaults = dict(caller_id=caller_instance.id, cnv_id=cnv_instance.id, sample_id=sample_instance.id)
        Queries.update_or_create(models.CalledCNV, self.session, defaults=called_cnv_defaults, json_data=json_data)
        self.session.commit()

    def upload_run_data(self, sample_names):
        gene_instance = (
            self.session.query(models.Gene)
            .filter_by(name=self.gene, capture=self.capture, genome_build=self.settings["genome_build_name"])
            .first()
        )
        caller_instance = self.session.query(models.Caller).filter_by(name=self.run_type).first()
        sample_ids = []
        for sample_name in sample_names:
            sample_instance = (
                self.session.query(models.Sample).filter_by(name=sample_name, gene_id=gene_instance.id).first()
            )
            sample_ids.append(sample_instance.id)

        duration = self.settings["end_datetime"] - self.settings["start_datetime"]
        if self.run_type.endswith("case"):
            duration += self.get_normal_panel_duration()

        run_defaults = {"gene_id": gene_instance.id, "caller_id": caller_instance.id}
        upload_data = {"samples": json.dumps(sample_ids), "duration": duration}
        run_instance, created = Queries.update_or_create(models.Run, self.session, defaults=run_defaults, **upload_data)
        # self.upload_all_md5sums(run_instance.id)
        self.session.commit()

    @logger.catch(reraise=True)
    def main(self):
        previous_run_settings_path = (
            f"{cnv_pat_dir}/successful-run-settings/{self.capture}/{self.run_type}/{self.gene}.toml"
        )
        if self.run_required(previous_run_settings_path):
            if self.run_type.endswith("cohort"):
                self.bam_headers = self.prerun_steps(self.sample_sheet, cnv_pat_settings["genome_fasta_path"])
                self.settings["start_datetime"] = datetime.datetime.now()
                self.run_workflow()
                self.settings["end_datetime"] = datetime.datetime.now()
            else:
                self.bam_headers = self.prerun_steps(self.sample_sheet, cnv_pat_settings["genome_fasta_path"])
                self.settings["start_datetime"] = datetime.datetime.now()
                output_paths, sample_ids = self.run_workflow()
                self.settings["end_datetime"] = datetime.datetime.now()
                self.upload_all_known_data()
                self.upload_all_called_cnvs(output_paths, sample_ids)
                self.upload_run_data(sample_ids)
            self.write_settings_toml()

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


class Queries:
    @staticmethod
    def get_or_create(model, session, defaults=None, **kwargs):
        params = {k: v for k, v in kwargs.items() if not isinstance(v, sql.ClauseElement)}
        params.update(defaults or {})
        instance = session.query(model).filter_by(**params).first()
        if instance:
            return instance, False
        else:
            instance = model(**params)
            session.add(instance)
            return instance, True

    @staticmethod
    def update_or_create(model, session, defaults=None, **kwargs):
        params = {k: v for k, v in kwargs.items() if not isinstance(v, sql.ClauseElement)}
        params.update(defaults or {})
        query = session.query(model).filter_by(**defaults)
        if query.first():
            query.update(kwargs)
            instance = query.first()
            return instance, False
        else:
            instance = model(**params)
            session.add(instance)
            return instance, True
