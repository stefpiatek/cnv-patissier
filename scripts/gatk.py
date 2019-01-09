"""
Defaults used from suggested post by Sheila (gatk staff) on 10th April 2018
with alterations for targetted sequencing: padding 50bp and no interval splitting
https://gatkforums.broadinstitute.org/gatk/discussion/comment/47431/
"""

import glob
import subprocess
import os

import toml

from . import utils, base_classes


class GATKBase(base_classes.BaseCNVTool):
    def __init__(self, capture, gene, start_time, normal_panel):
        super().__init__(capture, gene, start_time, normal_panel)

        self.settings = {
            **self.settings,
            "num_intervals_per_scatter": "10",  # currently not used
            "padding": "50",
            "ref_copy_number_autosomal_contigs": "2",
            "contig-ploidy-priors": f"/mnt/cnv-caller-resources/gatk/contig-ploidy-priors.tsv",
            "docker_image": "broadinstitute/gatk:4.0.11.0",
        }

    def run_gatk_command(self, args):
        """Create dir for output and runs a GATK command in docker"""
        try:
            os.makedirs(f"{self.output_base}/{args[0]}")
        except FileExistsError:
            base_classes.logger.info(f"Folder {args[0]} already exists")

        base_classes.logger.info(f"Running  GATK: {args[0]} \n output: {args[-1]}")
        self.run_docker_subprocess(
            [
                "java",
                f"-Xmx{self.settings['max_mem']}g",
                f"-XX:ConcGCThreads={self.settings['max_cpu']}",
                "-jar",
                "gatk.jar",
                *args,
            ]
        )
        base_classes.logger.info(f"Completed  GATK: {args[0]} {args[-1]}")


class GATKCase(GATKBase):
    def __init__(self, capture, gene, start_time):
        super().__init__(capture, gene, start_time, normal_panel=False)
        self.run_type = "gatk_case"
        self.output_base, self.docker_output_base = self.base_output_dirs()

        normal_path = f"{base_classes.cnv_pat_dir}/successful-run-settings/{self.capture}/gatk_cohort/{self.gene}.toml"
        with open(normal_path) as handle:
            normal_config = toml.load(handle)
        self.normal_path_base = f"/mnt/output/{self.capture}/{normal_config['start_time']}/gatk_cohort/{self.gene}"

        self.settings = {
            **self.settings,
            "normal_panel_start_time": normal_config["start_time"],
            "pre_process_interval": (f"{self.normal_path_base}/PreprocessIntervals/intervals.interval_list"),
            "contig_ploidy_model": (f"{self.normal_path_base}/DetermineGermlineContigPloidy/normal-cohort-model"),
            "gcnv_model": (f"{self.normal_path_base}/GermlineCNVCaller/normal-cohort-run-model"),
        }

    def parse_output_file(self, file_path, sample_id):
        with open(file_path) as handle:
            cnvs = self.parse_vcf_4_2(handle, sample_id)
        return cnvs

    def run_workflow(self):
        collect_read_counts = []
        for bam in self.settings["bams"]:
            sample_name = self.bam_to_sample[bam]
            hdf5_name = f"{self.docker_output_base}/CollectReadCounts/{sample_name}.hdf5"
            collect_read_counts.append(hdf5_name)
            self.run_gatk_command(
                [
                    "CollectReadCounts",
                    "-I",
                    bam,
                    "-L",
                    self.settings["pre_process_interval"],
                    "--interval-merging-rule",
                    "OVERLAPPING_ONLY",
                    "-O",
                    hdf5_name,
                ]
            )

        determine_germline_ploidy_dir = f"{self.docker_output_base}/DetermineGermlineContigPloidy"

        input_flags = []
        for sample in collect_read_counts:
            input_flags.append("--input")
            input_flags.append(sample)

        self.run_gatk_command(
            [
                "DetermineGermlineContigPloidy",
                *input_flags,
                "--model",
                self.settings["contig_ploidy_model"],
                "--output-prefix",
                "case",
                "--output",
                determine_germline_ploidy_dir,
            ]
        )

        germline_cnv_caller_dir = f"{self.docker_output_base}/GermlineCNVCaller"
        self.run_gatk_command(
            [
                "GermlineCNVCaller",
                "--run-mode",
                "CASE",
                "--contig-ploidy-calls",
                f"{determine_germline_ploidy_dir}/case-calls",
                "--model",
                self.settings["gcnv_model"],
                *input_flags,
                "--output-prefix",
                "case-run",
                "--output",
                germline_cnv_caller_dir,
            ]
        )

        post_germline_cnv_caller_dir = f"{self.docker_output_base}/PostprocessGermlineCNVCalls"

        sample_names = []

        for sample in glob.glob(f"{self.output_base}/GermlineCNVCaller/case-run-calls/SAMPLE_*"):
            sample_index = sample.split("_")[-1]
            with open(f"{sample}/sample_name.txt") as handle:
                sample_name = handle.readline().strip()
                sample_names.append(sample_name)
            self.run_gatk_command(
                [
                    "PostprocessGermlineCNVCalls",
                    "--calls-shard-path",
                    f"{germline_cnv_caller_dir}/case-run-calls",
                    "--model-shard-path",
                    f"{self.settings['gcnv_model']}",
                    "--contig-ploidy-calls",
                    f"{determine_germline_ploidy_dir}/case-calls",
                    "--sample-index",
                    sample_index,
                    "--autosomal-ref-copy-number",
                    self.settings["ref_copy_number_autosomal_contigs"],
                    "--allosomal-contig",
                    f"{self.settings['chromosome_prefix']}X",
                    "--allosomal-contig",
                    f"{self.settings['chromosome_prefix']}Y",
                    "--sequence-dictionary",
                    self.settings["ref_fasta"].replace(".fa", ".dict"),
                    "--output-genotyped-intervals",
                    f"{post_germline_cnv_caller_dir}/{sample_name}_intervals.vcf",
                    "--output-genotyped-segments",
                    f"{post_germline_cnv_caller_dir}/{sample_name}_segments.vcf",
                ]
            )


class GATKCohort(GATKBase):
    def __init__(self, cohort, gene, start_time):
        super().__init__(cohort, gene, start_time, normal_panel=True)
        self.run_type = "gatk_cohort"
        self.output_base, self.docker_output_base = self.base_output_dirs()

    def run_workflow(self):
        pre_process_interval_out = f"{self.docker_output_base}/PreprocessIntervals/intervals.interval_list"
        self.run_gatk_command(
            [
                "PreprocessIntervals",
                "-R",
                self.settings["ref_fasta"],
                "-L",
                self.settings["capture_path"],
                "--padding",
                self.settings["padding"],
                "--bin-length",
                "0",
                "--interval-merging-rule",
                "OVERLAPPING_ONLY",
                "-O",
                pre_process_interval_out,
            ]
        )

        collect_read_counts = []
        for bam in self.settings["bams"]:
            sample_name = self.bam_to_sample[bam]
            collect_read_count_out = f"{self.docker_output_base}/CollectReadCounts/{sample_name}.hd5f"

            collect_read_counts.append(collect_read_count_out)
            self.run_gatk_command(
                [
                    "CollectReadCounts",
                    "-I",
                    bam,
                    "-L",
                    pre_process_interval_out,
                    "--interval-merging-rule",
                    "OVERLAPPING_ONLY",
                    "-O",
                    collect_read_count_out,
                ]
            )

        input_flags = []
        for sample in collect_read_counts:
            input_flags.append("--input")
            input_flags.append(sample)

        determine_germline_ploidy_dir = f"{self.docker_output_base}/DetermineGermlineContigPloidy"
        self.run_gatk_command(
            [
                "DetermineGermlineContigPloidy",
                *input_flags,
                "--contig-ploidy-priors",
                self.settings["contig-ploidy-priors"],
                "--output-prefix",
                "normal-cohort",
                "--output",
                determine_germline_ploidy_dir,
            ]
        )

        germline_cnv_caller_dir = f"{self.docker_output_base}/GermlineCNVCaller"
        self.run_gatk_command(
            [
                "GermlineCNVCaller",
                "--run-mode",
                "COHORT",
                "-L",
                pre_process_interval_out,
                "--interval-merging-rule",
                "OVERLAPPING_ONLY",
                "--contig-ploidy-calls",
                f"{determine_germline_ploidy_dir}/normal-cohort-calls",
                *input_flags,
                "--output-prefix",
                "normal-cohort-run",
                "--output",
                germline_cnv_caller_dir,
            ]
        )

        post_germline_cnv_caller_dir = f"{self.docker_output_base}/PostprocessGermlineCNVCalls"

        for sample in glob.glob(f"{self.output_base}/GermlineCNVCaller/normal-cohort-run-calls/SAMPLE_*"):
            sample_index = sample.split("_")[-1]
            with open(f"{sample}/sample_name.txt") as handle:
                sample_name = handle.readline().strip()
            self.run_gatk_command(
                [
                    "PostprocessGermlineCNVCalls",
                    "--calls-shard-path",
                    f"{germline_cnv_caller_dir}/normal-cohort-run-calls",
                    "--model-shard-path",
                    f"{germline_cnv_caller_dir}/normal-cohort-run-model",
                    "--contig-ploidy-calls",
                    f"{determine_germline_ploidy_dir}/normal-cohort-calls",
                    "--sample-index",
                    sample_index,
                    "--autosomal-ref-copy-number",
                    self.settings["ref_copy_number_autosomal_contigs"],
                    "--allosomal-contig",
                    f"{self.settings['chromosome_prefix']}X",
                    "--allosomal-contig",
                    f"{self.settings['chromosome_prefix']}Y",
                    "--sequence-dictionary",
                    self.settings["ref_fasta"].replace(".fa", ".dict"),
                    "--output-genotyped-intervals",
                    f"{post_germline_cnv_caller_dir}/{sample_name}_intervals.vcf",
                    "--output-genotyped-segments",
                    f"{post_germline_cnv_caller_dir}/{sample_name}_segments.vcf",
                ]
            )
