import glob
import json
import subprocess
import os

from . import utils

cnv_pat_dir = utils.get_cnv_patissier_dir()
genome_fasta_path = "/var/reference_sequences/MiSeq/genome.fa"

"""
Defaults used from suggested post bu Sheila (gatk staff) on 10th April 2018
with alterations for targetted sequencing
https://gatkforums.broadinstitute.org/gatk/discussion/comment/47431/
"""


class GATKCohort(utils.BaseCNVNormalPanel):
    def __init__(self, cohort, gene, start_time):
        super().__init__(cohort, gene, start_time)
        self.run_type = "gatk-cohort"
        self.output_base, self.docker_output_base = self.base_output_dirs()

        self.settings = {
            **self.settings,
            "num_intervals_per_scatter": "10",  # WES 5000 so for single gene, 10?
            "padding": "50",
            "ref_copy_number_autosomal_contigs": "2",
            "contig-ploidy-priors": (
                f"/mnt/cnv-caller-resources/gatk/contig-ploidy-priors.tsv"
            ),
            "docker_image": "broadinstitute/gatk:4.0.11.0",
            "start_time": start_time,
        }

    def run_gatk_command(self, args):
        """Create dir for output and runs a GATK command in docker"""
        try:
            os.makedirs(f"{self.output_base}/{args[0]}")
        except FileExistsError:
            print(f"*** Folder {args[0]} already exists ***")

        print(f"*** Running  GATK: {args[0]} \n output: {args[-1]} ***")
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
        print(f"*** Completed  GATK: {args[0]} {args[-1]} ***")

    def run_workflow(self):
        pre_process_interval_out = (
            f"{self.docker_output_base}/PreprocessIntervals/intervals.interval_list"
        )
        self.run_gatk_command(
            [
                "PreprocessIntervals",
                "-R",
                self.settings["ref_fasta"],
                "-L",
                self.settings["intervals"],
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
            hdf5_name = bam.split("/")[-1].replace(".bam", ".hdf5")
            collect_read_count_out = (
                pre_process_interval_out.replace(
                    "PreprocessIntervals/", "CollectReadCounts/"
                )
                .replace("intervals.interval_list", hdf5_name)
                .replace(".sorted", "")
            )
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

        # # scattered intervals not using but kept in for a little bit
        # split_intervals_out = f"{self.docker_output_base}/SplitIntervals/"

        # self.run_gatk_command(
        #     [
        #         "SplitIntervals",
        #         "-R",
        #         self.settings["ref_fasta"],
        #         "-L",
        #         pre_process_interval_out,
        #         "--scatter-count",
        #         self.settings["num_intervals_per_scatter"],
        #         "-O",
        #         split_intervals_out,
        #     ]
        # )

        input_flags = []
        for sample in collect_read_counts:
            input_flags.append("--input")
            input_flags.append(sample)

        determine_germline_ploidy_dir = (
            f"{self.docker_output_base}/DetermineGermlineContigPloidy"
        )
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

        post_germline_cnv_caller_dir = (
            f"{self.docker_output_base}/PostprocessGermlineCNVCalls"
        )

        for sample in glob.glob(
            f"{self.output_base}/GermlineCNVCaller/normal-cohort-run-calls/SAMPLE_*"
        ):
            sample_index = sample.split("_")[-1]
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
                    "chrX",
                    "--allosomal-contig",
                    "chrY",
                    "--sequence-dictionary",
                    self.settings["ref_fasta"].replace(".fa", ".dict"),
                    "--output-genotyped-intervals",
                    f"{post_germline_cnv_caller_dir}/sample_{sample_index}_genotyped_intervals.vcf",
                    "--output-genotyped-segments",
                    f"{post_germline_cnv_caller_dir}/sample_{sample_index}_genotyped_segments.vcf",
                ]
            )

    def create_json(self):
        """Write case json data for successful run"""

        output_path = f"{cnv_pat_dir}/run-config-files/gatk-cohort-{self.gene}.json"
        with open(output_path, "w") as out_file:
            json.dump(self.settings, out_file)


class GATKCase:
    def __init__(self, cohort, gene, start_time, sample_bam):
        super(cohort, gene, start_time).__init__()
        # config
        max_cpu = 30
        max_mem = 16
        bam_dir = "/home/chris/stef/icr-mlpa/aligned"
        bed_file = "/home/chris/stef/bed-files/icr-filtered.bed"
        genome_fasta_path = "/var/reference_sequences/MiSeq/genome.fa"
        self.cromwell_jar_path = "/home/chris/stef/cromwell-32.jar"
        self.cohort_wdl_path = (
            "/home/chris/stef/cnv-patissier/workflows/cnv_germline_case_workflow.wdl"
        )

        # cohort job
        cromwell_dir = "/home/chris/stef/cromwell-executions"
        cohort_job_id = "7d5650f7-ed5d-4367-bbaf-cab25b8c1fca"
        model_prefix = "ICR"
        self.sample_name = sample_bam.split("/")[-1].replace(".bam", "")
        self.case_json_data = {
            "CNVGermlineCaseWorkflow.bam": f"{sample_bam}",
            "CNVGermlineCaseWorkflow.bam_idx": f"{sample_bam}.bai",
            "CNVGermlineCaseWorkflow.ref_fasta": genome_fasta_path,
            "CNVGermlineCaseWorkflow.ref_fasta_fai": f"{genome_fasta_path}.fai",
            "CNVGermlineCaseWorkflow.ref_fasta_dict": genome_fasta_path.replace(
                ".fa", ".dict"
            ),
            "CNVGermlineCaseWorkflow.intervals": bed_file,
            "CNVGermlineCaseWorkflow.contig_ploidy_model_tar": (
                f"{cromwell_dir}/CNVGermlineCohortWorkflow/{cohort_job_id}/"
                f"call-DetermineGermlineContigPloidyCohortMode/execution/{model_prefix}-contig-ploidy-model.tar.gz"
            ),
            "CNVGermlineCaseWorkflow.gcnv_model_tars": glob.glob(
                f"{cromwell_dir}/CNVGermlineCohortWorkflow/{cohort_job_id}/"
                f"call-GermlineCNVCallerCohortMode/shard-*/execution/{model_prefix}-model*.tar.gz"
            ),
            "CNVGermlineCaseWorkflow.num_intervals_per_scatter": 5000,
            "CNVGermlineCaseWorkflow.padding": 50,
            "CNVGermlineCaseWorkflow.ref_copy_number_autosomal_contigs": 2,
            "CNVGermlineCaseWorkflow.cpu_for_determine_germline_contig_ploidy": max_cpu,
            "CNVGermlineCaseWorkflow.CollectCounts.cpu": max_cpu,
            "CNVGermlineCaseWorkflow.PreprocessIntervals.cpu": max_cpu,
            "CNVGermlineCaseWorkflow.cpu_for_germline_cnv_caller": max_cpu,
            "CNVGermlineCaseWorkflow.mem_gb_for_germline_cnv_caller": max_mem,
            "CNVGermlineCaseWorkflow.PreprocessIntervals.mem_gb": max_mem,
            "CNVGermlineCaseWorkflow.ScatterIntervals.mem_gb": max_mem,
            "CNVGermlineCaseWorkflow.PostprocessGermlineCNVCalls.mem_gb": max_mem,
            "CNVGermlineCaseWorkflow.gatk_docker": "broadinstitute/gatk:4.0.7.0",
        }

    def remove_successful_run(self):
        """Remove cromwell execution and docker container"""
        pass

    def run_gcnv(self, json_path):
        """Run GATK gCNV CNVGermlineCaseWorkflow"""
        print(f"Running GATK gCNV for {self.sample_name}")
        gcnv = subprocess.run(
            [
                "java",
                f"-Xmx{self.max_mem}g" "-jar",
                self.cromwell_jar_path,
                "run",
                self.cohort_wdl_path,
                "-i",
                json_path,
            ],
            stdout=subprocess.PIPE,
            universal_newlines=True,
            check=True,
        )

        gatk_out_lines = [
            line.strip().rstrip(",")
            for line in gcnv.stdout.splitlines()
            if line.strip().startswith(
                ('"id"', '"CNVGermlineCaseWorkflow.genotyped_segments_vcf"')
            )
        ]
        gatk_run_info = {}
        for line in gatk_out_lines:
            line_dict = "{" + line + "}"
            gatk_run_info.update(eval(line_dict))
        return gatk_run_info

    def save_cnv_calls(self, run_info):
        """Write CNV calls to database of results"""
        pass

    def main(self):
        """Create json, run gCNV and clean up if successful"""
        json_path = self.write_case_json()
        gatk_run_info = self.run_gcnv(json_path)
        gatk_run_info = {
            "CNVGermlineCaseWorkflow.genotyped_segments_vcf": "/home/chris/stef/cnv-patissier/cromwell-executions/CNVGermlineCaseWorkflow/63a2f7e5-44bc-44c0-b23f-a8a76cdc6c02/call-PostprocessGermlineCNVCalls/execution/genotyped-segments-17296.sorted.vcf.gz",
            "id": "63a2f7e5-44bc-44c0-b23f-a8a76cdc6c02",
        }

    def write_case_json(self):
        """Write case json data and returns the path"""
        output_path = f"run-config-files/{self.sample_name}-gCNV.json"
        with open(output_path, "w") as out_file:
            json.dump(self.case_json_data, out_file)
        return output_path


if __name__ == "__main__":
    for bam in glob.glob("/home/chris/stef/icr-mlpa/aligned/*.bam"):
        case_runner = GATKCase(bam)
        case_runner.main()
