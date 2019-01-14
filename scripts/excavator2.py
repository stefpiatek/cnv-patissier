"""
Excavator2 run with defaults https://github.com/matheuscburger/Excavator2/blob/master/docs/EXCAVATOR2%20Manual.pdf

"""

import os
import subprocess

from . import utils, base_classes


class Excavator2(base_classes.BaseCNVTool):
    def __init__(self, capture, gene, start_time, normal_panel=True):
        super().__init__(capture, gene, start_time, normal_panel)
        self.extra_db_fields = ["id", "ref", "qual", "filter", "format_data", "info_data"]

        self.run_type = "excavator2"

        self.output_base, self.docker_output_base = self.base_output_dirs()

        sample_ids, bams = utils.SampleUtils.select_samples(self.gene_list, normal_panel=False)
        self.bam_mount = utils.SampleUtils.get_mount_point(bams)
        docker_bams = [f"/mnt/bam-input/{bam.split(self.bam_mount)[-1]}" for bam in bams]

        self.settings = {
            **self.settings,
            "contig-ploidy-priors": f"/mnt/cnv-caller-resources/gatk/contig-ploidy-priors.tsv",
            "docker_image": "stefpiatek/excavator2:1.1.2",
            "unknown_bams": docker_bams,
        }
        self.settings["normal_bams"] = self.settings.pop("bams")

    def parse_output_file(self, file_path, sample_id):
        with open(file_path) as handle:
            cnvs = self.parse_vcf_4_2(handle, sample_id)
        return cnvs

    def run_workflow(self):
        try:
            os.makedirs(f"{self.output_base}/")
        except FileExistsError:
            base_classes.logger.info(f"Excavator2 output folder already exists")
        prepared_bed_file = f"{self.output_base}/capture.bed"
        extra_chroms_bed = f"{self.output_base}/capture_with_mock.bed"
        docker_extra_chroms_bed = f"{self.output_base}/capture_with_mock.bed"
        chromosomes = [f"chr{chrom}" for chrom in list(range(1, 23)) + ["X"]]
        with open(extra_chroms_bed, "w") as bed_output:
            for index, chromosome in enumerate(chromosomes):
                with open(self.settings["capture_path"].replace("/mnt", base_classes.cnv_pat_dir), "r") as bed_input:
                    if not any([line.startswith(chromosome) for line in bed_input]):
                        bed_output.write("\t".join([chromosome, "10000", "10010", "dummy\n"]))
            with open(self.settings["capture_path"].replace("/mnt", base_classes.cnv_pat_dir), "r") as bed_input:
                for line in bed_input:
                    bed_output.write(f"{line}")

        docker_prepared_bed_file = prepared_bed_file.replace(self.output_base, self.docker_output_base)
        with (open(prepared_bed_file.replace(".bed", ".sorted"), "w")) as handle:
            subprocess.run(f"sort -k1,1 -k2,2n {docker_extra_chroms_bed}", shell=True, check=True, stdout=handle)

        source_target = f"{self.output_base}/SourceTarget.txt"
        with open(source_target, "w") as handle:
            handle.write(f"/usr/EXCAVATOR2_Package_v1.1.2/data/ucsc.hg19.bw {self.settings['ref_fasta']}")

        with open(prepared_bed_file, "w") as handle:
            self.run_docker_subprocess(
                ["bedtools", "merge", "-i", f"{docker_prepared_bed_file.replace('.bed', '.sorted')}"],
                docker_image="stefpiatek/bedtools:latest",
                stdout=handle,
            )

        experimental_file_prep = f"{self.output_base}/ExperimentalFilePrepare.txt"
        docker_experimental_file_prep = experimental_file_prep.replace(self.output_base, self.docker_output_base)
        with open(experimental_file_prep, "w") as handle:
            for bam in self.settings["unknown_bams"] + self.settings["normal_bams"]:
                sample = self.bam_to_sample[bam]
                handle.write(f"{bam} {self.docker_output_base}/DataPrep/{sample} {sample}\n")

        exp_file_data = f"{self.output_base}/ExperimentalFileData.txt"
        with open(exp_file_data, "w") as handle:
            t_counter = 0
            for bam in self.settings["unknown_bams"]:
                t_counter += 1
                sample = self.bam_to_sample[bam]
                handle.write(f"T{t_counter} {self.docker_output_base}/DataPrep/{sample} {sample}\n")
        with open(exp_file_data, "a") as handle:
            c_counter = 0
            for bam in self.settings["normal_bams"]:
                c_counter += 1
                sample = self.bam_to_sample[bam]
                handle.write(f"C{c_counter} {self.docker_output_base}/DataPrep/{sample} {sample}\n")

        self.run_docker_subprocess(
            [
                "python3.6",
                f"/mnt/cnv-caller-resources/excavator2/excavator_runner.py",
                "--output-base",
                self.docker_output_base,
                "--max-mem",
                self.settings["max_mem"],
                "--max-cpu",
                self.settings["max_cpu"],
            ]
        )
