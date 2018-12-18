"""
Excavator2 run with defaults https://github.com/matheuscburger/Excavator2/blob/master/docs/EXCAVATOR2%20Manual.pdf

"""

import os
import subprocess

from . import utils, base_classes

cnv_pat_dir = utils.get_cnv_patissier_dir()

class Excavator2(base_classes.BaseCNVTool):
    def __init__(self, cohort, gene, start_time, normal_panel=True):
        super().__init__(cohort, gene, start_time, normal_panel)

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

    def development(self):
        experimental_file_prep = f"{self.output_base}/ExperimentalFilePrepare.txt"
        docker_experimental_file_prep = experimental_file_prep.replace(self.output_base, self.docker_output_base)
        with open(experimental_file_prep, "w") as handle:
            for bam in self.settings['unknown_bams'] + self.settings['normal_bams']:
                sample = self.bam_to_sample[bam]
                handle.write(f"{bam} {self.docker_output_base}/DataPrep/{sample} {sample}\n")

        prepared_bed_file = f"{self.output_base}/capture.bed"
        extra_chroms_bed = f"{self.output_base}/capture_with_mock.bed"
        docker_extra_chroms_bed = f"{self.output_base}/capture_with_mock.bed"
        chromosomes = [f"chr{chrom}" for chrom in list(range(1, 23)) + ["X"]]
        with open(extra_chroms_bed, "w") as bed_output:
            for index, chromosome in enumerate(chromosomes):
               with open(self.settings['capture_path'].replace("/mnt", cnv_pat_dir), "r") as bed_input:
                    if not any([line.startswith(chromosome) for line in bed_input]):
                        bed_output.write("\t".join([chromosome, "10000", "10010", "dummy\n"]))
            with open(self.settings['capture_path'].replace("/mnt", cnv_pat_dir), "r") as bed_input:
                for line in bed_input:
                    bed_output.write(f"{line}")


        docker_prepared_bed_file = prepared_bed_file.replace(self.output_base, self.docker_output_base)
        with(open(prepared_bed_file.replace(".bed", ".sorted"), "w")) as handle:
            subprocess.run(
                f"sort -k1,1 -k2,2n {docker_extra_chroms_bed}",
                shell=True,
                check=True,
                stdout=handle)

        source_target = f"{self.output_base}/SourceTarget.txt"
        docker_source_target = source_target.replace(self.output_base, self.docker_output_base)
        with open(source_target, "w") as handle:
            handle.write(f"/usr/EXCAVATOR2_Package_v1.1.2/ucsc.hg19.bw {self.settings['ref_fasta']}")

        with open(prepared_bed_file, "w") as handle:
            self.run_docker_subprocess(
                ["bedtools", "merge", "-i", f"{docker_prepared_bed_file.replace('.bed', '.sorted')}"],
                docker_image="stefpiatek/bedtools:latest",
                stdout=handle)

        data_file_prep = f"{self.output_base}/ExperimentalFilePrepare.txt"
        docker_data_file_prep = data_file_prep.replace(self.output_base, self.docker_output_base)

        with open(data_file_prep, "w") as handle:
            for bam in self.settings['unknown_bams']:
                sample = self.bam_to_sample[bam]
                handle.write(f"'T' {self.docker_output_base}/{sample} {sample}\n")
        with open(data_file_prep, "a") as handle:
            for bam in self.settings['normal_bams']:
                sample = self.bam_to_sample[bam]
                handle.write(f"'C' {self.docker_output_base}/{sample} {sample}\n")

        self.run_docker_subprocess(["python3.6", f"/mnt/cnv-caller-resources/excavator2/excavator_runner.py", 
            "--output-base", self.docker_output_base,
            "--max-mem", self.settings['max_mem'],
            "--max-cpu", self.settings['max_cpu'],
            ])

    def run_workflow(self):
        try:
            os.makedirs(f"{self.output_base}/")
        except FileExistsError:
            print(f"*** Excavator2 output folder already exists ***")


