import glob
import os
import subprocess

from . import utils

cnv_pat_dir = utils.get_cnv_patissier_dir()
genome_fasta_path = "/var/reference_sequences/MiSeq/genome.fa"


class BaseCNVTool:
    def __init__(self, cohort, gene, start_time, normal_panel=True):
        self.start_time = start_time
        self.cohort = cohort
        self.gene = gene
        gene_list = f"{cnv_pat_dir}/input/{cohort}/sample-lists/{gene}_samples.txt"

        sample_ids, bams = utils.SampleUtils.select_samples(gene_list, normal_panel=normal_panel)

        self.bam_mount = utils.SampleUtils.get_mount_point(bams)

        docker_bams = [f"/mnt/bam-input/{bam.split(self.bam_mount)[-1]}" for bam in bams]

        max_cpu = "30"
        max_mem = "50"

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
        }

    def base_output_dirs(self):
        """Returns base directory for output: (system_base, docker_base)"""
        output_base = (
            f"{cnv_pat_dir}/output/{self.cohort}/{self.start_time}/" f"{self.run_type}/{self.gene}"
        )
        docker_output_base = output_base.replace(cnv_pat_dir, "/mnt")

        return (output_base, docker_output_base)

    def run_docker_subprocess(self, args):
        """Run docker subprocess as non root user, mounting input and reference genome dir"""
        ref_genome_dir = os.path.dirname(genome_fasta_path)
        subprocess.run(
            [
                "docker",
                "run",
                # "--user",
                # "1000",
                "-v",
                f"{ref_genome_dir}:/mnt/ref_genome/:ro",
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
        )

    def parse_vcf_4_2(self, vcf_path):
        """Parses VCF v4.2, if positive cnv, returns dict of information"""
        cnvs = []
        with open(vcf_path) as handle:
            for line in vcf_path:
                if line.startswith("#"):
                    continue
                chrom, pos, var_id, ref, alt, qual, var_filter, info, var_format, data = (
                    line.split()
                )

                call_data = {
                    key: value for (key, value) in zip(var_format.split(":"), data.split(":"))
                }
                if call_data["CN"] != "0":
                    call_data["chrom"] = chrom
                    call_data["id"] = var_id
                    call_data["start"] = pos
                    call_data["end"] = info.replace("END=", "")

                    cnvs.append(call_data)
