import glob
import subprocess
import csv
import os

genome_fasta_path = "/var/reference_sequences/MiSeq/genome.fa"


class BaseCNVNormalPanel:
    def __init__(self, cohort, gene, start_time):
        self.start_time = start_time
        self.cohort = cohort
        self.gene = gene
        gene_list = f"{cnv_pat_dir}/input/{cohort}/sample-lists/{gene}_samples.txt"

        bams = SampleUtils.select_samples(gene_list, "normal-panel")
        self.bam_mount = SampleUtils.get_mount_point(bams)

        docker_bams = [
            f"/mnt/bam-input/{bam.split(self.bam_mount )[-1]}" for bam in bams
        ]

        max_cpu = "30"
        max_mem = "50"

        self.settings = {
            "bams": docker_bams,
            "ref_fasta": f"/mnt/ref_genome/{genome_fasta_path.split('/')[-1]}",
            "intervals": f"/mnt/input/{cohort}/bed/{gene}.bed",
            "max_cpu": max_cpu,
            "max_mem": max_mem,
            "docker_image": None,
        }

    def base_output_dirs(self):
        output_base = (
            f"{cnv_pat_dir}/output/{self.cohort}/{self.start_time}/"
            f"{self.run_type}/{self.gene}"
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


def get_cnv_patissier_dir():
    """Returns the base directory of the project from the scripts folder"""
    return os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


class SampleUtils:
    @classmethod
    def select_samples(cls, gene_list, cnv_status):
        """returns list of sample paths for normal, normal-panel and postive samples"""
        output_samples = []
        with open(gene_list) as handle:
            samples = csv.DictReader(handle, delimiter="\t")
            for sample in samples:
                if sample["result_type"] == cnv_status:
                    output_samples.append(sample["sample_path"])
        return output_samples

    @classmethod
    def get_mount_point(cls, paths):
        """Returns common directory for a set of paths"""
        common_prefix = os.path.commonprefix(paths)
        end_of_path = common_prefix.rfind("/")
        common_path = common_prefix[: end_of_path + 1]
        return common_path


cnv_pat_dir = get_cnv_patissier_dir()
