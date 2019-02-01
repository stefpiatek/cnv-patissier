"""
WISExome 

Atrribution as per  Attribution-NonCommercial-ShareAlike CC BY-NC-SA License:
Copyright (C) 2017 VU University Medical Center Amsterdam
Author: Roy Straver (github.com/rstraver)
https://github.com/VUmcCGP/wisexome/blob/master/LICENSE.md

WISExome not implemented because of bugs in code
"""

import csv
import subprocess
import pathlib

from . import base_classes


class WISExome(base_classes.BaseCNVTool):
    def __init__(self, capture, gene, start_time, normal_panel=True):
        self.run_type = "wisexome"
        super().__init__(capture, gene, start_time, normal_panel=normal_panel)
        self.extra_db_fields = []
        self.settings = {**self.settings, "docker_image": "stefpiatek/wisexome:latest"}

    def parse_output_file(self, file_path, sample_id):
        cnvs = []
        return cnvs

    def run_command(self, args, stdout=None):
        pathlib.Path(self.output_base, args[0].replace(".py", "")).mkdir(exist_ok=True, parents=True)
        base_classes.logger.info(f"Running {args[0]} step of {self.run_type} for {self.gene}\n output:{args[-1]}")
        self.run_docker_subprocess(["python2.7", *args], stdout=stdout)
        base_classes.logger.info(f"Completed {args[0]} step of {self.run_type} for {self.gene}\n output:{args[-1]}")

    def run_workflow(self):
        raise Exception("WISExome not implemented because of bugs in code")

        # sorted bed required
        extra_chroms_bed = f"{self.output_base}/capture_with_mock.bed"
        docker_extra_chroms_bed = f"{self.docker_output_base}/capture_with_mock.bed"
        all_chromosomes = [f"chr{chrom}" for chrom in list(range(1, 23)) + ["X", "Y", "M"]]

        pathlib.Path(self.output_base).mkdir(exist_ok=True, parents=True)
        with open(extra_chroms_bed, "w") as bed_output:
            for index, chromosome in enumerate(all_chromosomes):
                with open(self.settings["capture_path"].replace("/mnt", base_classes.cnv_pat_dir), "r") as bed_input:
                    if not any([line.startswith(chromosome) for line in bed_input]):
                        bed_output.write("\t".join([chromosome, "1000000", "1100000", "dummy\n", chromosome, "2000000", "2100000", "dummy\n"]))
            with open(self.settings["capture_path"].replace("/mnt", base_classes.cnv_pat_dir), "r") as bed_input:
                for line in bed_input:
                    bed_output.write(f"{line}")

        for data_type in ["unknown", "normal"]:
            local_hit_dir = pathlib.Path(self.output_base, "consam", data_type)
            local_hit_dir.mkdir(exist_ok=True, parents=True)
            hit_dir = pathlib.Path(self.docker_output_base, "consam", data_type)            
            for bam_file in self.settings[f"{data_type}_bams"]:
                sample_name = self.bam_to_sample[bam_file]
                self.run_command(["consam.py", bam_file, docker_extra_chroms_bed, 
                    f"{hit_dir / sample_name}.hits"])

            local_normalised_hits = pathlib.Path(self.output_base, "lennormalize", data_type)
            local_normalised_hits.mkdir(exist_ok=True, parents=True)
            normalised_hits = pathlib.Path(self.docker_output_base, "lennormalize", data_type)

            self.run_command(["lennormalize.py", f"{hit_dir}", docker_extra_chroms_bed, f"{normalised_hits}"])
      
        # This step doesn't work if a chromosome doesn't have a target on it, 
        autosomes = list(range(1, 23))
        x_and_autosomes = autosomes + ["X"]

        pathlib.Path(f"{self.output_base}/prepref/").mkdir(exist_ok=True, parents=True)
        reference_data_dir = pathlib.Path(f"{self.docker_output_base}/prepref")
        for target in autosomes[::-1]:
            for ref in autosomes:
                if target != ref:
                    self.run_command(["prepref.py",  
                        f"{self.docker_output_base}/lennormalize/normal/", 
                        f"{reference_data_dir}/{target}.{ref}.ref",
                        f"chr{target}", f"chr{ref}"])

        ref_base = f"{self.docker_output_base}/takeref/refname"
        # step fails here because X chromosome isn't in reference set (prepref step)
        # step above fails if you put X chromosome
        # beyond scope of this project to fix the code for this project
        self.run_command(["takeref.py", f"{self.docker_output_base}/prepref", ref_base])

