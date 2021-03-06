"""
Canvas run using targeted capture somatic settings
Version 1.11.0 (Released 3rd June 2016) used because WES workflow broken in later versions: 
    - https://github.com/Illumina/canvas/issues/86 
    - https://github.com/Illumina/canvas/issues/62
Can't reuse normal panel because this is broken and never tested to work, so only using 5 normal samples 
    - https://github.com/Illumina/canvas/issues/96

used filter13 for hg19 from canvas data https://github.com/Illumina/canvas/

This requires GenomeSize.xml in genome dir and kmer.fa genome.fa, though I don't expect to run this 
in large comparison set so don't add to requirements until initial comparison


Have not implemented "chr" prefixed chromosomes at this point as not used
"""

import subprocess
import os

from . import utils, base_classes
from settings import cnv_pat_settings


class Canvas(base_classes.BaseCNVTool):
    def __init__(self, capture, gene, start_time, normal_panel=True):
        self.run_type = "canvas"
        super().__init__(capture, gene, start_time, normal_panel=normal_panel)

        self.settings = {**self.settings, "docker_image": "stefpiatek/canvas:1.11.0"}

    def run_canvas_command(self, args):
        """Creates dir for output and runs a GATK command in docker"""

        base_classes.logger.info(f"Running canvas: {args[0]} \n output: {args[-1]}")
        self.run_docker_subprocess(
            ["mono", "./Canvas.exe", *args],
            docker_genome="/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/",
        )
        base_classes.logger.info(f"Completed canvas: {args[0]} {args[-1]}")

    def run_workflow(self):
        raise Exception("Canvas only gives output for large exome capture only and isn't used in this project")
        try:
            os.makedirs(f"{self.output_base}")
        except FileExistsError:
            pass

        nextera_manifest_out = f"{self.output_base}/nextera_manifest.txt".replace(f"/{self.gene}/", "/")
        docker_nextera_manifest_out = f"{self.docker_output_base}/nextera_manifest.txt".replace(f"/{self.gene}/", "/")
        ploidy_bed_out = f"{self.output_base}/capture_ploidy.bed".replace(f"/{self.gene}/", "/")
        docker_ploidy_bed_out = f"{self.docker_output_base}/capture_ploidy.bed".replace(f"/{self.gene}/", "/")
        parsed_bed_out = f"{self.output_base}/sorted-capture.bed".replace(f"/{self.gene}/", "/")
        docker_parsed_bed_out = f"{self.docker_output_base}/sorted-capture.bed".replace(f"/{self.gene}/", "/")
        capture_vcf_out = f"{self.output_base}/pop_af.vcf".replace(f"/{self.gene}/", "/")
        docker_capture_vcf_out = f"{self.docker_output_base}/pop_af.vcf".replace(f"/{self.gene}/", "/")

        with open(nextera_manifest_out, "w") as output_handle:
            header = [
                "#Nextera Rapid Capture Expanded Exome targeted regions manifest\t\t\t\t\t",
                "[Header]\t\t\t\t",
                "Manifest Version\t1\t\t\t\t",
                "ReferenceGenome\t/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta\t\t\t\t",
                "\t\t\t\t\t",
                "[Regions]\t\t\t\t\t",
                "Name\tChromosome\tStart\tEnd\tUpstream Probe Length\tDownstream Probe Length",
            ]
            for line in header:
                output_handle.write(f"{line}\n")
            with open(f"{base_classes.cnv_pat_dir}/{self.settings['capture_path'][5:]}", "r") as input_handle:
                for line in input_handle:
                    chrom, start, end, name = line.split()
                    output_handle.write("\t".join([name, chrom, start, end, "0", f"0\n"]))

        with open(parsed_bed_out, "w") as handle:
            subprocess.run(
                [
                    # remove 'chr' from bed file chromosome column
                    f"sed -e 's/chr//' {base_classes.cnv_pat_dir}/{self.settings['capture_path'][5:]} | "
                    # sort bed file
                    "sort -k1,1n -k2,2n"
                ],
                check=True,
                shell=True,
                stdout=handle,
            )

        with open(ploidy_bed_out, "w") as ploidy_bed:
            with open(parsed_bed_out, "r") as capture:
                for line in capture:
                    # write in default ploidy of 2, assuming only autosomes
                    ploidy_bed.write(f"chr{line.strip()}\t2\n")

        if not os.path.exists(f"{capture_vcf_out}"):

            with open(capture_vcf_out, "w") as handle:
                base_classes.logger.info("Started downloading population frequencies for canvas")
                proxy = []
                if cnv_pat_settings["http_proxy"]:
                    proxy += ["-e", f"http_proxy={cnv_pat_settings['http_proxy']}"]
                subprocess.run(
                    [
                        "docker",
                        "run",
                        "--rm",
                        "-v",
                        f"{base_classes.cnv_pat_dir}/output/:/mnt/output/:rw",
                        *proxy,
                        "lethalfang/tabix:1.7",
                        "tabix",
                        "http://storage.googleapis.com/gnomad-public/release/2.1/vcf/"
                        "genomes/gnomad.genomes.r2.1.sites.vcf.bgz",
                        "--print-header",
                        "-R",
                        docker_parsed_bed_out,
                    ],
                    check=True,
                    stdout=handle,
                )

            # add chr back on to chromosomes: if line doesn't start with #, add chr to start of line
            subprocess.run(["sed", "-i", "'/^#/! s/^/chr/'", capture_vcf_out])

            # subprocess.run(["gzip", capture_vcf_out], check=True)
            base_classes.logger.info("Finished downloading population frequencies for canvas")

        for index, bam in enumerate(self.settings["unknown_bams"]):
            error_count = 0
            sample = self.bam_to_sample[bam]
            sample_out = f"{self.docker_output_base}/{sample}"
            try:
                os.makedirs(f"{self.output_base}/{sample}")
            except FileExistsError:
                pass
            # Turns out buggy as anything and can't reuse precomputed reference
            # https://github.com/Illumina/canvas/issues/96
            if True:
                # if index == 0:
                # only use 4 samples based on comments in issue 96 above
                control_list = [f"--control-bam={control_bam}" for control_bam in self.settings["normal_bams"]][0:5]
                control_binned_file = f"{sample_out}/TempCNV_{sample}/{sample}.normal.binned"
                bin_size_file = f"{self.output_base}/{sample}/TempCNV_{sample}/{sample}.binsize"
            # else:
            #     with open(bin_size_file, "r") as handle:
            #         bin_size = handle.readline()
            #     control_list = [f"--control-binned={control_binned_file}", f"--control-bin-size={bin_size}"]

            self.run_canvas_command(
                [
                    "Somatic-Enrichment",
                    f"--manifest={docker_nextera_manifest_out}",
                    f"--bam={bam}",
                    *control_list,
                    f"--b-allele-vcf={docker_capture_vcf_out}",
                    "--exclude-non-het-b-allele-sites",  # As not matched from sample
                    f"--sample-name={sample}",  # TODO
                    f"--reference=/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/kmer.fa",
                    f"--genome-folder=/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta",
                    f"--filter-bed=/mnt/cnv-caller-resources/canvas/filter13.bed",
                    f"--ploidy-bed={docker_ploidy_bed_out}",
                    "-s=2",
                    "--output",
                    sample_out,
                ]
            )

            base_classes.logger.info("Manually running checkpoint 3")
            self.run_docker_subprocess(
                [
                    "mono",
                    "CanvasClean.exe",
                    "-i",
                    f"{sample_out}/TempCNV_{sample}/{sample}.ratio.binned",
                    "-o",
                    f"{sample_out}/TempCNV_{sample}/{sample}.cleaned",
                    "-g",
                    "-t",
                    f"{sample_out}/TempCNV_{sample}/manifest.txt",
                ]
            )
            self.run_docker_subprocess(
                [
                    "python3.7",
                    f"/mnt/cnv-caller-resources/canvas/canvas_clean_writer.py",
                    "--output-base",
                    sample_out,
                    "--sample",
                    sample,
                ]
            )
            try:
                self.run_canvas_command(
                    [
                        "Somatic-Enrichment",
                        f"--manifest={docker_nextera_manifest_out}",
                        f"--bam={bam}",
                        *control_list,
                        f"--b-allele-vcf={docker_capture_vcf_out}",
                        "--exclude-non-het-b-allele-sites",  # As not matched from sample
                        f"--sample-name={sample}",
                        f"--reference=/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/kmer.fa",
                        f"--genome-folder=/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta",
                        f"--filter-bed=/mnt/cnv-caller-resources/canvas/filter13.bed",
                        f"--ploidy-bed={docker_ploidy_bed_out}",
                        "-c=4",
                        "-s=4",
                        "--output",
                        sample_out,
                    ]
                )
            except subprocess.CalledProcessError as e:
                # swallow this error as it seems like this step actually finishes sucesfully?
                pass

            self.run_canvas_command(
                [
                    "Somatic-Enrichment",
                    f"--manifest={docker_nextera_manifest_out}",
                    f"--bam={bam}",
                    *control_list,
                    f"--b-allele-vcf={docker_capture_vcf_out}",
                    "--exclude-non-het-b-allele-sites",  # As not matched from sample
                    f"--sample-name={sample}",
                    f"--reference=/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/kmer.fa",
                    f"--genome-folder=/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta",
                    f"--filter-bed=/mnt/cnv-caller-resources/canvas/filter13.bed",
                    f"--ploidy-bed={docker_ploidy_bed_out}",
                    "-c=5",
                    "--output",
                    sample_out,
                ]
            )
