"""
Canvas run using targeted capture somatic settings
Version 1.11.0 (Released 3rd June 2016) used because WES workflow broken in later versions: 
    - https://github.com/Illumina/canvas/issues/86 
    - https://github.com/Illumina/canvas/issues/62

used filter13 for hg19
"""

import subprocess
import os

from . import utils, base_classes

cnv_pat_dir = utils.get_cnv_patissier_dir()

class Canvas(base_classes.BaseCNVTool):
    def __init__(self, cohort, gene, start_time, normal_panel=True):
        super().__init__(cohort, gene, start_time, normal_panel)

        self.run_type = "canvas"

        self.output_base, self.docker_output_base = self.base_output_dirs()

        sample_ids, bams = utils.SampleUtils.select_samples(self.gene_list, normal_panel=False)
        self.bam_mount = utils.SampleUtils.get_mount_point(bams)
        docker_bams = [f"/mnt/bam-input/{bam.split(self.bam_mount)[-1]}" for bam in bams]

        self.settings = {
            **self.settings,
            "contig-ploidy-priors": f"/mnt/cnv-caller-resources/gatk/contig-ploidy-priors.tsv",
            "docker_image": "stefpiatek/canvas:1.11.0",
            "unknown_bams": docker_bams,
        }
        self.settings["normal_bams"] = self.settings.pop("bams")

    def run_canvas_command(self, args):
        """Creates dir for output and runs a GATK command in docker"""

        print(f"*** Running canvas: {args[0]} \n output: {args[-1]} ***")
        self.run_docker_subprocess(
            [
                "mono",
                "./Canvas.exe",
                *args,
            ]
        )
        print(f"*** Completed canvas: {args[0]} {args[-1]} ***")

    def development(self):
        try:
            os.makedirs(f"{self.output_base}")
        except FileExistsError:
            pass

        parsed_bed_out = f"{self.output_base}/sorted-capture.bed".replace(f"/{self.gene}/", "/")
        docker_parsed_bed_out = f"{self.docker_output_base}/sorted-capture.bed".replace(f"/{self.gene}/", "/")
        nextera_manifest_out = f"{self.output_base}/nextera_manifest.txt".replace(f"/{self.gene}/", "/")
        docker_nextera_manifest_out = f"{self.docker_output_base}/nextera_manifest.txt".replace(f"/{self.gene}/", "/")       
        capture_vcf_out = f"{self.output_base}/pop_af.vcf".replace(f"/{self.gene}/", "/")
        docker_capture_vcf_out = f"{self.docker_output_base}/pop_af.vcf".replace(f"/{self.gene}/", "/")

        if not os.path.exists(f"{capture_vcf_out}"):

            with open(nextera_manifest_out, "w") as output_handle:
                header = [
                    "#Nextera Rapid Capture Expanded Exome targeted regions manifest\t\t\t\t\t",
                    "[Header]\t\t\t\t",
                    "Manifest Version\t1\t\t\t\t",
                    "ReferenceGenome\t/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFASTA\t\t\t\t",
                    "\t\t\t\t\t",
                    "[Regions]" "\t\t\t\t\t",
                    "Name\tChromosome\tStart\tEnd\tUpstream Probe Length\tDownstream Probe Length",
                ]
                for line in header:
                    output_handle.write(f"{line}\n")
                with open(f"{cnv_pat_dir}/{self.settings['capture_path'][5:]}", "r") as input_handle:
                    for line in input_handle:
                        chrom, start, end, name = line.split()
                        output_handle.write('\t'.join([name, chrom, start, end, "0", f"0\n"]))

            with open(parsed_bed_out, "w") as handle:
                subprocess.run([
                    # remove 'chr' from bed file chromosome column 
                    f"sed -e 's/chr//' {cnv_pat_dir}/{self.settings['capture_path'][5:]} | "
                    # sort bed file 
                    "sort -k1,1n -k2,2n"  
                    ],
                    check=True,
                    shell=True,
                    stdout=handle
                    )

            with open(capture_vcf_out, "w") as handle:
                print("*** Started downloading population frequencies for canvas ***")
                subprocess.run(
                [
                    "docker",
                    "run",
                    "--rm",
                    "-v",
                    f"{cnv_pat_dir}/output/:/mnt/output/:rw",
                    "-e",
                    "http_proxy=http://10.101.112.70:8080/",  # TODO proxy here
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

            #subprocess.run(["gzip", capture_vcf_out], check=True)
            print("*** Finished downloading population frequencies for canvas ***")

        
        for index, bam in enumerate(self.settings['unknown_bams']):
            error_count = 0
            sample = bam.split('/')[-1].replace('.bam', '').replace('.sorted', '')
            if index == 0:
                control_list = [f"--control-bam={control_bam}" for control_bam in self.settings['normal_bams']][0:4]  # TODO not limit?
                index_bam_file = bam
            else:
                control_list = [f"--control-binned="]  # TODO
            self.run_canvas_command([
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
                "-s=2",
                "--output", 
                self.docker_output_base,
                ])
            try:    
                self.run_canvas_command([
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
                    "-c=3",
                    "-s=3",
                    "--output", 
                    self.docker_output_base,
                    ])
            except subprocess.CalledProcessError:
                print("*** Canvas CheckPoint3 failed, running manually ***")
                with open(f"{self.output_base}/Logging/commands.tsv", "r") as handle:
                    last_command = handle.readlines()[-1]
                # first part of command is the optput file so remove that
                last_command = last_command.replace('"', '').split()[1:]
                self.run_docker_subprocess(last_command)
                self.run_docker_subprocess(
                    ["echo '\"${Output}" f"/TempCNV_{sample}/{sample}.cleaned\"' "
                    f" > {self.docker_output_base}/Checkpoints/03-CanvasClean.json"]
                )

            self.run_canvas_command([
                    "Somatic-Enrichment",
                    f"--manifest={docker_nextera_manifest_out}",
                    f"--bam={bam}",
                    *control_list,
                    f"--b-allele-vcf={docker_capture_vcf_out}",
                    "--exclude-non-het-b-allele-sites",  # As not matched from sample
                    f"--sample-name={bam.split('/')[-1].replace('.bam', '').replace('.sorted', '')}",  # TODO
                    f"--reference=/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/kmer.fa",
                    f"--genome-folder=/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta",
                    f"--filter-bed=/mnt/cnv-caller-resources/canvas/filter13.bed",
                    "-c=4",
                    "--output", 
                    self.docker_output_base,
                    ])


        # TODO: require GenomeSize.xml in genome dir and kmer.fa genome.fa
        # TODO: use sample name from sheet
        # TODO make cohort be capture and then rewrite all of panel or normals to be written above line