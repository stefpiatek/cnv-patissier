"""
Defaults used from XHMM documentation https://atgu.mgh.harvard.edu/xhmm/tutorial.shtml

- masks regions with extreme GC content 
- no parallism at the moment
- no repeat masking using PLINK/Seq
- Altered filtering to allow for high read depth
   - Max target read depth: 1000
   - Max sample mean read depth: 1000
   - Max sample SD read depth: 300
"""

import subprocess
import os

from . import utils, base_classes

cnv_pat_dir = utils.get_cnv_patissier_dir()


class XHMM(base_classes.BaseCNVTool):
    def __init__(self, capture, gene, start_time, normal_panel=True):
        super().__init__(capture, gene, start_time, normal_panel)

        self.run_type = "xhmm"

        self.output_base, self.docker_output_base = self.base_output_dirs()

        sample_ids, bams = utils.SampleUtils.select_samples(self.gene_list, normal_panel=False)
        self.bam_mount = utils.SampleUtils.get_mount_point(bams)
        docker_bams = [f"/mnt/bam-input/{bam.split(self.bam_mount)[-1]}" for bam in bams]

        self.settings = {
            **self.settings,
            "contig-ploidy-priors": f"/mnt/cnv-caller-resources/gatk/contig-ploidy-priors.tsv",
            "docker_image": "stefpiatek/xhmm:1.0",
            "unknown_bams": docker_bams,
        }
        self.settings["normal_bams"] = self.settings.pop("bams")

    def run_gatk_command(self, args):
        """Create dir for output and runs a GATK command in docker"""
        try:
            os.makedirs(f"{self.output_base}/{args[0]}")
        except FileExistsError:
            print(f"*** Folder {args[0]} already exists ***")

        print(f"*** Running  GATK step of XHMM: {args[0]} \n output: {args[-1]} ***")
        self.run_docker_subprocess(
            [
                "java",
                f"-Xmx{self.settings['max_mem']}g",
                f"-XX:ConcGCThreads={self.settings['max_cpu']}",
                "-jar",
                "GenomeAnalysisTK.jar",
                "-T",
                *args,
            ]
        )
        print(f"*** Completed  GATK step of XHMM: {args[0]} {args[-1]} ***")

    def run_xhmm_command(self, args):
        """Runs xhmm command in docker"""
        print(f"*** Running  XHMM: {args[0]} \n output: {args[-1]} ***")
        self.run_docker_subprocess(["xhmm", *args])
        print(f"*** Completed  XHMM: {args[0]} {args[-1]} ***")

    def run_workflow(self):
        os.makedirs(f"{self.output_base}")

        all_bams = self.settings["unknown_bams"] + self.settings["normal_bams"]
        with open(f"{self.output_base}/bam.list", "w") as handle:
            for bam in all_bams:
                handle.write(f"{bam}\n")

        depth_out = f"{self.docker_output_base}/DepthOfCoverage/coverage.DATA"

        self.run_gatk_command(
            [
                "DepthOfCoverage",
                "-I",
                f"{self.docker_output_base}/bam.list",
                "-L",
                self.settings["capture_path"],
                "-R",
                self.settings["ref_fasta"],
                "-dt",
                "BY_SAMPLE",
                "-dcov",
                "5000",
                "-l",
                "INFO",
                "--omitDepthOutputAtEachBase",
                "--omitLocusTable",
                "--minBaseQuality",
                "0",
                "--minMappingQuality",
                "20",
                "--start",
                "1",
                "--stop",
                "5000",
                "--nBins",
                "200",
                "--includeRefNSites",
                "--countType",
                "COUNT_FRAGMENTS",
                "-o",
                depth_out,
            ]
        )

        xhmm_read_depth_out = f"{self.docker_output_base}/DATA.RD.txt"
        self.run_xhmm_command(
            ["--mergeGATKdepths", "--GATKdepths", f"{depth_out}.sample_interval_summary", "-o", xhmm_read_depth_out]
        )

        gc_content_out = f"{self.docker_output_base}/GCContentByInterval/locus_GC.txt"
        self.run_gatk_command(
            [
                "GCContentByInterval",
                "-L",
                self.settings["capture_path"],
                "-R",
                self.settings["ref_fasta"],
                "-o",
                gc_content_out,
            ]
        )

        extreme_gc_out = f"{self.docker_output_base}/extreme_gc_targets.txt"
        with open(extreme_gc_out.replace(self.docker_output_base, self.output_base), "w") as handle:
            subprocess.run(
                [
                    f"cat {gc_content_out.replace(self.docker_output_base, self.output_base)} | "
                    "awk '{if ($2 < 0.1 || $2 > 0.9) print $1}' "
                ],
                shell=True,
                check=True,
                stdout=handle,
            )

        centered_filtered_out = f"{self.docker_output_base}/DATA.filtered_centered.RD.txt"

        self.run_xhmm_command(
            [
                "--matrix",
                "-r",
                xhmm_read_depth_out,
                "--centerData",
                "--centerType",
                "target",
                "--outputExcludedTargets",
                f"{centered_filtered_out}.filtered-targets.txt",
                "--outputExcludedSamples",
                f"{centered_filtered_out}.filtered-samples.txt",
                "--excludeTargets",
                extreme_gc_out,
                "--minTargetSize",
                "10",
                "--maxTargetSize",
                "10000",
                "--minMeanTargetRD",
                "10",
                "--maxMeanTargetRD",
                "1000",
                "--minMeanSampleRD",
                "25",
                "--maxMeanSampleRD",
                "1000",
                "--maxSdSampleRD",
                "300",
                "-o",
                centered_filtered_out,
            ]
        )

        PCA_out = f"{self.docker_output_base}/DATA.RD_PCA"
        self.run_xhmm_command(["--PCA", "-r", centered_filtered_out, "--PCAfiles", PCA_out])

        normalised_out = f"{self.docker_output_base}/DATA.PCA_normalized.txt"

        self.run_xhmm_command(
            [
                "--normalize",
                "-r",
                centered_filtered_out,
                "--PCAfiles",
                PCA_out,
                "--normalizeOutput",
                normalised_out,
                "--PCnormalizeMethod",
                "PVE_mean",
                "--PVE_mean_factor",
                "0.7",
            ]
        )

        sample_zscores = f"{self.docker_output_base}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt"

        self.run_xhmm_command(
            [
                "--matrix",
                "-r",
                normalised_out,
                "--centerData",
                "--centerType",
                "sample",
                "--zScoreData",
                "--outputExcludedTargets",
                f"{sample_zscores}.filtered-targets.txt",
                "--outputExcludedSamples",
                f"{sample_zscores}.filtered-samples.txt",
                "--maxSdTargetRD",
                "30",
                "-o",
                sample_zscores,
            ]
        )

        matched_filter_read_depth = f"{self.docker_output_base}/DATA.same_filtered.RD.txt"

        self.run_xhmm_command(
            [
                "--matrix",
                "-r",
                xhmm_read_depth_out,
                "--excludeTargets",
                f"{centered_filtered_out}.filtered-targets.txt",
                "--excludeTargets",
                f"{sample_zscores}.filtered-targets.txt",
                "--excludeSamples",
                f"{centered_filtered_out}.filtered-samples.txt",
                "--excludeSamples",
                f"{sample_zscores}.filtered-samples.txt",
                "-o",
                matched_filter_read_depth,
            ]
        )

        discover_cnv_out = f"{self.docker_output_base}/DATA.xcnv"
        self.run_xhmm_command(
            [
                "--discover",
                "-p",
                "/mnt/cnv-caller-resources/xhmm/params.txt",
                "-r",
                sample_zscores,
                "-R",
                matched_filter_read_depth,
                "-a",
                discover_cnv_out.replace(".xcnv", ".aux_xcnv"),
                "-s",
                f"{self.docker_output_base}/DATA",
                "-c",
                discover_cnv_out,
            ]
        )

        vcf_out = f"{self.docker_output_base}/DATA.vcf"

        self.run_xhmm_command(
            [
                "--genotype",
                "-p",
                "/mnt/cnv-caller-resources/xhmm/params.txt",
                "-r",
                sample_zscores,
                "-R",
                matched_filter_read_depth,
                "-g",
                discover_cnv_out,
                "-F",
                self.settings["ref_fasta"],
                "-v",
                vcf_out,
            ]
        )
