"""
Defaults used from XHMM documentation https://atgu.mgh.harvard.edu/xhmm/tutorial.shtml

- masks regions with extreme GC content 
- no parallism at the moment
- no repeat masking using PLINK/Seq
- Altered filtering to allow for high read depth
   - Max target read depth: 100,000
   - Max sample mean read depth: 100,000
   - Max sample SD read depth: 300
"""

import subprocess
import time
import os

from . import utils, base_classes


class XHMM(base_classes.BaseCNVTool):
    def __init__(self, capture, gene, start_time, normal_panel=True):
        self.run_type = "xhmm"
        super().__init__(capture, gene, start_time, normal_panel=normal_panel)
        self.extra_db_fields = ["id", "ref", "qual", "filter", "format_data", "info_data"]
        self.settings = {**self.settings, "docker_image": "stefpiatek/xhmm:1.0"}

    def run_gatk_command(self, args):
        """Create dir for output and runs a GATK command in the XHMM docker"""
        try:
            os.makedirs(f"{self.output_base}/{args[0]}")
        except FileExistsError:
            base_classes.logger.info(f"Folder {args[0]} already exists")

        base_classes.logger.info(f"Running  GATK step of XHMM: {args[0]} \n output: {args[-1]}")
        self.run_docker_subprocess(
            [
                "java",
                f"-Xmx{self.max_mem}g",
                f"-XX:ConcGCThreads={self.max_cpu}",
                "-jar",
                "GenomeAnalysisTK.jar",
                "-T",
                *args,
            ]
        )
        base_classes.logger.info(f"Completed  GATK step of XHMM: {args[0]} {args[-1]}")

    def parse_output_file(self, file_path, sample_id):
        if "PYTEST_CURRENT_TEST" in os.environ.keys():
            docker_file_path = file_path.replace(f"{base_classes.cnv_pat_dir}/tests", "/mnt")
        else:
            docker_file_path = file_path.replace(base_classes.cnv_pat_dir, "/mnt")

        if not os.path.exists(f"{file_path}.gz.tbi"):
            self.run_docker_subprocess(["bgzip", docker_file_path], docker_image="lethalfang/tabix:1.7")
            self.run_docker_subprocess(
                ["tabix", "-p", "vcf", f"{docker_file_path}.gz"], docker_image="lethalfang/tabix:1.7"
            )

        sample_vcf = self.run_docker_subprocess(
            ["bcftools", "view", "-c1", "-Ov", "--force-samples", "-s", sample_id, f"{docker_file_path}.gz"],
            docker_image="halllab/bcftools:v1.9",
            stdout=subprocess.PIPE,
        )
        # avoid Error response from daemon: containerd: container did not start before the specified timeout.
        time.sleep(5)

        vcf_data = str(sample_vcf.stdout, "utf-8").split("\n")
        cnvs = self.parse_vcf(vcf_data, sample_id)

        return cnvs

    def run_xhmm_command(self, args):
        """Runs xhmm command in docker"""
        base_classes.logger.info(f"Running  XHMM: {args[0]} \n output: {args[-1]}")
        self.run_docker_subprocess(["xhmm", *args])
        base_classes.logger.info(f"Completed  XHMM: {args[0]} {args[-1]}")

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
                "100000",
                "--minMeanSampleRD",
                "25",
                "--maxMeanSampleRD",
                "100000",
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
        sample_names = [f"{self.bam_to_sample[unknown_bam]}" for unknown_bam in self.settings["unknown_bams"]]
        output_paths = [f"{self.output_base}/DATA.vcf" for sample_name in sample_names]

        return output_paths, sample_names
