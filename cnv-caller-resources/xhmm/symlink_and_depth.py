"""Painful hack to symlink bam files and then run depth of coverage"""
import argparse
import os
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument("--bams", nargs="*")
parser.add_argument("--symlink_bams", nargs="*")
parser.add_argument("--genome")
parser.add_argument("--capture_path")
parser.add_argument("--mem")
parser.add_argument("--cpu")
parser.add_argument("--output")
args = parser.parse_args()

# couldn't get list to parse nicely with nargs, 
original_bams = args.bams[0].split(" ")
symlink_bams = args.symlink_bams[0].split(" ")



for original_bam, symlink_bam in zip(original_bams, symlink_bams):
    os.symlink(original_bam, symlink_bam)

gatk_command = subprocess.Popen(
    [
        "java",
        "-" + args.mem,
        "-" + args.cpu,
        "-jar",
        "GenomeAnalysisTK-3.8-0.jar",
        "-T",
        "DepthOfCoverage",
        "-I",
        " ".join(symlink_bams),
        "-L",
        args.capture_path,
        "-R",
        args.genome,
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
        args.output,
    ],
)

gatk_command.wait()
