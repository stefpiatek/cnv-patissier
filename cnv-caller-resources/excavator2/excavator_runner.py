"""
Excavator2 alters its directory with the target files, so need to be run in one session

"""
import argparse
import os
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument("--output-base")
parser.add_argument("--max-mem")
parser.add_argument("--max-cpu")
args = parser.parse_args()

subprocess.run(
    [
        "perl",
        "TargetPerla.pl",
        f"{args.output_base}/SourceTarget.txt",
        f"{args.output_base}/capture.bed",
        "target",
        "50000",
        "hg19",
    ],
    check=True,
)


subprocess.run(
    [
        "perl",
        "EXCAVATORDataPrepare.pl",
        f"{args.output_base}/ExperimentalFilePrepare.txt",
        "--processors",
        args.max_cpu,
        "--target",
        "target",
        "--assembly",
        "hg19",
    ],
    check=True,
)

subprocess.run(
    [
        "perl",
        "EXCAVATORDataAnalysis.pl",
        f"{args.output_base}/ExperimentalFileData.txt",
        "--processors",
        args.max_cpu,
        "--assembly",
        "hg19",
        "--target",
        "target",
        "--mode",
        "pooling",
        "--output",
        f"{args.output_base}/results",
    ],
    check=True,
)
