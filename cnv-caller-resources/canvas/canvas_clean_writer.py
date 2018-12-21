"""
Writes a json for manually running canvas clean step
"""

import argparse
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument("--output-base")
parser.add_argument("--sample")
args = parser.parse_args()

json_input = [
    "{",
    '  "CleanedPath": "${Output}' f"/TempCNV_{args.sample}/{args.sample}." 'cleaned",',
    '  "FfpePath": null',
    "}",
]

with open(f"{args.output_base}/Checkpoints/03-CanvasClean.json", "w") as handle:
    for line in json_input:
        handle.write(f"{line}\n")
