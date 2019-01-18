import pathlib

import pytest

from scripts.copywriter import Copywriter
from scripts import models


class TestParseOutputFile:
    def setup(self):
        self.caller = Copywriter("capture", "gene_1", "time")
        self.caller.bam_to_sample = {f"/mnt/data/sample_{i}_sorted.bam": f"sample_{i}" for i in [1, 2, 3, "control"]}
        self.caller.settings = {"chromosome_prefix": "chr"}
        self.output = pathlib.Path("tests/test_files/output_parsing/copywriter/results.txt")
        self.del_expected_output = [
            {
                "chrom": "chr2",
                "start": "5700",
                "end": "5750",
                "num.mark": "4301",
                "seg.mean": -2.6262,
                "alt": "DEL",
                "sample_id": "sample_1",
                "control_id": "sample_control",
                "unknown": "sample_1_sorted.bam",
            },
            {
                "chrom": "chr17",
                "start": "100",
                "end": "4000",
                "num.mark": "4575",
                "seg.mean": -2.6551,
                "alt": "DEL",
                "sample_id": "sample_1",
                "control_id": "sample_control",
                "unknown": "sample_1_sorted.bam",
            },
        ]

    def test_del(self):
        parsed = self.caller.parse_output_file(self.output, "sample_1")
        assert parsed == self.del_expected_output

    def test_dup(self):
        expected_output = list(self.del_expected_output)
        for row in expected_output:
            row["seg.mean"] *= -1
            row["alt"] = "DUP"
            row["sample_id"] = "sample_2"
            row["unknown"] = "sample_2_sorted.bam"
        parsed = self.caller.parse_output_file(self.output, "sample_2")
        assert parsed == expected_output

    def test_normal(self):
        parsed = self.caller.parse_output_file(self.output, "sample_3")
        assert parsed == []
