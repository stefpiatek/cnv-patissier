import pathlib

import pytest

from scripts.cnv_kit import CNVKit
from scripts import models


class TestParseOutputFile:
    def setup(self):
        self.caller = CNVKit("capture", "gene_1", "time")
        self.output_base = pathlib.Path("tests/test_files/output_parsing/cnvkit")
        self.del_expected_output = [
            dict(
                chrom="chr2",
                start="5700",
                end="5750",
                log2="-0.683746",
                depth="214.013",
                probes="1363",
                weight="113.171",
                cn=1,
                alt="DEL",
                sample_id="test",
            ),
            dict(
                chrom="chr17",
                start="100",
                end="4000",
                log2="-0.272278",
                depth="498.159",
                probes="436",
                weight="75.3287",
                cn=1,
                alt="DEL",
                sample_id="test",
            ),
        ]

    def test_del(self):
        parsed = self.caller.parse_output_file(self.output_base / "sample_1_del_sorted_calls.ci", "test")
        assert parsed == self.del_expected_output

    def test_dup(self):
        expected_output = list(self.del_expected_output)
        for row in expected_output:
            row["cn"] = 3
            row["alt"] = "DUP"
        parsed = self.caller.parse_output_file(self.output_base / "sample_1_dup_sorted_calls.ci", "test")
        assert parsed == expected_output

    def test_normal(self):
        parsed = self.caller.parse_output_file(self.output_base / "sample_1_normal_sorted_calls.ci", "test")
        assert parsed == []
