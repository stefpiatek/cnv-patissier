import pathlib

import pytest

from scripts import models
from scripts.base_classes import cnv_pat_dir
from scripts.xhmm import XHMM


@pytest.mark.usefixtures("db", "db_session", "cleanup_after_xhmm")
class TestParseOutputFile:
    def setup(self):
        self.caller = XHMM("capture", "gene_1", "time")
        # because needs to load docker - use dummy direcetory for bam mount
        self.caller.bam_mount = f"{cnv_pat_dir}/logs"
        self.caller.settings = {"chromosome_prefix": "chr"}
        self.output = f"{cnv_pat_dir}/tests/test_files/output_parsing/xhmm/DATA.vcf"
        self.del_expected_output = [
            {
                "chrom": "chr2",
                "start": "5700",
                "end": "5750",
                "id": "chr2:5700-5750",
                "ref": "<DIP>",
                "alt": "DEL",
                "qual": ".",
                "filter": ".",
                "info_data": {
                    "calling": "IMPRECISE",
                    "AC": "1,0",
                    "AF": "0.01,0",
                    "AN": "1",
                    "END": "5750",
                    "SVLEN": "3496",
                    "SVTYPE": "CNV",
                    "TPOS": "5700",
                    "TEND": "5750",
                    "NUMT": "2",
                    "GQT": "78",
                    "PREVTARGSTART": "47643435",
                    "PREVTARGEND": "47643569",
                    "POSTTARGSTART": "48030559",
                    "POSTTARGEND": "48030825",
                },
                "format_data": {
                    "GT": "1",
                    "NDQ": "85",
                    "DQ": "0",
                    "EQ": "78,0",
                    "SQ": "85,0",
                    "NQ": "0,99",
                    "LQ": "24,0",
                    "RQ": "39,0",
                    "PL": "85,0,255",
                    "RD": "-7.86",
                    "ORD": "257.96",
                    "DSCVR": "Y",
                },
                "sample_id": "sample_1",
            },
            {
                "chrom": "chr17",
                "start": "100",
                "end": "4000",
                "id": "chr17:100-4000",
                "ref": "<DIP>",
                "alt": "DEL",
                "qual": ".",
                "filter": ".",
                "info_data": {
                    "calling": "IMPRECISE",
                    "AC": "1,0",
                    "AF": "0,0.01",
                    "AN": "1",
                    "END": "4000",
                    "SVLEN": "42028",
                    "SVTYPE": "CNV",
                    "TPOS": "100",
                    "TEND": "4000",
                    "NUMT": "4",
                    "GQT": "54",
                    "POSTTARGSTART": "52436304",
                    "POSTTARGEND": "52436438",
                },
                "format_data": {
                    "GT": "1",
                    "NDQ": "85",
                    "DQ": "0",
                    "EQ": "0,54",
                    "SQ": "0,99",
                    "NQ": "99,0",
                    "LQ": "0,54",
                    "RQ": "0,71",
                    "PL": "137,255,0",
                    "RD": "5.76",
                    "ORD": "586.59",
                    "DSCVR": "Y",
                },
                "sample_id": "sample_1",
            },
        ]

    def test_del(self):
        parsed = self.caller.parse_output_file(self.output, "sample_1")
        assert parsed == self.del_expected_output

    def test_dup(self):
        expected_output = list(self.del_expected_output)
        for row in expected_output:
            row["alt"] = "DUP"
            row["sample_id"] = "sample_2"
            row["format_data"]["GT"] = "2"
            row["format_data"]["NDQ"] = "99"
            row["info_data"]["AC"] = "0,1"
        parsed = self.caller.parse_output_file(self.output, "sample_2")
        assert parsed == expected_output

    def test_normal(self):
        parsed = self.caller.parse_output_file(self.output, "sample_3")

        assert parsed == []
