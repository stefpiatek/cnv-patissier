import pathlib

from scripts.excavator2 import Excavator2


class TestParseOutputFile:
    def setup(self):
        self.caller = Excavator2("capture", "gene_1", "time")
        self.output_base = pathlib.Path("tests/test_files/output_parsing/excavator2/")

        self.del_expected_output = [
            {
                "chrom": "chr7",
                "start": "513245",
                "end": "514245",
                "id": ".",
                "ref": "c",
                "alt": "DEL",
                "qual": ".",
                "filter": "PASS",
                "info_data": {"calling": "IMPRECISE", "SVTYPE": "CNV", "END": "514245", "SVLEN": "8850000"},
                "format_data": {"GT": "1/1", "CN": "1", "CNF": "2.94", "FCL": "1", "FCP": "1"},
                "sample_id": "sample_1",
            },
            {
                "chrom": "chr2",
                "start": "1250",
                "end": "1350",
                "id": ".",
                "ref": "c",
                "alt": "DEL",
                "qual": ".",
                "filter": "PASS",
                "info_data": {"calling": "IMPRECISE", "SVTYPE": "CNV", "END": "1350", "SVLEN": "8850000"},
                "format_data": {"GT": "1/1", "CN": "1", "CNF": "2.94", "FCL": "1", "FCP": "1"},
                "sample_id": "sample_1",
            },
        ]

    def test_del(self):
        parsed = self.caller.parse_output_file(self.output_base / "sample_1.vcf", "sample_1")
        assert parsed == self.del_expected_output

    def test_dup(self):
        expected_output = list(self.del_expected_output)
        for row in expected_output:
            row["format_data"]["CN"] = "3"
            row["alt"] = "DUP"
            row["sample_id"] = "sample_2"

        parsed = self.caller.parse_output_file(self.output_base / "sample_2.vcf", "sample_2")
        assert parsed == expected_output

    def test_normal(self):
        parsed = self.caller.parse_output_file(self.output_base / "sample_3.vcf", "sample_3")
        assert parsed == []
