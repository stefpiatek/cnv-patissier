import pathlib

from scripts.exome_depth import ExomeDepthBase


class TestParseOutputFile:
    def setup(self):
        self.caller = ExomeDepthBase("capture", "gene_1", "time", normal_panel=True)
        self.caller.settings["chromosome_prefix"] = "chr"
        self.output_base = pathlib.Path("tests/test_files/output_parsing/exome_depth/")
        self.del_expected_output = [
            {
                "chrom": "chr2",
                "start": "5700",
                "end": "5750",
                "alt": "DEL",
                "sample_id": "del",
                "start.p": "34",
                "end.p": "42",
                "nexons": "9",
                "id": "chr7:5700-5750",
                "BF": "28.3",
                "reads.expected": "7791",
                "reads.observed": "5101",
                "reads.ratio": "0.655",
            },
            {
                "chrom": "chr17",
                "start": "100",
                "end": "4000",
                "alt": "DEL",
                "sample_id": "del",
                "start.p": "62",
                "end.p": "62",
                "nexons": "1",
                "id": "chr2:100-4000",
                "BF": "4.56",
                "reads.expected": "378",
                "reads.observed": "220",
                "reads.ratio": "0.582",
            },
        ]

    def test_del(self):
        parsed = self.caller.parse_output_file(self.output_base / "del.txt", "del")
        assert parsed == self.del_expected_output

    def test_dup(self):
        expected_output = list(self.del_expected_output)
        for row in expected_output:
            row["alt"] = "DUP"
            row["sample_id"] = "dup"

        parsed = self.caller.parse_output_file(self.output_base / "dup.txt", "dup")
        assert parsed == expected_output

    def test_normal(self):
        parsed = self.caller.parse_output_file(self.output_base / "normal.txt", "normal")
        assert parsed == []
