import pathlib

import pytest

from scripts.decon import DECoN


@pytest.mark.usefixtures("db", "db_session")
class TestParseOutputFile:
    def setup(self):
        self.caller = DECoN("capture", "gene_1", "time", normal_panel=True)
        self.caller.settings["chromosome_prefix"] = "chr"
        self.caller.bam_to_sample = {
            "/mnt/data/del_sorted.bam": "del",
            "/mnt/data/dup_sorted.bam": "dup",
            "/mnt/data/normal_sorted.bam": "normal",
        }
        self.output_file = pathlib.Path("tests/test_files/output_parsing/decon/cnv_calls_all.txt")
        self.del_expected_output = [
            {
                "CNV.ID": "1",
                "Correlation": "0.9972",
                "N.comp": "22",
                "chrom": "chr2",
                "start": "5700",
                "end": "5750",
                "alt": "DEL",
                "sample_id": "del",
                "Start.b": "34",
                "End.b": "42",
                "N.exons": "9",
                "Genomic.ID": "chrchr2:5700-5750",
                "BF": "28.3",
                "Reads.expected": "7791",
                "Reads.observed": "5101",
                "Reads.ratio": "0.655",
                "Gene": "example_gene",
            },
            {
                "CNV.ID": "2",
                "Correlation": "0.9971",
                "N.comp": "14",
                "chrom": "chr17",
                "start": "100",
                "end": "4000",
                "alt": "DEL",
                "sample_id": "del",
                "Start.b": "62",
                "End.b": "62",
                "N.exons": "1",
                "Genomic.ID": "chrchr17:100-4000",
                "BF": "4.56",
                "Reads.expected": "378",
                "Reads.observed": "220",
                "Reads.ratio": "0.582",
                "Gene": "example_gene",
            },
        ]

    def test_del(self):
        parsed = self.caller.parse_output_file(self.output_file, "del")
        assert parsed == self.del_expected_output

    def test_dup(self):
        expected_output = list(self.del_expected_output)
        for row in expected_output:
            row["alt"] = "DUP"
            row["sample_id"] = "dup"

        parsed = self.caller.parse_output_file(self.output_file, "dup")
        assert parsed == expected_output

    def test_normal(self):
        parsed = self.caller.parse_output_file(self.output_file, "normal")
        assert parsed == []
