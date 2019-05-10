import pathlib

from scripts.panelcn_mops import panelcnMOPS


class TestParseOutputFile:
    def setup(self):
        self.caller = panelcnMOPS("capture", "gene_1", "time")
        self.caller.settings = {"chromosome_prefix": "chr"}
        self.output = pathlib.Path("tests/test_files/output_parsing/panelcn_mops/results.txt")
        self.del_expected_output = [
            {
                "chrom": "chr2",
                "start": "5700",
                "end": "5750",
                "gene": "test",
                "exon": "test.NA.chr2.5700.5750",
                "alt": "DEL",
                "sample_id": "sample_1",
                "rc": "1430",
                "medrc": "2546",
                "rc.norm": "1537",
                "medrc.norm": "2508",
                "lowqual": "",
                "cn": "CN1",
            },
            {
                "chrom": "chr17",
                "start": "100",
                "end": "4000",
                "gene": "test",
                "exon": "test.NA.chr17.100.4000",
                "alt": "DEL",
                "sample_id": "sample_1",
                "rc": "1430",
                "medrc": "2546",
                "rc.norm": "1537",
                "medrc.norm": "2508",
                "lowqual": "",
                "cn": "CN1",
            },
        ]

    def test_del(self):
        parsed = self.caller.parse_output_file(self.output, "sample_1")
        assert parsed == self.del_expected_output

    def test_dup(self):
        expected_output = list(self.del_expected_output)
        for row in expected_output:
            row["rc"] = "5430"
            row["rc.norm"] = "3537"
            row["alt"] = "DUP"
            row["sample_id"] = "sample_2"
            row["cn"] = "CN3"
        parsed = self.caller.parse_output_file(self.output, "sample_2")
        assert parsed == expected_output

    def test_normal(self):
        parsed = self.caller.parse_output_file(self.output, "sample_3")
        assert parsed == []
