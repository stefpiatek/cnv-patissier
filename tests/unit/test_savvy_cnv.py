import pathlib

import pytest

from scripts.savvy_cnv import SavvyCNV


@pytest.mark.usefixtures("db", "db_session")
class TestParseOutputFile:
    def setup(self):
        self.caller = SavvyCNV("capture", "gene_1", "time", normal_panel=True)
        self.caller.docker_output_base = "/mnt/output/capture/time/savvycnv/gene_1"
        self.output_file = pathlib.Path("tests/test_files/output_parsing/savvy_cnv/cnv_calls.txt")
        self.del_expected_output = [
            {
                "chrom": "chr2",
                "start": "5700",
                "end": "5750",
                "alt": "DEL",
                "supporting_chunks": "10",
                "width_of_cnv": "618",
                "phred_score": "708.2041717",
                "normalised_phred": "1.145961443",
                "relative_dosage": "0.555971128",
                "coverage_filename": "/mnt/output/capture/time/savvycnv/gene_1/CoverageBinner/del.coverageBinner",
                "sample_id": "del",
            },
            {
                "chrom": "chr17",
                "start": "100",
                "end": "4000",
                "alt": "DEL",
                "supporting_chunks": "2",
                "width_of_cnv": "5",
                "phred_score": "115.2891442",
                "normalised_phred": "23.05782883",
                "relative_dosage": "0.555971128",
                "coverage_filename": "/mnt/output/capture/time/savvycnv/gene_1/CoverageBinner/del.coverageBinner",
                "sample_id": "del",
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
            row["relative_dosage"] = "1.339272793"
            row["coverage_filename"] = "/mnt/output/capture/time/savvycnv/gene_1/CoverageBinner/dup.coverageBinner"

        parsed = self.caller.parse_output_file(self.output_file, "dup")
        assert parsed == expected_output

    def test_normal(self):
        parsed = self.caller.parse_output_file(self.output_file, "normal")
        assert parsed == []
