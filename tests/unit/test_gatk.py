from scripts.gatk import GATKBase


class TestParseOutputFile:
    def setup(self):
        self.del_vcf = f"tests/test_files/output_parsing/vcf/15384-del_segments.vcf"
        self.caller = GATKBase("capture", "gene_1", "time", normal_panel=True)
        self.del_expected_vcf_dict = {
            "chrom": "chr7",
            "start": "513245",
            "end": "514245",
            "id": "CNV_chr7_513245_514245",
            "ref": "N",
            "alt": "DEL",
            "qual": ".",
            "filter": ".",
            "info_data": {"END": "514245"},
            "format_data": {"GT": "1", "CN": "1", "NP": "4", "QA": "118", "QS": "473", "QSE": "67", "QSS": "118"},
            "sample_id": "15384-del",
        }

    def test_del(self):
        parsed = self.caller.parse_output_file(self.del_vcf, "15384-del")
        assert parsed[0] == self.del_expected_vcf_dict

    def test_dup(self):
        expected = dict(self.del_expected_vcf_dict)
        expected["alt"] = "DUP"
        expected["format_data"]["GT"] = "2"
        expected["format_data"]["CN"] = "3"
        expected["sample_id"] = "15384-dup"
        parsed = self.caller.parse_output_file(self.del_vcf.replace("15384-del", "15384-dup"), "15384-dup")

        assert parsed[0] == expected

    def test_no_cnv(self):
        parsed = self.caller.parse_output_file(self.del_vcf.replace("15384-del", "15384-none"), "15384-none")
        assert parsed == []
