from scripts.base_classes import BaseCNVTool


class TestGetBAMHeader:
    pass


class TestParseVCF:
    def setup(self):
        self.del_vcf = f"tests/test_files/output_parsing/vcf/15384-del_segments.vcf"
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
        with open(self.del_vcf) as handle:
            vcf_dict = BaseCNVTool.parse_vcf_4_2(handle, "15384-del")
        assert vcf_dict[0] == self.del_expected_vcf_dict

    def test_dup(self):
        expected = dict(self.del_expected_vcf_dict)
        expected["alt"] = "DUP"
        expected["format_data"] = {"GT": "2", "CN": "3", "NP": "4", "QA": "118", "QS": "473", "QSE": "67", "QSS": "118"}
        expected["sample_id"] = "15384-dup"
        with open(self.del_vcf.replace("15384-del", "15384-dup")) as handle:
            vcf_dict = BaseCNVTool.parse_vcf_4_2(handle, "15384-dup")
        assert vcf_dict[0] == expected

    def test_no_cnv(self):
        with open(self.del_vcf.replace("15384-del", "15384-none")) as handle:
            vcf_dict = BaseCNVTool.parse_vcf_4_2(handle, None)
        assert vcf_dict == []


class TestUploadCNVCaller:
    pass


class TestUploadGene:
    pass


class TestUploadSamples:
    pass


class TestUploadPositiveCNVs:
    pass


class TestUploadCalledCNVs:
    pass
