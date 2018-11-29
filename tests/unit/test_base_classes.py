from scripts.utils import get_cnv_patissier_dir
from scripts.base_classes import BaseCNVTool

cnv_pat_dir = get_cnv_patissier_dir()


class TestBaseCNVToolParseVCF:
    def setup(self):
        self.del_vcf = (
            f"{cnv_pat_dir}/tests/test_files/output/cohort/date-run/caller/gene/extra-layers/17340-del_segments.vcf"
        )
        self.del_expected_vcf_dict = {
            "chrom": "chr17",
            "start": "41251742",
            "end": "41258601",
            "id": "CNV_chr17_41251742_41258601",
            "ref": "N",
            "alt": "DEL",
            "qual": ".",
            "filter": ".",
            "info": "END=41258601",
            "format": "GT:CN:NP:QA:QS:QSE:QSS",
            "data": "1:1:4:118:473:67:118",
            "cnv_caller": "caller",
            "gene": "gene",
            "sample_id": "17340-del",
            "cohort": "cohort",
        }

    def test_del(self):
        vcf_dict = BaseCNVTool.parse_vcf_4_2(self.del_vcf)
        assert vcf_dict[0] == self.del_expected_vcf_dict

    def test_dup(self):
        expected = dict(self.del_expected_vcf_dict)
        expected["alt"] = "DUP"
        expected["data"] = "2:3:4:118:473:67:118"
        expected["sample_id"] = "17340-dup"
        vcf_dict = BaseCNVTool.parse_vcf_4_2(self.del_vcf.replace("17340-del", "17340-dup"))
        assert vcf_dict[0] == expected

    def test_no_cnv(self):
        vcf_dict = BaseCNVTool.parse_vcf_4_2(self.del_vcf.replace("17340-del", "17340-none"))
        assert vcf_dict == []
