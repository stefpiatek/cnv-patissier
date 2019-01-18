import pytest

from scripts.base_classes import BaseCNVTool
from scripts import models


def instance_data(instance):
    return {k: v for k, v in instance.__dict__.items() if k != "_sa_instance_state"}


@pytest.mark.usefixtures("db", "db_session")
class TestFilterCNVs:
    pass


@pytest.mark.usefixtures("db", "db_session")
class TestGetBAMHeader:
    def setup(self):
        self.caller = BaseCNVTool("capture", "gene", "time")
        self.caller.run_type = "example_type"
        self.caller.sample_to_bam = {"sample_1": "/mnt/test_files/input/bam_header.bam"}
        self.caller.bam_mount = "/mnt/data/"

    def test_basic(self):
        with open("tests/test_files/input/bam_header.sam") as handle:
            header_list = handle.readlines()
        expected_header = "".join(header_list).replace("\n", "\r\n")
        output_header = self.caller.get_bam_header("sample_1")
        assert expected_header == output_header


@pytest.mark.usefixtures("db", "db_session", "populate_db")
class TestUploadCNVCaller:
    def setup(self):
        self.caller = BaseCNVTool("capture", "gene", "time")
        self.caller.run_type = "example_type"

    def test_new_caller(self):
        no_caller = self.caller.session.query(models.Caller).filter_by(name="example_type").first()
        self.caller.upload_cnv_caller()
        self.caller.session.commit()
        existing_caller = self.caller.session.query(models.Caller).filter_by(name="example_type").first()
        assert no_caller is None
        assert existing_caller

    def test_reupload_caller(self):
        self.caller.upload_cnv_caller()
        self.caller.session.commit()
        existing_caller = self.caller.session.query(models.Caller).filter_by(name="example_type").all()
        assert len(existing_caller) == 1


@pytest.mark.usefixtures("db", "db_session", "populate_db")
class TestUploadGene:
    def setup(self):
        self.caller = BaseCNVTool("ICR", "gene_1", "time")
        self.caller.run_type = "example_type"
        self.caller.settings = {"capture_path": "/mnt/tests/test_files/input/capture.bed", "genome_build_name": "hg19"}

    def test_new_gene(self):
        self.caller.gene = "gene_3"
        no_instance = self.caller.session.query(models.Gene).filter_by(name="gene_3").first()
        self.caller.upload_gene()
        gene_instance = self.caller.session.query(models.Gene).filter_by(name="gene_3").first()
        assert no_instance is None
        assert gene_instance.chrom == "chr3"
        assert gene_instance.start == 1500
        assert gene_instance.end == 1950

    def test_existing_gene(self):
        self.caller.gene = "gene_1"
        previous_genes = self.caller.session.query(models.Gene).filter_by(name="gene_1").all()
        prevous_data = [instance_data(previous_gene) for previous_gene in previous_genes]
        self.caller.upload_gene()
        current_genes = self.caller.session.query(models.Gene).filter_by(name="gene_1").all()
        current_data = [instance_data(current_gene) for current_gene in current_genes]
        assert prevous_data == current_data

    def test_updating_gene(self):
        self.caller.gene = "gene_2"
        previous_gene = self.caller.session.query(models.Gene).filter_by(name="gene_2").first()
        previous_data = instance_data(previous_gene)
        self.caller.upload_gene()
        current_gene = self.caller.session.query(models.Gene).filter_by(name="gene_2").first()
        assert previous_data["end"] == 9000
        assert current_gene.end == 5900


@pytest.mark.usefixtures("db", "db_session", "populate_db")
class TestUploadSamples:
    def setup(self):
        self.caller = BaseCNVTool("ICR", "gene_1", "time")
        self.caller.run_type = "example_type"
        self.caller.sample_to_bam = {
            "12S13548": "/mnt/test_files/input/bam_header.bam",
            "10S21354": "/mnt/test_files/input/bam_header.bam",
        }
        self.caller.bam_mount = "/mnt/data/"
        self.caller.settings = {"genome_build_name": "hg19"}
        self.header_file = "tests/test_files/input/bam_header.sam"
        with open("tests/test_files/input/bam_header.sam") as handle:
            header_list = handle.readlines()
        self.expected_header = "".join(header_list).replace("\n", "\r\n")
        self.sample_sheet = "tests/test_files/input/gene_1.txt"
        self.expected_output = [
            {
                "cnv": {"alt": "DUP", "genome_build": "hg19", "chrom": "chr17", "start": "1200", "end": "1500"},
                "sample_defaults": {
                    "name": "12S13548",
                    "gene_id": 1,
                    "path": "/mnt/data/181225_NB503215_run/analysis/Alignments/12S13548_sorted.bam",
                },
            }
        ]

    def test_basic(self):
        """Existing instance should stay the same, 10S shouldn't exist, then be uploaded
        Data returned should be a list of dictionaries for upload of known cnv information from sample sheet"""

        existing_1 = self.caller.session.query(models.Sample).filter_by(name="12S13548").first()
        existing_data_1 = instance_data(existing_1)
        no_instance = self.caller.session.query(models.Sample).filter_by(name="10S21354").first()

        output_table = self.caller.upload_samples(self.sample_sheet)
        uploaded_1 = self.caller.session.query(models.Sample).filter_by(name="12S13548").first()
        uploaded_2 = self.caller.session.query(models.Sample).filter_by(name="10S21354").first()
        assert existing_data_1 == instance_data(uploaded_1)
        assert no_instance is None
        assert uploaded_2.name == "10S21354"
        assert uploaded_2.bam_header == existing_data_1["bam_header"]
        assert uploaded_2.result_type == "normal-panel"
        assert output_table == self.expected_output


@pytest.mark.usefixtures("db", "db_session", "populate_db")
class TestUploadPositiveCNVs:
    def setup(self):
        self.caller = BaseCNVTool("ICR", "gene_1", "time")
        self.caller.run_type = "example_type"
        self.known_cnv = [
            {
                "cnv": {"alt": "DUP", "genome_build": "hg19", "chrom": "chr17", "start": "1200", "end": "1500"},
                "sample_defaults": {
                    "name": "12S13548",
                    "gene_id": 1,
                    "path": "/mnt/data/181225_NB503215_run/analysis/Alignments/12S13548_sorted.bam",
                },
            }
        ]

    def test_existing(self):
        existing_cnv = self.caller.session.query(models.CNV).filter_by(start=1200, end=1500, chrom="chr17").first()
        existing_cnv_data = instance_data(existing_cnv)
        existing_known_cnv = self.caller.session.query(models.KnownCNV).filter_by(cnv_id=1, sample_id=1).first()
        existing_known_cnv_data = instance_data(existing_known_cnv)

        self.caller.upload_positive_cnvs(self.known_cnv)
        uploaded_cnv = self.caller.session.query(models.CNV).filter_by(start=1200, end=1500, chrom="chr17").first()
        uploaded_known_cnv = self.caller.session.query(models.KnownCNV).filter_by(cnv_id=1, sample_id=1).first()

        assert existing_cnv_data == instance_data(uploaded_cnv)
        assert existing_known_cnv_data == instance_data(uploaded_known_cnv)

    def test_create_new(self):
        self.known_cnv[0]["cnv"]["start"] = "1700"
        self.known_cnv[0]["cnv"]["end"] = "1900"
        self.caller.upload_positive_cnvs(self.known_cnv)
        uploaded_cnv = self.caller.session.query(models.CNV).filter_by(start=1700, end=1900, chrom="chr17").first()
        uploaded_known_cnv = (
            self.caller.session.query(models.KnownCNV).filter_by(cnv_id=uploaded_cnv.id, sample_id=1).first()
        )

        assert uploaded_cnv
        assert uploaded_cnv.id != 1
        assert uploaded_known_cnv.id > 2


@pytest.mark.usefixtures("db", "db_session", "populate_db")
class TestUploadCalledCNV:
    def setup(self):
        self.caller = BaseCNVTool("ICR", "gene_1", "time")
        self.caller.run_type = "first_caller"
        self.caller.settings = {"genome_build_name": "hg19"}
        self.cnv_call = {
            "start": "10",
            "end": "120",
            "chrom": "chr1",
            "sample_id": "12S13548",
            "alt": "DEL",
            "json_data": {"extra_field1": "extra_data1", "extra_field2": "extra_data2"},
        }

    def test_new_cnv(self):
        self.caller.upload_called_cnv(self.cnv_call)
        uploaded_cnv = self.caller.session.query(models.CNV).filter_by(start=10, end=120, chrom="chr1").first()
        uploaded_called_cnv = (
            self.caller.session.query(models.CalledCNV)
            .filter_by(caller_id=1, sample_id=1, cnv_id=uploaded_cnv.id)
            .first()
        )
        assert uploaded_cnv.id > 2
        assert uploaded_cnv.alt == "DEL"
        assert uploaded_called_cnv.id > 2
        assert uploaded_called_cnv.json_data == '{"extra_field1": "extra_data1", "extra_field2": "extra_data2"}'
