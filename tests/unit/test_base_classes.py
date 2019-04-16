import pathlib
import subprocess

import pytest

from scripts.base_classes import BaseCNVTool
from scripts import models, utils

cnv_pat_dir = utils.get_cnv_patissier_dir()


def instance_data(instance):
    return {k: v for k, v in instance.__dict__.items() if k != "_sa_instance_state"}


@pytest.mark.usefixtures("db", "db_session")
class TestFilterCNVs:
    def setup(self):
        self.caller = BaseCNVTool("capture", "gene", "time")
        self.caller.run_type = "example_type"
        self.caller.extra_db_fields = ["extra"]
        self.gene_bed = [
            "chr1\t1000\t1500\tgene",
            "chr1\t2000\t2500\tgene",
            "chr1\t5000\t6000\tgene",
            "chr17\t1000\t1500\tgene_2",
        ]

    def test_filter_out(self):
        chrom_mismatch = dict(chrom="chr21", start="1000", end="1500", extra="dummy")
        outside_start_end = dict(chrom="chr17", start="2000", end="2500", extra="dummy")
        filtered = self.caller.filter_cnvs([chrom_mismatch, outside_start_end], self.gene_bed)
        assert filtered == []

    def test_filter_within(self):

        within = dict(chrom="chr17", start="1200", end="1300", extra="dummy")
        within_filtered = self.caller.filter_cnvs([within], self.gene_bed)

        start_within = dict(chrom="chr17", start="1200", end="2000", extra="dummy")
        start_within_filtered = self.caller.filter_cnvs([start_within], self.gene_bed)

        end_within = dict(chrom="chr17", start="500", end="1200", extra="dummy")
        end_within_filtered = self.caller.filter_cnvs([end_within], self.gene_bed)

        assert within_filtered[0]["start"] == "1200"
        assert start_within_filtered[0]["start"] == "1200"
        assert end_within_filtered[0]["end"] == "1200"

    def test_filter_span(self):
        span = dict(chrom="chr1", start="500", end="1600", extra="dummy")
        span_filtered = self.caller.filter_cnvs([span], self.gene_bed)

        span_multiple = dict(chrom="chr1", start="500", end="4000", extra="dummy")
        span_multiple_filtered = self.caller.filter_cnvs([span_multiple], self.gene_bed)

        assert span_filtered[0]["start"] == "500"
        assert len(span_multiple_filtered) == 1
        assert span_multiple_filtered[0]["start"] == "500"

    def test_multiple_cnvs(self):
        span_multiple = dict(chrom="chr1", start="500", end="4000", extra="dummy")
        within = dict(chrom="chr17", start="1200", end="1300", extra="dummy")

        filtered = self.caller.filter_cnvs([span_multiple, within], self.gene_bed)
        assert len(filtered) == 2
        assert filtered[0]["chrom"] == "chr1"
        assert filtered[1]["chrom"] == "chr17"

    def test_extra_fields(self):
        dummy = {
            "key1": "val1",
            "dict_key": {"sub_dict_key1": "sub_dict_val1", "sub_dict_key2": "sub_dict_val2"},
            "key2": "val2",
            "key3": "val3",
        }
        within = dict(chrom="chr17", start="1200", end="1300", extra=dummy)
        filtered = self.caller.filter_cnvs([within], self.gene_bed)
        print(filtered[0])
        assert filtered[0]["json_data"]["extra"] == dummy


@pytest.mark.usefixtures("db", "db_session")
class TestGetBAMHeader:
    def setup(self):
        self.caller = BaseCNVTool("capture", "gene", "time")
        self.caller.run_type = "example_type"
        self.caller.sample_to_bam = {
            "12S13548": "/mnt/test_files/input/bam_header.bam",
            "sample_mismatch": "/mnt/test_files/input/bam_header.bam",
        }
        self.caller.bam_mount = "/mnt/data/"

    def test_basic(self):
        with open("tests/test_files/input/bam_header.sam") as handle:
            header_list = handle.readlines()
        expected_header = "".join(header_list).replace("\n", "\r\n")
        output_header = self.caller.get_bam_header("12S13548")
        assert expected_header == output_header

    def test_sample_mismatch(self):
        with pytest.raises(AssertionError):
            self.caller.get_bam_header("sample_mismatch")


@pytest.mark.usefixtures("db", "db_session")
class TestGetMd5sum:
    def setup(self):
        self.caller = BaseCNVTool("capture", "gene", "time")
        self.caller.run_type = "example_type"

    def test_simple(self):
        expected = ("d41d8cd98f00b204e9800998ecf8427e", "cnv-pat/input/.gitkeep")
        output = self.caller.get_md5sum(pathlib.Path("input/.gitkeep").resolve())
        assert expected == output

    def test_missing_file(self):
        with pytest.raises(subprocess.CalledProcessError):
            output = self.caller.get_md5sum("does_not_exist.txt")


@pytest.mark.usefixtures("db", "db_session", "populate_db")
class TestPreRunSteps:
    def setup(self):
        self.test_file_prefix = f"{cnv_pat_dir}/tests/test_files/input/checks"

        with open(f"{self.test_file_prefix}/sample_sheet_working.txt", "w") as handle:
            handle.write("sample_id\tsample_path\n")
            handle.write(f"12S13548\t{cnv_pat_dir}/tests/test_files/input/bam_header.bam\n")

        self.caller = BaseCNVTool("ICR", "gene_1", "time")
        self.caller.bam_mount = "/mnt/data/"
        self.caller.run_type = "example_type"
        self.caller.sample_to_bam = {"12S13548": "/mnt/test_files/input/bam_header.bam"}
        self.header_file = "tests/test_files/input/bam_header.sam"
        with open("tests/test_files/input/bam_header.sam") as handle:
            header_list = handle.readlines()

        self.expected_header = {"12S13548": "".join(header_list).replace("\n", "\r\n")}

    def test_header(self):
        header = self.caller.prerun_steps(f"{self.test_file_prefix}/sample_sheet_working.txt")
        assert header == self.expected_header


@pytest.mark.usefixtures("db", "db_session", "populate_db")
class TestUploadAllMd5sums:
    def setup(self):
        self.caller = BaseCNVTool("capture", "gene", "time")
        self.caller.run_type = "example_type"

    def test_working(self):
        before_upload = self.caller.session.query(models.File).all()
        self.caller.upload_all_md5sums(1)
        after_upload = self.caller.session.query(models.File).all()
        assert len(before_upload) < len(after_upload)


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
        self.caller.bam_headers = {"12S13548": "header1", "10S21354": "header2"}
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
        """Existing instance should stay the same, 12S13548 shouldn't exist, then be uploaded
        Data returned should be a list of dictionaries for upload of known cnv information from sample sheet"""

        no_instance = self.caller.session.query(models.Sample).filter_by(name="12S13548").first()

        output_table = self.caller.upload_samples(self.sample_sheet)
        uploaded_1 = self.caller.session.query(models.Sample).filter_by(name="12S13548").first()
        assert not no_instance
        assert uploaded_1.name == "12S13548"
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
                    "name": "10S21354",
                    "gene_id": 1,
                    "path": "/mnt/data/181225_NB503215_run/analysis/Alignments/10S21354_sorted.bam",
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
            "sample_id": "10S21354",
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
