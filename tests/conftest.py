import datetime
import json
import os
import subprocess

import pytest

from scripts import models
from scripts.db_session import DbSession
from scripts.base_classes import Queries, BaseCNVTool

# global application scope.  create Session class, engine


@pytest.yield_fixture(scope="class", autouse=True)
def cleanup_after_xhmm():
    vcf_path = "tests/test_files/output_parsing/xhmm/DATA.vcf"
    yield
    # teardown
    subprocess.run(["gunzip", f"{vcf_path}.gz"], check=True)
    subprocess.run(["rm", f"{vcf_path}.gz.tbi"], check=True)


@pytest.yield_fixture(scope="session")
def db():
    """Session-wide test database."""

    db_path = "tests/test.sqlite"
    DbSession.create_test_db(db_path)

    yield

    os.remove(db_path)


@pytest.yield_fixture(scope="class", autouse=True)
def db_session(db):
    """
    Creates a new database session for a test.
    https://docs.sqlalchemy.org/en/latest/orm/session_transaction.html#joining-a-session-into-an-external-transaction-such-as-for-test-suites
    """
    # create db_session
    DbSession.create_test_session()

    yield
    # not sure about rollback, maybe look into later but not really bothered for now


@pytest.fixture(scope="class", autouse=True)
def populate_db(db):
    caller = BaseCNVTool("capture", "gene", "time")
    session = caller.session

    # callers
    Queries.update_or_create(models.Caller, session, defaults={"id": 1}, name="first_caller")
    Queries.update_or_create(models.Caller, session, defaults={"id": 2}, name="second_caller")

    # cnvs
    cnv_1 = {"alt": "DUP", "chrom": "chr17", "start": "1200", "end": "1500", "genome_build": "hg19"}
    Queries.update_or_create(models.CNV, session, defaults={"id": 1}, **cnv_1)
    cnv_2 = {"alt": "DEL", "chrom": "chr2", "start": "5700", "end": "5900", "genome_build": "hg19"}
    Queries.update_or_create(models.CNV, session, defaults={"id": 2}, **cnv_2)

    # genes
    gene_1 = {
        "capture": "ICR",
        "chrom": "chr17",
        "start": "1000",
        "end": "3000",
        "genome_build": "hg19",
        "name": "gene_1",
    }
    Queries.update_or_create(models.Gene, session, defaults={"id": 1}, **gene_1)
    gene_2 = {
        "capture": "ICR",
        "chrom": "chr2",
        "start": "5000",
        "end": "9000",
        "genome_build": "hg19",
        "name": "gene_2",
    }
    Queries.update_or_create(models.Gene, session, defaults={"id": 2}, **gene_2)

    # samples
    bam_header = "\r\n".join(
        [
            "@HD\tVN:1.3\tSO:coordinate",
            "@SQ\tSN:chrM\tLN:16571",
            "@SQ\tSN:chr1\tLN:249250621",
            "@SQ\tSN:chr2\tLN:243199373",
            "@SQ\tSN:chr3\tLN:198022430",
            "@SQ\tSN:chr4\tLN:191154276",
            "@SQ\tSN:chr5\tLN:180915260",
            "@SQ\tSN:chr6\tLN:171115067",
            "@SQ\tSN:chr7\tLN:159138663",
            "@SQ\tSN:chr8\tLN:146364022",
            "@SQ\tSN:chr9\tLN:141213431",
            "@SQ\tSN:chr10\tLN:135534747",
            "@SQ\tSN:chr11\tLN:135006516",
            "@SQ\tSN:chr12\tLN:133851895",
            "@SQ\tSN:chr13\tLN:115169878",
            "@SQ\tSN:chr14\tLN:107349540",
            "@SQ\tSN:chr15\tLN:102531392",
            "@SQ\tSN:chr16\tLN:90354753",
            "@SQ\tSN:chr17\tLN:81195210",
            "@SQ\tSN:chr18\tLN:78077248",
            "@SQ\tSN:chr19\tLN:59128983",
            "@SQ\tSN:chr20\tLN:63025520",
            "@SQ\tSN:chr21\tLN:48129895",
            "@SQ\tSN:chr22\tLN:51304566",
            "@SQ\tSN:chrX\tLN:155270560",
            "@SQ\tSN:chrY\tLN:59373566",
            "@RG\tID:18\tCN:GOSH\tDS:2018-12-25\tDT:2019-01-01\tLB:L001\tPG:PipelineV1\tPL:NB503215\tSM:12S13548",
            "@PG\tID:18\tPN:bwa\tCL:bwa\tmem\t-M\t-t\t25\t-R\t@RG\tVN:0.7.15-r1140",
            "@PG\tID:SAMBLASTER\tCL:samblaster\t-i\tstdin\t-o\tstdout\tVN:0.1.24\r\n",
        ]
    )
    sample_1 = {
        "bam_header": bam_header,
        "gene_id": 1,
        "name": "12S13548",
        "path": "/mnt/data/181225_NB503215_run/analysis/Alignments/12S13548_sorted.bam",
        "result_type": "positive",
    }
    Queries.update_or_create(models.Sample, session, defaults={"id": 1}, **sample_1)
    sample_2 = {
        "bam_header": bam_header.replace("12S13548", "92S13548"),
        "gene_id": 1,
        "name": "92S13548",
        "path": "/mnt/data/181225_NB503215_run/analysis/Alignments/92S13548_sorted.bam",
        "result_type": "normal",
    }
    Queries.update_or_create(models.Sample, session, defaults={"id": 2}, **sample_2)
    sample_3 = {
        "bam_header": bam_header.replace("12S13548", "02S13548"),
        "gene_id": 1,
        "name": "02S13548",
        "path": "/mnt/data/181225_NB503215_run/analysis/Alignments/02S13548_sorted.bam",
        "result_type": "normal-panel",
    }
    Queries.update_or_create(models.Sample, session, defaults={"id": 3}, **sample_3)
    sample_4 = {
        "bam_header": bam_header.replace("12S13548", "2S13548"),
        "gene_id": 2,
        "name": "2S13548",
        "path": "/mnt/data/181225_NB503215_run/analysis/Alignments/2S13548_sorted.bam",
        "result_type": "positive",
    }
    Queries.update_or_create(models.Sample, session, defaults={"id": 3}, **sample_4)

    # runs
    run_1 = {
        "caller_id": 1,
        "gene_id": 1,
        "start_time": datetime.datetime(2019, 1, 1, 11, 34, 59),
        "end_time": datetime.datetime(2019, 1, 1, 13, 4, 1),
    }
    Queries.update_or_create(models.Run, session, defaults={"id": 1}, **run_1)
    run_2 = {
        "caller_id": 1,
        "gene_id": 2,
        "start_time": datetime.datetime(2019, 1, 1, 14, 34, 59),
        "end_time": datetime.datetime(2019, 1, 1, 16, 4, 1),
    }
    Queries.update_or_create(models.Run, session, defaults={"id": 1}, **run_1)

    # known_cnvs
    Queries.update_or_create(models.KnownCNV, session, defaults={"id": 1}, cnv_id=1, sample_id=1)
    Queries.update_or_create(models.KnownCNV, session, defaults={"id": 2}, cnv_id=2, sample_id=4)

    # called_cnvs
    called_cnv_1 = {"caller_id": 1, "cnv_id": 1, "sample_id": 1, "json_data": json.dumps({"depth": "12"})}
    Queries.update_or_create(models.CalledCNV, session, defaults={"id": 1}, **called_cnv_1)
    called_cnv_2 = {"caller_id": 1, "cnv_id": 2, "sample_id": 2, "json_data": json.dumps({"depth": "120"})}
    Queries.update_or_create(models.CalledCNV, session, defaults={"id": 2}, **called_cnv_2)

    # files
    file_1 = {
        "caller_id": 1,
        "gene_id": 1,
        "relative_path": "scripts/base_classes.py",
        "md5sum": "cef8890c1c8051d0c87919cf5e30fb54",
    }
    Queries.update_or_create(models.File, session, defaults={"id": 1}, **file_1)
    file_2 = {
        "caller_id": 1,
        "gene_id": 1,
        "relative_path": "scripts/__init__.py",
        "md5sum": "d41d8cd98f00b204e9800998ecf8427e",
    }
    Queries.update_or_create(models.File, session, defaults={"id": 2}, **file_2)
    session.commit()
