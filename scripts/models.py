from sqlalchemy import create_engine, ForeignKey, Column, Date, Integer, String, UniqueConstraint, Text, Time, sql
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker

from . import base_classes, utils

cnv_pat_dir = utils.get_cnv_patissier_dir()

# set up models
Base = declarative_base()


class CalledCNV(Base):
    __tablename__ = "called_cnvs"
    id = Column(Integer, primary_key=True)
    caller_id = Column(Integer, ForeignKey("callers.id", ondelete="CASCADE"), nullable=False)
    capture = Column(String)
    cnv_id = Column(Integer, ForeignKey("cnvs.id", ondelete="CASCADE"), nullable=False)
    gene_id = Column(Integer, ForeignKey("genes.id", ondelete="CASCADE"), nullable=False)
    json_data = Column(Text)
    sample_id = Column(Integer, ForeignKey("samples.id", ondelete="CASCADE"), nullable=False)

    def __repr__(self):
        return f"{self.id}"


class Caller(Base):
    __tablename__ = "callers"
    id = Column(Integer, primary_key=True)
    name = Column(String)

    def __repr__(self):
        return f"{self.name}"


class CNV(Base):
    __tablename__ = "cnvs"
    id = Column(Integer, primary_key=True)
    alt = Column(String)
    chrom = Column(String)
    end = Column(Integer)
    genome_build = Column(String)
    start = Column(Integer)

    def __repr__(self):
        return f"{self.build } {self.chrom}:{self.start}-{self.end} {self.alt}"


class File(Base):
    __tablename__ = "files"
    id = Column(Integer, primary_key=True)
    caller_id = Column(Integer, ForeignKey("callers.id", ondelete="CASCADE"), nullable=False)
    gene_id = Column(Integer, ForeignKey("genes.id", ondelete="CASCADE"), nullable=False)
    relative_path = Column(String)
    md5sum = Column(Text)

    __table_args__ = (UniqueConstraint(gene_id, relative_path),)

    def __repr__(self):
        return f"{self.relative_path}"


class Gene(Base):
    __tablename__ = "genes"
    id = Column(Integer, primary_key=True)
    genome_build = Column(String)
    capture = Column(String)
    chrom = Column(String)
    end = Column(Integer)
    start = Column(Integer)
    name = Column(String)

    __table_args__ = (UniqueConstraint(chrom, start, end, capture, genome_build),)

    def __repr__(self):
        return f"{self.name}"


class KnownCNV(Base):
    __tablename__ = "known_cnvs"
    id = Column(Integer, primary_key=True)
    cnv_id = Column(Integer, ForeignKey("cnvs.id", ondelete="CASCADE"), nullable=False)
    gene_id = Column(Integer, ForeignKey("genes.id", ondelete="CASCADE"), nullable=False)
    sample_id = Column(Integer, ForeignKey("samples.id", ondelete="CASCADE"), nullable=False)

    def __repr__(self):
        return f"{self.id}"


class Run(Base):
    __tablename__ = "runs"
    id = Column(Integer, primary_key=True)
    caller_id = Column(Integer, ForeignKey("callers.id", ondelete="CASCADE"), nullable=False)
    end_time = Column(Time)
    gene_id = Column(Integer, ForeignKey("genes.id", ondelete="CASCADE"), nullable=False)
    start_time = Column(Time)

    __table_args__ = (UniqueConstraint(caller_id, gene_id),)

    def __repr__(self):
        return f"{self.id}"


class Sample(Base):
    __tablename__ = "samples"
    id = Column(Integer, primary_key=True)
    bam_header = Column(Text)
    genome_build = Column(String)
    name = Column(String)
    path = Column(Text)
    result_type = Column(String)

    def __repr__(self):
        return f"{self.name}"


def get_or_create(model, defaults=None, **kwargs):  # could add session in arguments?
    params = {k: v for k, v in kwargs.items() if not isinstance(v, sql.ClauseElement)}
    params.update(defaults or {})
    instance = session.query(model).filter_by(**params).first()
    if instance:
        return instance, False
    else:
        instance = model(**params)
        session.add(instance)
        return instance, True


def update_or_create(model, defaults=None, **kwargs):
    params = {k: v for k, v in kwargs.items() if not isinstance(v, sql.ClauseElement)}
    params.update(defaults or {})
    query = session.query(model).filter_by(**defaults)
    if query.first():
        query.update(kwargs)
        instance = query.first()
        return instance, False
    else:
        instance = model(**params)
        session.add(instance)
        return instance, True


# setup db and session
engine = create_engine(f"sqlite:///{cnv_pat_dir}/output/results.db")
Base.metadata.create_all(engine)
Session = sessionmaker()
Session.configure(bind=engine)
session = Session()
