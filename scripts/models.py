from sqlalchemy import create_engine, ForeignKey, Column, Date, Integer, String, UniqueConstraint, Text, Time
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker

from . import base_classes

# setup db and session
engine = create_engine(f"sqlite:///{base_classes.cnv_pat_dir}/output/results.db")

# set up models
Base = declarative_base()


class CalledCNV(Base):
    __tablename__ = "called_cnvs"
    id = Column(Integer, primary_key=True)
    caller_id = Column(Integer, ForeignKey("callers.id", ondelete='CASCADE'), nullable=False)
    capture = Column(String)
    cnv_id = Column(Integer, ForeignKey("cnvs.id", ondelete='CASCADE'), nullable=False)
    gene_id = Column(Integer, ForeignKey("genes.id", ondelete='CASCADE'), nullable=False)
    json_data = Column(Text)
    sample_id = Column(Integer, ForeignKey("samples.id", ondelete='CASCADE'), nullable=False)

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
    build = Column(String)
    chrom = Column(String)
    end = Column(Integer)
    start = Column(Integer)

    def __repr__(self):
        return f"{self.build }{self.chrom}:{self.start}-{self.end} {self.alt}"


class File(Base):
    __tablename__ = "files"
    id = Column(Integer, primary_key=True)
    caller_id = Column(Integer, ForeignKey("callers.id", ondelete='CASCADE'), nullable=False)
    gene_id = Column(Integer, ForeignKey("genes.id", ondelete='CASCADE'), nullable=False)
    relative_path = Column(String)
    md5sum = Column(Text)

    __table_args__ = (UniqueConstraint(gene_id, relative_path),)

    def __repr__(self):
        return f"{self.relative_path}"


class Gene(Base):
    __tablename__ = "genes"
    id = Column(Integer, primary_key=True)
    build = Column(String)
    cohort = Column(String)
    chrom = Column(String)
    end = Column(Integer)
    start = Column(Integer)
    name = Column(String)

    __table_args__ = (UniqueConstraint(chrom, start, end, cohort, build),)

    def __repr__(self):
        return f"{self.name}"


class KnownCNV(Base):
    __tablename__ = "known_cnvs"
    id = Column(Integer, primary_key=True)
    cnv_id = Column(Integer, ForeignKey("cnvs.id", ondelete='CASCADE'), nullable=False)
    gene_id = Column(Integer, ForeignKey("genes.id", ondelete='CASCADE'), nullable=False)
    sample_id = Column(Integer, ForeignKey("samples.id", ondelete='CASCADE'), nullable=False)

    def __repr__(self):
        return f"{self.id}"


class Run(Base):
    __tablename__ = "runs"
    id = Column(Integer, primary_key=True)
    caller_id = Column(Integer, ForeignKey("callers.id", ondelete='CASCADE'), nullable=False)
    end_time = Column(Time)
    gene_id = Column(Integer, ForeignKey("genes.id", ondelete='CASCADE'), nullable=False)
    start_time = Column(Time)

    __table_args__ = (UniqueConstraint(caller_id, gene_id),)

    def __repr__(self):
        return f"{self.id}"


class Sample(Base):
    __tablename__ = "samples"
    id = Column(Integer, primary_key=True)
    bam_header = Column(Text)
    build = Column(String)
    name = Column(String)
    path = Column(Text)

    def __repr__(self):
        return f"{self.name}"


Base.metadata.create_all(engine)
Session = sessionmaker()
Session.configure(bind=engine)
session = Session()
