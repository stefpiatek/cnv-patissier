from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base


Base = declarative_base()


class DbSession:
    factory = None
    engine = None
    connection = None

    @staticmethod
    def global_init(sqlite_path):
        if DbSession.factory:
            return

        if not sqlite_path or not sqlite_path.strip():
            raise Exception("No file given for database")

        connection = f"sqlite:///{sqlite_path}"
        engine = create_engine(connection)
        DbSession.engine = engine
        DbSession.factory = sessionmaker(bind=engine)
        DbSession.connection = Base.metadata.create_all(engine)

    @staticmethod
    def create_test_db(sqlite_path):
        connection = f"sqlite:///{sqlite_path}"
        engine = create_engine(connection)
        DbSession.engine = engine
        DbSession.connection = Base.metadata.create_all(engine)

    @staticmethod
    def create_test_session():
        DbSession.factory = sessionmaker(bind=DbSession.engine)


# added at end to avoid circular import
from . import models
