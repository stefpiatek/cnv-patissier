from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base


Base = declarative_base()


class DbSession:
    factory = None
    engine = None

    @staticmethod
    def global_init(sqlite_path):
        if DbSession.factory:
            return

        if not sqlite_path or not sqlite_path.strip():
            raise Exception("You must specify a data file.")

        connection = f"sqlite:///{sqlite_path}"
        engine = create_engine(connection)
        DbSession.engine = engine
        DbSession.factory = sessionmaker(bind=engine)
        Base.metadata.create_all(engine)


# added at end to avoid circular import
from . import models
