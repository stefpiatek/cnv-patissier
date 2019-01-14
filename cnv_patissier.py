from scripts.db_session import DbSession


def init_db():
    db_path = f"{utils.get_cnv_patissier_dir()}/output/results.sqlite"
    DbSession.global_init(db_path)

if __name__ == "__main__":
    init_db()
    # now run the caller scripts
