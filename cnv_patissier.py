from argparse import ArgumentParser
import datetime
import pathlib

from scripts.db_session import DbSession
from scripts import utils, copywriter, codex2, cnv_kit, decon, excavator2, exome_depth, gatk, savvy_cnv, xhmm



def init_db(capture_name):
    db_path = f"{utils.get_cnv_patissier_dir()}/output/{capture_name}.sqlite"
    DbSession.global_init(db_path)

if __name__ == "__main__":
    parser = ArgumentParser(description="Ochrestrating your CNV-caller bakeoff")
    parser.add_argument("capture_name", help="After following the setup in the README.md, please give the capture name")
    args = parser.parse_args()

    start_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%m-%S")

    capture_name = args.capture_name 
    init_db(capture_name)
    cnv_pat_dir = utils.get_cnv_patissier_dir()
    sample_sheet_path = pathlib.Path(cnv_pat_dir, "input", capture_name, "sample-sheets")

    genes = [path.stem for path in list(sample_sheet_path.glob("*.txt"))]
    assert genes, "No genes found in path!"

    for gene in genes:
        # cnv_caller = cnv_kit.CNVKit(capture_name, gene, start_time)
        # cnv_caller.main()
        cnv_caller = codex2.CODEX2(capture_name, gene, start_time)
        cnv_caller.main()        
        cnv_caller = copywriter.Copywriter(capture_name, gene, start_time)
        cnv_caller.main()
        cnv_caller = decon.DECoN(capture_name, gene, start_time)
        cnv_caller.main()
        cnv_caller = excavator2.Excavator2(capture_name, gene, start_time)
        cnv_caller.main()
        cnv_caller = exome_depth.ExomeDepthCohort(capture_name, gene, start_time)
        cnv_caller.main()
        cnv_caller = exome_depth.ExomeDepthCase(capture_name, gene, start_time)
        cnv_caller.main()
        cnv_caller = gatk.GATKCohort(capture_name, gene, start_time)
        cnv_caller.main()
        cnv_caller = gatk.GATKCase(capture_name, gene, start_time)
        cnv_caller.main()
        cnv_caller = savvy_cnv.SavvyCNV(capture_name, gene, start_time)
        cnv_caller.main()
        cnv_caller = xhmm.XHMM(capture_name, gene, start_time)
        cnv_caller.main()


print("Congrats, you're all done")