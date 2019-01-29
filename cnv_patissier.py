import glob

from scripts.db_session import DbSession
from scripts import utils, copywriter, cnv_kit, decon, excavator2, exome_depth, gatk, xhmm

start_time = "now"

def init_db():
    db_path = f"{utils.get_cnv_patissier_dir()}/output/results.sqlite"
    DbSession.global_init(db_path)

if __name__ == "__main__":
    init_db()
    cnv_pat_dir = utils.get_cnv_patissier_dir()

    for file in glob.glob(f"{cnv_pat_dir}/input/ICR/sample-sheets/*"):
        gene = file.replace("_samples.txt", "").split("/")[-1]
        cnv_caller = exome_depth.ExomeDepthCohort("ICR", gene, start_time)
        cnv_caller.main()
        cnv_caller = exome_depth.ExomeDepthTransitionCase("ICR", gene, start_time)
        cnv_caller.main()

    for file in glob.glob(f"{cnv_pat_dir}/input/ICR/sample-sheets/*"):
        gene = file.replace("_samples.txt", "").split("/")[-1]
        cnv_caller = cnv_kit.CNVKit("ICR", gene, start_time)
        cnv_caller.main()        
        cnv_caller = copywriter.Copywriter("ICR", gene, start_time)
        cnv_caller.main()
        cnv_caller = decon.DECoN("ICR", gene, start_time)
        cnv_caller.main()
        cnv_caller = excavator2.Excavator2("ICR", gene, start_time)
        cnv_caller.main()
        cnv_caller = exome_depth.ExomeDepthCohort("ICR", gene, start_time)
        cnv_caller.main()
        cnv_caller = exome_depth.ExomeDepthCase("ICR", gene, start_time)
        cnv_caller.main()
        cnv_caller = exome_depth.ExomeDepthTransitionCase("ICR", gene, start_time)
        cnv_caller.main()
        cnv_caller = gatk.GATKCohort("ICR", gene, start_time)
        cnv_caller.main()
        cnv_caller = gatk.GATKCase("ICR", gene, start_time)
        cnv_caller.main()
        cnv_caller = xhmm.XHMM("ICR", gene, start_time)
        cnv_caller.main()
    print("Congrats, you're all done")
