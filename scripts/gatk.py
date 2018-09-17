import glob
import json
import subprocess

class GATKCase:
    def __init__(self, sample_bam):
        # config
        max_cpu = 16
        max_mem = 16
        bam_dir = '/home/chris/stef/icr-mlpa/aligned'
        bed_file = '/home/chris/stef/bed-files/icr-filtered.bed'
        genome_fasta_path = '/var/reference_sequences/MiSeq/genome.fa'
        self.cromwell_jar_path = '/home/chris/stef/cromwell-32.jar'
        self.cohort_wdl_path = '/home/chris/stef/cnv-patissier/workflows/cnv_germline_case_workflow.wdl'

        # cohort job
        cromwell_dir = '/home/chris/stef/cromwell-executions'
        cohort_job_id = '7d5650f7-ed5d-4367-bbaf-cab25b8c1fca'
        model_prefix = 'ICR'


        self.sample_name = sample_bam.split("/")[-1].replace('.bam', '')


        self.case_json_data = {
            'CNVGermlineCaseWorkflow.bam': f'{bam_dir}/{sample_bam}',
            'CNVGermlineCaseWorkflow.bam_idx': f'{bam_dir}/{sample_bam}.bai',
            'CNVGermlineCaseWorkflow.ref_fasta': genome_fasta_path,
            'CNVGermlineCaseWorkflow.ref_fasta_fai': f'{genome_fasta_path}.fai',
            'CNVGermlineCaseWorkflow.ref_fasta_dict': genome_fasta_path.replace('.fa', '.dict'),
            'CNVGermlineCaseWorkflow.intervals': bed_file,
            'CNVGermlineCaseWorkflow.contig_ploidy_model_tar': (
                f'{cromwell_dir}/CNVGermlineCohortWorkflow/{cohort_job_id}/'
                f'call-DetermineGermlineContigPloidyCohortMode/execution/{model_prefix}-contig-ploidy-model.tar.gz'),
            'CNVGermlineCaseWorkflow.gcnv_model_tars':
                glob.glob(
                    f'{cromwell_dir}/CNVGermlineCohortWorkflow/{cohort_job_id}/'
                    f'call-GermlineCNVCallerCohortMode/shard-0/execution/{model_prefix}-model*.tar.gz'),
            'CNVGermlineCaseWorkflow.num_intervals_per_scatter': 5000,
            'CNVGermlineCaseWorkflow.padding': 50,
            'CNVGermlineCaseWorkflow.ref_copy_number_autosomal_contigs': 2,
            'CNVGermlineCaseWorkflow.cpu_for_determine_germline_contig_ploidy': max_cpu,
            'CNVGermlineCaseWorkflow.CollectCounts.cpu': max_cpu,
            'CNVGermlineCaseWorkflow.PreprocessIntervals.cpu': max_cpu,
            'CNVGermlineCaseWorkflow.cpu_for_germline_cnv_caller': max_cpu,
            'CNVGermlineCaseWorkflow.mem_gb_for_germline_cnv_caller': max_mem,
            'CNVGermlineCaseWorkflow.PreprocessIntervals.mem_gb': max_mem,
            'CNVGermlineCaseWorkflow.ScatterIntervals.mem_gb': max_mem,
            'CNVGermlineCaseWorkflow.PostprocessGermlineCNVCalls.mem_gb': max_mem,
            'CNVGermlineCaseWorkflow.gatk_docker': 'broadinstitute/gatk:4.0.7.0'
        }

    def remove_successful_run(self):
        """Remove cromwell execution and docker container"""
        pass

    def run_gcnv(self, json_path):
        """Run GATK gCNV CNVGermlineCaseWorkflow"""
        subprocess.run([
            'java', '-jar', self.cromwell_jar_path,
            'run', self.cohort_wdl_path,
            '-i', json_path],
        check=True)

    def save_cnv_calls(self):
        """Write CNV calls to database of results"""
        pass

    def main(self):
        """Create json, run gCNV and clean up if successful"""
        json_path = self.write_case_json()
        self.run_gcnv(json_path)

    def write_case_json(self):
        """Write case json data and returns the path"""
        output_path = f'{self.sample_name}-gCNV.json'
        with open(output_path, 'w') as out_file:
            json.dump(self.case_json_data, out_file)
        return output_path


if __name__ == '__main__':
    case_runner = GATKCase('17296.sorted.bam')
    case_runner.main()
    for bam in glob.glob('/home/chris/stef/icr-mlpa/aligned/*.bam'):
        case_runner = GATKCase(bam)
