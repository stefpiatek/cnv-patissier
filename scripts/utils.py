import csv
import os


def get_cnv_patissier_dir():
    """Returns the base directory of the project"""
    return os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


class SampleUtils:
    @classmethod
    def select_samples(cls, gene_list, normal_panel):
        """returns (sample_ids, sample_paths) from a gene's sample sheet"""
        if normal_panel:
            cnv_statuses = ["normal-panel"]
        else:
            cnv_statuses = ["normal", "positive"]
        output_ids = []
        output_paths = []
        with open(gene_list) as handle:
            samples = csv.DictReader(handle, delimiter="\t")
            for sample in samples:
                for cnv_status in cnv_statuses:
                    if sample["result_type"] == cnv_status:
                        output_ids.append(sample["sample_id"].strip())
                        output_paths.append(sample["sample_path"].strip())
        assert len(set(output_ids)) == len(output_ids), "sample sheet sample_ids must be unique"
        assert len(set(output_paths)) == len(
            output_paths
        ), "sample sheet sample_paths must be unique"
        return output_ids, output_paths

    @classmethod
    def get_mount_point(cls, paths):
        """Returns common root path for a list of paths"""
        abs_paths = [os.path.abspath(path) for path in paths]
        common_prefix = os.path.commonprefix(abs_paths)
        end_of_path = common_prefix.rfind("/")
        common_path = common_prefix[: end_of_path + 1]
        assert common_path, "all bams must be on the same drive/mount point"
        return common_path


cnv_pat_dir = get_cnv_patissier_dir()
