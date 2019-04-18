import csv
import os
import pathlib


def get_cnv_patissier_dir():
    """Returns the base directory of the project"""
    return os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


class SampleUtils:
    @classmethod
    def check_files(cls, paths):
        """ Takes in a list of paths and raises Exception if:
         - File doesn't exist
         - Has invalid character in filename
         - Paths are not unique

        """
        cls.check_unique(paths, "sample_path")

        files = [pathlib.Path(path).resolve() for path in paths]
        for file in files:
            # CODEX2 automatically replaces '-' with '.'
            if "-" in file.name:
                raise Exception(
                    f"File {file} has a '-' in which is not allowed, please rename (or make a temporary copy of) "
                    "this file and the bam index, then update sample sheet field"
                )
            if not file.exists():
                raise Exception(f"File {file} does not exist")

    @classmethod
    def check_unique(cls, items, data_type):
        """If items are not unique, returns exception with duplicate items listed """
        non_unique = set(item for item in items if items.count(item) > 1)
        if non_unique:
            non_unique_out = "\n ".join(non_unique)
            raise Exception(f"The the following {data_type}(s) are not unique:\n {non_unique_out}")

    @classmethod
    def get_bam_to_id(cls, sample_sheet):
        """Returns dictionary of bam paths to sample ids from the sample sheet """
        normal_id, normal_path = cls.select_samples(sample_sheet, normal_panel=True)
        unknown_id, unknown_path = cls.select_samples(sample_sheet, normal_panel=False)
        paths = normal_path + unknown_path
        sample_ids = normal_id + unknown_id
        return {path: sample_id for (path, sample_id) in zip(paths, sample_ids)}

    @classmethod
    def select_samples(cls, sample_sheet, normal_panel):
        """returns (sample_ids, sample_paths) from a gene's sample sheet"""
        if normal_panel:
            cnv_statuses = ["normal-panel"]
        else:
            cnv_statuses = ["normal", "positive"]
        output_ids = []
        output_paths = []
        with open(sample_sheet) as handle:
            samples = csv.DictReader(handle, delimiter="\t")
            for sample in samples:
                if sample["sample_id"] and sample["sample_path"]:
                    for cnv_status in cnv_statuses:
                        if sample["result_type"] == cnv_status:
                            output_ids.append(sample["sample_id"].strip())
                            output_paths.append(sample["sample_path"].strip())
        if normal_panel:
            assert len(output_ids) >= 30, "There must be 30 normal-panel samples in the sample sheet"
        return output_ids, output_paths

    @classmethod
    def get_mount_point(cls, paths):
        """Returns common root path for a list of paths"""
        abs_paths = [os.path.abspath(path) for path in paths]
        common_prefix = os.path.commonprefix(abs_paths)
        end_of_path = common_prefix.rfind("/")
        common_path = common_prefix[: end_of_path + 1]
        assert common_path != "/", "all bams must be on the same drive/mount point"
        return common_path


cnv_pat_dir = get_cnv_patissier_dir()
