import pytest

from scripts import utils

cnv_pat_dir = utils.get_cnv_patissier_dir()


class TestSampleUtilsSelectSamples:
    def setup(self):
        self.sample_path_prefix = f"{cnv_pat_dir}/tests/test_files/input/cohort/sample-sheets/gene"

    def test_normal_panel(self):
        expected_ids = [str(number) for number in range(17327, 17331)]
        expected_paths = [f"/path/{number}.sorted.bam" for number in range(17327, 17331)]
        ids, paths = utils.SampleUtils.select_samples(f"{self.sample_path_prefix}_samples.txt", normal_panel=True)
        assert ids == expected_ids
        assert paths == expected_paths

    def test_unknown_cases(self):
        expected_ids = [str(number) for number in range(17331, 17333)]
        expected_paths = [f"/path/{number}.sorted.bam" for number in range(17331, 17333)]
        ids, paths = utils.SampleUtils.select_samples(f"{self.sample_path_prefix}_samples.txt", normal_panel=False)
        assert ids == expected_ids
        assert paths == expected_paths

    def test_dup_id(self):
        with pytest.raises(AssertionError):
            utils.SampleUtils.select_samples(f"{self.sample_path_prefix}-dup-sample-id_samples.txt", True)

    def test_dup_path(self):
        with pytest.raises(AssertionError):
            utils.SampleUtils.select_samples(f"{self.sample_path_prefix}-dup-sample-path_samples.txt", True)


class TestSampleUtilsGetMountPoint:
    def setup(self):
        self.paths = ["/mnt/data/runs/run1/bams/", "/mnt/data/runs/run2/bams/", "/mnt/data/validation/run1/bams/"]

    def test_simple(self):
        assert utils.SampleUtils.get_mount_point(self.paths) == "/mnt/data/"

    def test_different_mount(self):
        different_root = [*self.paths, "/home/user/this/wont/work/"]
        with pytest.raises(AssertionError):
            utils.SampleUtils.get_mount_point(different_root)
