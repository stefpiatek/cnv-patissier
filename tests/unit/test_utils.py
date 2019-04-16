import pathlib

import pytest

from scripts import utils

cnv_pat_dir = utils.get_cnv_patissier_dir()


class TestCheckFiles:
    def setup(self):
        self.test_file_prefix = f"{cnv_pat_dir}/tests/test_files/input/checks"

        for sample in range(40):
            with open(f"{self.test_file_prefix}/sample{sample}.txt", "w") as handle:
                handle.write("dummy")

    def test_no_issues(self):
        paths = [pathlib.Path(f"{self.test_file_prefix}/sample{sample}.txt") for sample in range(40)]
        utils.SampleUtils.check_files(paths)

    def test_dash_in_name(self):
        paths = [pathlib.Path(f"{self.test_file_prefix}/sample-{sample}.txt") for sample in range(1)]
        with pytest.raises(Exception):
            utils.SampleUtils.check_files(paths)

    def test_file_doesnt_exist(self):
        paths = [pathlib.Path(f"{self.test_file_prefix}/sample{sample}.txt") for sample in range(50, 60)]
        with pytest.raises(Exception):
            utils.SampleUtils.check_files(paths)

    def test_not_unique(self):
        path_1 = [pathlib.Path(f"{self.test_file_prefix}/sample{sample}.txt") for sample in range(30)]
        path_2 = [pathlib.Path(f"{self.test_file_prefix}/sample{sample}.txt") for sample in range(10, 12)]
        with pytest.raises(Exception):
            utils.SampleUtils.check_files([*path_1, *path_2])


class TestCheckUnique:
    def test_no_issues(self):
        samples = [f"sample_{number}" for number in range(50)]
        utils.SampleUtils.check_unique(samples, "samples")

    def test_not_unique(self):
        sample_1 = [f"sample_{number}" for number in range(50)]
        sample_2 = [f"sample_{number}" for number in range(10, 12)]
        with pytest.raises(Exception):
            utils.SampleUtils.check_files([*sample_1, *sample_2], "samples")


class TestSampleUtilsSelectSamples:
    def setup(self):
        self.sample_path_prefix = f"{cnv_pat_dir}/tests/test_files/input/capture/sample-sheets/gene"

    def test_normal_panel(self):
        expected_ids = [str(number) for number in range(0, 30)]
        expected_paths = [f"/path/{number}.sorted.bam" for number in range(0, 30)]
        ids, paths = utils.SampleUtils.select_samples(f"{self.sample_path_prefix}_samples.txt", normal_panel=True)
        assert ids == expected_ids
        assert paths == expected_paths

    def test_unknown_cases(self):
        expected_ids = [str(number) for number in range(30, 32)]
        expected_paths = [f"/path/{number}.sorted.bam" for number in range(30, 32)]
        ids, paths = utils.SampleUtils.select_samples(f"{self.sample_path_prefix}_samples.txt", normal_panel=False)
        assert ids == expected_ids
        assert paths == expected_paths

    def test_lt_30_normal_panel(self):
        with pytest.raises(AssertionError):
            utils.SampleUtils.select_samples(f"{self.sample_path_prefix}-lt_30.txt", True)


class TestSampleUtilsGetMountPoint:
    def setup(self):
        self.paths = ["/mnt/data/runs/run1/bams/", "/mnt/data/runs/run2/bams/", "/mnt/data/validation/run1/bams/"]

    def test_simple(self):
        assert utils.SampleUtils.get_mount_point(self.paths) == "/mnt/data/"

    def test_different_mount(self):
        different_root = [*self.paths, "/home/user/this/wont/work/"]
        with pytest.raises(AssertionError):
            utils.SampleUtils.get_mount_point(different_root)
