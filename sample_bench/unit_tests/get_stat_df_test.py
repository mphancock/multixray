import unittest
import pandas as pd
import sys
from pathlib import Path

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from get_stat_df import get_stat_df


class TestGetStatDF(unittest.TestCase):
    def setUp(self):
        # Create sample input values for testing.
        self.data_path = Path(Path.home(), "xray/sample_bench/data/test")
        self.log_file_groups = [[Path(self.data_path, "log_0.csv"), Path(self.data_path, "log_1.csv")], [Path(self.data_path, "log_2.csv")]]
        self.fields = ['xray', 'r_free', 'ff']
        self.stats = ['mean', 'std', 'min']
        self.N = 2
        self.offset = 1
        self.equil = 0

    def test_get_stat_df_returns_dataframe(self):
        # Ensure that the function returns a pandas dataframe.
        stat_df = get_stat_df(
            log_file_groups=self.log_file_groups,
            main_field="xray",
            main_stat="mean"
        )
        self.assertIsInstance(stat_df, pd.DataFrame)

    def test_get_stat_df_returns_expected_columns(self):
        # Ensure that the columns of the returned dataframe match the expected columns.
        stat_df = get_stat_df(
            log_file_groups=self.log_file_groups,
            main_field="xray",
            main_stat="max",
            N=2,
            bonus_fields=["pdb"]
        )
        columns = stat_df.columns.tolist()
        exp_cols = ["xray_max_0", "xray_max_0_pdb", "xray_max_1", "xray_max_1_pdb"]
        self.assertEqual(columns, exp_cols)

    def test_get_stat_df_returns_expected_index(self):
        # Ensure that the index of the returned dataframe matches the expected index.
        stat_df = get_stat_df(
            log_file_groups=self.log_file_groups,
            main_field="xray",
            main_stat="mean"
        )
        index = stat_df.index.tolist()
        expected_index = ["/wynton/home/sali/mhancock/xray/sample_bench/data/test/log_0.csv,/wynton/home/sali/mhancock/xray/sample_bench/data/test/log_1.csv", "/wynton/home/sali/mhancock/xray/sample_bench/data/test/log_2.csv"]
        self.assertEqual(index, expected_index)

    def test_get_stat_df_with_invalid_stats(self):
        # Ensure that an error is raised if an invalid stat is passed.
        with self.assertRaises(RuntimeError):
            get_stat_df(
                log_file_groups=self.log_file_groups,
                main_field="xray",
                main_stat="invalid"
            )

    def test_get_stat_df_with_invalid_field(self):
        # Ensure that an error is raised if an invalid field is passed.
        with self.assertRaises(RuntimeError):
            get_stat_df(
                log_file_groups=self.log_file_groups,
                main_field="invalid",
                main_stat="min"
            )

    def test_get_stat_df_when_empty(self):
        # Ensure that if the dataframe is empty or equil is greater than size of the dataframe that returned values are np.nan.
        stat_df = get_stat_df(
            log_file_groups=self.log_file_groups,
            main_field="xray",
            main_stat="mean",
            equil=100
        )
        for i in range(stat_df.shape[0]):
            for j in range(stat_df.shape[1]):
                self.assertTrue(pd.isna(stat_df.iloc[i, j]))

    def test_get_stat_df_when_corrupt(self):
        # Ensure that stat_df contains np.nan if the log file is empty.
        stat_df = get_stat_df(
            log_file_groups=[[Path(self.data_path, "corrupt.csv")]],
            main_field="xray",
            main_stat="mean"
        )
        self.assertTrue(pd.isna(stat_df["xray_mean"].iloc[0]))

    def test_get_stat_df(self):
        correct_stats = [5.625173381290343,0.02690603303822287]
        field_stats = [("xray", "mean"), ("r_free", "std")]
        for i in range(len(field_stats)):
            field, stat = field_stats[i]
            stat_df = get_stat_df(
                log_file_groups=self.log_file_groups,
                main_field=field,
                main_stat=stat
            )
            self.assertAlmostEqual(stat_df["{}_{}".format(field, stat)].iloc[0], correct_stats[i])

    def test_get_stat_df_multi(self):
        correct_stats = [179.245577,193.378256]
        correct_pdbs = ["/wynton/home/sali/mhancock/xray/sample_bench/data/test/pdbs/1.pdb", "/wynton/home/sali/mhancock/xray/sample_bench/data/test/pdbs/7.pdb"]
        stat_df = get_stat_df(
            log_file_groups=self.log_file_groups,
            main_field="ff",
            main_stat="min",
            N=2,
            bonus_fields=["pdb"],
            offset=10
        )
        for i in range(len(correct_stats)):
            self.assertAlmostEqual(stat_df["ff_min_{}".format(i)].iloc[1], correct_stats[i], delta=1e-5)
            self.assertEqual(stat_df["ff_min_{}_pdb".format(i)].iloc[1], correct_pdbs[i])

    def test_get_stat_df_fast(self):
        stat_df = get_stat_df(
            log_file_groups=[[Path(self.data_path, "log_0.csv"), Path(self.data_path, "log_1.csv"), Path(self.data_path, "log_2.csv")]],
            main_field="xray",
            main_stat="min",
            bonus_fields=["rmsd"],
            equil=10
        )
        self.assertAlmostEqual(stat_df["xray_min_0"].iloc[0], 5.589055, delta=1e-5)
        self.assertAlmostEqual(stat_df["xray_min_0_rmsd"].iloc[0], 0.119362, delta=1e-5)


if __name__ == '__main__':
    unittest.main()