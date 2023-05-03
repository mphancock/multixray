import unittest
import pandas as pd
import sys
from pathlib import Path

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from get_stat_df import get_stat_df


class TestGetStatDF(unittest.TestCase):

    def setUp(self):
        # Create sample input values for testing
        data_dir = Path(Path.home(), "xray/sample_bench/unit_tests/data")
        self.log_file_groups = [[Path(data_dir, "log_0.csv"), Path(data_dir, "log_1.csv")], [Path(data_dir, "log_2.csv")]]
        self.fields = ['xray', 'r_free', 'ff']
        self.stats = ['mean', 'std', 'min']
        self.N = 5
        self.offset = 5
        self.equil = 10

    def test_get_stat_df_returns_dataframe(self):
        # Ensure that the function returns a pandas dataframe
        stat_df = get_stat_df(self.log_file_groups, self.fields, self.stats, self.N, self.offset, self.equil)
        self.assertIsInstance(stat_df, pd.DataFrame)

    def test_get_stat_df_returns_expected_columns(self):
        # Ensure that the columns of the returned dataframe match the expected columns
        stat_df = get_stat_df(self.log_file_groups, self.fields, self.stats, self.N, self.offset, self.equil)
        columns = stat_df.columns.tolist()
        exp_cols = ['xray_mean', 'r_free_std']
        for i in range(self.N):
            exp_cols.append("ff_min_{}".format(i))
            exp_cols.append("ff_min_{}_log".format(i))
            exp_cols.append("ff_min_{}_id".format(i))
        self.assertEqual(columns, exp_cols)

    def test_get_stat_df_returns_expected_index(self):
        # Ensure that the index of the returned dataframe matches the expected index
        stat_df = get_stat_df(self.log_file_groups, self.fields, self.stats, self.N, self.offset, self.equil)
        index = stat_df.index.tolist()
        expected_index = [str(log_file_group) for log_file_group in self.log_file_groups]
        self.assertEqual(index, expected_index)

    def test_get_stat_df_raises_error_with_invalid_number_of_fields(self):
        # Ensure that an error is raised if the number of fields is not equal to the number of stats
        with self.assertRaises(RuntimeError):
            get_stat_df(self.log_file_groups, ['xray', 'r_free'], self.stats, self.N, self.offset, self.equil)

    def test_get_stat_df_with_invalid_stats(self):
        # Ensure that if an invalid stat is passed that the returned value is np.nan
        stat_df = get_stat_df(self.log_file_groups, self.fields, ['mean', 'min', 'invalid_stat'], self.N, self.offset, self.equil)
        for i in range(stat_df.shape[0]):
            self.assertTrue(pd.isna(stat_df["ff_invalid_stat"].iloc[i]))

    def test_get_stat_df_with_invalid_field(self):
        # Ensure that if an invalid field is passed that the returned value is np.nan
        stat_df = get_stat_df(self.log_file_groups, ['xray', 'invalid_field', 'r_free'], self.stats, self.N, self.offset, self.equil)
        for i in range(stat_df.shape[0]):
            self.assertTrue(pd.isna(stat_df["invalid_field_std"].iloc[i]))

    def test_get_stat_df_when_equil_greater_than_length(self):
        # Ensure that if the dataframe is empty or equil is greater than size of the dataframe that returned values are np.nan
        stat_df = get_stat_df(self.log_file_groups, self.fields, self.stats, self.N, self.offset, 1000)
        for i in range(stat_df.shape[0]):
            for j in range(stat_df.shape[1]):
                self.assertTrue(pd.isna(stat_df.iloc[i, j]))

    def test_get_stat_df_when_empty(self):
        # Ensure that if the passed file is an empty file that returned values are np.nan
        stat_df = get_stat_df([[Path(Path.home(), "xray/sample_bench/unit_tests/data/log_empty.csv")]], self.fields, self.stats, self.N, self.offset, self.equil)
        for i in range(stat_df.shape[0]):
            for j in range(stat_df.shape[1]):
                self.assertTrue(pd.isna(stat_df.iloc[i, j]))

    def test_get_stat_df(self):
        stat_df = get_stat_df(self.log_file_groups, self.fields, self.stats, self.N, 1, 0)
        self.assertAlmostEqual(stat_df["xray_mean"].iloc[0], 5.625173381290343)
        self.assertAlmostEqual(stat_df["r_free_std"].iloc[0], 0.02690603303822287)
        self.assertAlmostEqual(stat_df["ff_min_0"].iloc[0], 160.939455, delta=1e-5)
        self.assertAlmostEqual(stat_df["ff_min_1"].iloc[0], 161.779053, delta=1e-5)
        self.assertEqual(stat_df["ff_min_0_log"].iloc[0], Path(Path.home(), "xray/sample_bench/unit_tests/data/log_1.csv"))
        self.assertEqual(stat_df["ff_min_0_id"].iloc[0], 6)

    def test_get_stat_df_equil_and_offset(self):
        stat_df = get_stat_df(self.log_file_groups, self.fields, self.stats, self.N, self.offset, self.equil)
        self.assertAlmostEqual(stat_df["xray_mean"].iloc[1], 5.695093119734579)

if __name__ == '__main__':
    unittest.main()