import unittest
import pandas as pd
import sys
from pathlib import Path

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from log_file_analysis import get_stat_df


class TestGetStatDF(unittest.TestCase):

    def setUp(self):
        # Create sample input values for testing
        data_dir = Path(Path.home(), "xray/sample_bench/unit_tests/data")
        self.log_file_groups = [[Path(data_dir, "log_0.csv"), Path(data_dir, "log_1.csv")], [Path(data_dir, "log_2.csv")]]
        self.fields = ['xray', 'r_free', 'ff']
        self.stats = ['mean', 'std', 'min']
        self.equil = 100
        self.include_id = True

    def test_get_stat_df_returns_dataframe(self):
        # Ensure that the function returns a pandas dataframe
        stat_df = get_stat_df(self.log_file_groups, self.fields, self.stats, self.equil, self.include_id)
        self.assertIsInstance(stat_df, pd.DataFrame)

    def test_get_stat_df_returns_expected_shape(self):
        # Ensure that the shape of the returned dataframe matches the expected shape
        stat_df = get_stat_df(self.log_file_groups, self.fields, self.stats, self.equil, self.include_id)
        expected_shape = (2, 5)
        self.assertEqual(stat_df.shape, expected_shape)

    def test_get_stat_df_with_no_id(self):
        # Ensure that the function returns a dataframe without ID columns when include_id is False
        stat_df = get_stat_df(self.log_file_groups, self.fields, self.stats, self.equil, include_id=False)
        columns = stat_df.columns.tolist()
        expected_columns = ['xray_mean', 'r_free_std', 'ff_min']
        self.assertEqual(columns, expected_columns)

    def test_get_stat_df_with_invalid_stats(self):
        # Ensure that the function raises a ValueError when an invalid stat is provided
        with self.assertRaises(RuntimeError):
            get_stat_df(self.log_file_groups, self.fields, ['mean', 'invalid_stat'], self.equil, self.include_id)

if __name__ == '__main__':
    unittest.main()