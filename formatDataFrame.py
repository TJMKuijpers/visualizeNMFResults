import pandas as pd


class formatDataFrame:

    def __init__(self, data_frame_name):
        self.data_frame_name = data_frame_name
        self.data_frame = None
        self.data_subset = None

    def set_data_frame(self, data):
        self.data_frame = data

    def column_to_take(self, column_names):
        data_subset_to_columns=self.data_frame.loc[:, column_names]
        self.data_subset = data_subset_to_columns

    def rows_to_take(self, row_names):
        data_subset_to_rows = self.data_frame.loc[row_names,:]
        self.data_subset = data_subset_to_rows

    def rows_and_columns_to_take(self, column_names, row_names):
        data_subset = self.data_frame.loc[row_names, column_names]
        self.data_subset = data_subset


