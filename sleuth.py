# File to handle running and parsing all sleuth operations
# Will need to run R scripts from this file since sleuth is done in R
# if graphs are a thing then need way to save those to png, jpeg etc.
import os
import csv


def default_conditions():
    return {'SRR5660030.1': '2dpi', 'SRR5660033.1': '2dpi',
            'SRR5660044.1': '6dpi', 'SRR5660045.1': '6dpi'}


def test_for_condition(condition_dict, path):
    for key in condition_dict:
        if key in path:
            return key, condition_dict[key]


def make_sleuth_table(kallisto_dirs, output_dir,
                      table_name='sleuth_table.csv', condition_dict=None):
    '''
    Given the list of output dirs created by individual kallisto runs
    creates an output csv table for sleuth to use. Names of file dirs are
    used for each entry.
    '''
    if condition_dict == None:
        condition_dict = default_conditions()

    sleuth_table = os.path.join(output_dir, table_name)
    with open(sleuth_table, 'w') as st:
        writer = csv.writer(st)
        writer.writerow(['sample', 'condition', 'path'])
        for k_dir in kallisto_dirs:
            sample, condition = test_for_condition(condition_dict, k_dir)
            writer.writerow([sample, condition, k_dir])
    return sleuth_table