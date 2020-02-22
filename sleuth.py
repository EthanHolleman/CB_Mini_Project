# File to handle running and parsing all sleuth operations
# Will need to run R scripts from this file since sleuth is done in R
# if graphs are a thing then need way to save those to png, jpeg etc.
import os
import csv


def default_conditions():
    return {'SRR5660030.1': '2dpi', 'SRR5660033.1': '2dpi',
            'SRR5660044.1': '6dpi', 'SRR5660045.1': '6dpi'}


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
            SRA_name = k_dir.split('_')[1].split('/')[1]
            if SRA_name in condition_dict:
                writer.writerow([SRA_name, condition_dict[SRA_name], k_dir])
            else:
                print(k_dir + ' not found making sleuth table')
    return sleuth_table