import os
import subprocess
import csv

from data import if_not_dir_make


def run_sleuth(sleuth_table, output_dir, log, results_file_name='Sleuth_results.txt'):
    '''
    Given a sleuth table and an output dir run the sleuth_R.R script
    to execute differntial expression analysis via the sleuth R package.
    The sleuth table should be created from the make_sleuth_table function
    if it does not already exist. Returns path to the sleuth results.
    '''
    output_dir = if_not_dir_make(output_dir, 'sleuth_results')
    sleuth_results = os.path.join(output_dir, results_file_name)
    cmd = ['Rscript', 'sleuth_R.R', '-f', sleuth_table, '-o', sleuth_results]
    subprocess.call(cmd)

    write_results_to_log(sleuth_results, log)

    return sleuth_results


# need to see the format these come out in
def write_results_to_log(sleuth_results, log):
    '''
    Given the path to the sleuth results writes a header plus
    specific columns from each row of results to the log file.
    Does not return anything.
    '''
    HEADER = 'TARGET_id, test_stat, pval, qval\n'
    log_sleuth = []
    with open(sleuth_results) as sr:
        reader = csv.reader(sr, delimiter=' ')
        next(reader)  # skip the header row
        for row in reader:
            log_sleuth.append([row[0], row[3], row[1], row[2]])
    writer = csv.writer(log, delimiter='\t')
    log.write(HEADER)  # write the header row
    writer.writerows(log_sleuth)  # write all rows each row is represented as a sublist in log_sleuth


def default_conditions():
    '''
    Returns a dictionary of sample conditions based on the project URLs.
    2dpi = group1 and 6dpi = group2.
    '''
    return {'SRR5660030.1': '2dpi', 'SRR5660033.1': '6dpi',
            'SRR5660044.1': '2dpi', 'SRR5660045.1': '6dpi'}


def test_for_condition(condition_dict, path):
    '''
    Helper function to test if the string corresponding to a condition is
    present in a file path using the condition dictionary. This is to ensure
    proper assignment of condtions to kallisto results so they are not
    mixed up when handed off to sleuth in the form of the sleuth table.
    '''
    for key in condition_dict:  # iterate all keys see if in path
        if key in path:
            return key, condition_dict[key]  # return key + condition type


def make_sleuth_table(kallisto_dirs, output_dir,
                      table_name='sleuth_table.csv', condition_dict=None):
    '''
    Given the list of output dirs created by individual kallisto runs
    creates an output csv table for sleuth to use. Names of file dirs are
    used for each entry.
    '''
    if condition_dict == None:  # allows for providing alt condition dict
        condition_dict = default_conditions()

    sleuth_table = os.path.join(output_dir, table_name)
    with open(sleuth_table, 'w') as st:
        writer = csv.writer(st)
        writer.writerow(['sample', 'condition', 'path'])  # write the header
        for k_dir in kallisto_dirs:  # iterate through individual result dirs
            sample, condition = test_for_condition(condition_dict, k_dir)
            writer.writerow([sample, condition, k_dir])
    return sleuth_table
