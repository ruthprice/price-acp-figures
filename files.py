# files.py
# Contains functions for manipulating files
# =====================================================================
import numpy as np
# =====================================================================
def get_csv_contents(filename, split_char=','):
    # Open a file and return its contents in an object array
    # Option to use different character for splitting rows
    # default is a comma

    with open(filename, 'r') as file_reader:
        lines = file_reader.readlines()

    n_rows = len(lines)
    n_cols = len(lines[0].split(split_char))

    file_contents = np.empty([n_rows, n_cols], dtype=object)

    # put file contents into 2D list
    for i in np.arange(0, n_rows):
        for j in np.arange(0, n_cols):
            file_contents[i,j] = lines[i].strip().split(split_char)[j]

    return [file_contents, n_rows, n_cols]
