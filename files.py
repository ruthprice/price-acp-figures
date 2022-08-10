# files.py
# Contains functions for manipulating files
# =====================================================================
import numpy as np
from glob import glob
import re
import iris
# =====================================================================
def get_csv_contents(filename, split_char=','):
    '''
    Open a file and return its contents in an object array
    Option to use different character for splitting rows
    default is a comma
    '''
    with open(filename, 'r') as file_reader:
        lines = file_reader.readlines()

    n_rows = len(lines)
    n_cols = len(lines[0].split(split_char))

    file_contents = np.empty([n_rows, n_cols], dtype=object)

    # put file contents into 2D list
    for i in np.arange(0, n_rows):
        for j in np.arange(0, n_cols):
            file_contents[i,j] = lines[i].strip().split(split_char)[j]

    return file_contents, n_rows, n_cols

def save_cube(cube, save_folder, save_name=None):
    """
    Saves cube as a netCDF file.
    """
    if save_name == None:
        if cube.standard_name:
            save_name=save_folder+'L1_'+cube.standard_name+'.nc'
        elif cube.long_name:
            save_name=save_folder+'L1_'+cube.long_name+'.nc'
        save_name = re.sub("\s+", "_", save_name)
    iris.save(cube, save_folder+save_name, netcdf_format="NETCDF4")
