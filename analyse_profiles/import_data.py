'''
File to load data from the acquisitons and format it for analysis.

Author: Birk Emil Karlsen-Bæck
'''

# Import
import numpy as np
import matplotlib.pyplot as plt
import os


# Functions
def get_profile_data(f, fdir):
    data = np.load(fdir + f)
    return data

def find_file_in_folder(f, fdir):
    file_name = None
    for file in os.listdir(fdir):
        if file.startswith(f):
            file_name = file

    return file_name


