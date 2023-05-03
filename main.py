import numpy as np
import pandas as pd
import glob

# bitte alle funktionen, die verwendet werden sollen einzeln so importieren, nicht ganze dateien importieren
from functions.functions import load_data_frame






if __name__ == '__main__':
    all_files = glob.glob(r'./DMS_data/*.csv')

    load_data_frame(all_files)
