import numpy as np
import pandas as pd
import glob

# bitte alle funktionen, die verwendet werden sollen einzeln so importieren, nicht ganze dateien importieren
import functions as func




if __name__ == '__main__':
    all_files = glob.glob(r'./DMS_data/*.csv')

    func.load_data_frame(all_files)
