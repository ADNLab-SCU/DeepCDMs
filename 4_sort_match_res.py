# -*- coding: utf-8 -*-

import os
import pandas as pd

def read_spectra_matched_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()
    match_result = content.split('\n\n')
    match_result.pop()
    return match_result

def get_cos_sim_res(line):
    parts = line.split()
    if parts[-2].replace(".", "", 1).isdigit():
        return float(parts[-2])
    return 0.0

def match_to_csv(match_result, write_path):
    match_res = []
    for match in match_result:
        lines = match.split("\n")
        sorted_data = sorted(lines[1:], key=get_cos_sim_res, reverse=True)
        try:    
            split_res = sorted_data[0].split()
            match_res.append([" ".join(i for i in split_res[:-3]), split_res[-2], split_res[-3], lines[0].strip("PEPMASS = ")])
        except:
            pass
    
    column_name = ["NAMES", "WCS", "SMILES", "EXACT_MASS"]
    match_res = pd.DataFrame(match_res, columns = column_name)
    match_res.to_csv(write_path, index=None)

def main():
    file = "_" # file of output of spectra_matching.py
    write_path = "_" # write out path
    
    match_result = read_spectra_matched_file(file)
    match_to_csv(match_result, write_path)

if __name__ == '__main__':
    main()