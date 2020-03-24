#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sample_ncbi_table.py: a helper script for downsampling of NCBI genome lists

@author: Tom C.C. HO 
"""
import random
import sys

if __name__ == '__main__':
    table_path = sys.argv[1]
    sample_size = int(sys.argv[2])
    random_seed = int(sys.argv[3])
    print('# Sampling from', table_path, 'a total of', sample_size, 'rows', file=sys.stderr)
    print('# Random seed:', random_seed, file=sys.stderr)
    output_array = []
    all_content_rows = []
    # read the genome list
    with open(table_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                output_array.append(line)
            elif line.rstrip('\r').rstrip('\n') == '':
                pass
            else:
                all_content_rows.append(line)
    # perform the sampling
    random.seed(random_seed)
    output_array += random.sample(all_content_rows, sample_size)
    # output the result
    for row in output_array:
        print(row, end='')
    

