#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse_ncbi_table.py: a helper script for retrieving accession numbers from
                     NCBI Genome database

@author: Tom C.C. HO 
"""
import csv
import os
import sys

def get_accessions(file_path, exclusion_string=None, flatten=True):
    """Return a list (of list) of accessions """
    assert os.path.isfile(file_path)
    print('\t# Loading accession list from', file_path + '...', file=sys.stderr)
    print('\t# Excluding species names containing "' + exclusion_string +'"', file=sys.stderr)
    accession_lists = list()
    with open(file_path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        replicons_col = None
        for row in csv_reader:
            if row[0].startswith('#'):
                print('\t# Presumptive header line detected.', file=sys.stderr)
                col_index = 0
                for cell in row:
                    if cell == 'Replicons':
                        replicons_col = col_index
                        print('\t# Replicons column detected at array index ' + str(col_index), file=sys.stderr)
                    col_index += 1
                continue
            if exclusion_string:
                if exclusion_string in row[0]:
                    continue
            accession_lists.append(parse_replicon_field(row[replicons_col]))
    if flatten:
        flat_list = list()
        for accession_list in accession_lists:
            for accession in accession_list:
                flat_list.append(accession)
        return flat_list
    return accession_list

def parse_replicon_field(accession_string):
    """Return a list of accessions"""
    accessions = list()
    replicons = accession_string.split('; ')
    for replicon in replicons:
        replicon_accession = replicon.split('/')[-1]
        replicon_accession = replicon_accession.split(':')[-1]
        accessions.append(replicon_accession)
    return accessions

def parse_simple_acclist(file_path):
    """Return a list of accessions"""
    assert os.path.isfile(file_path)
    accessions = list()
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.rstrip()
                if not line == '':
                    accessions.append(line)
    return accessions