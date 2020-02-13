#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse_ncbi_table.py: a helper script for retrieving accession numbers from
                     NCBI Genome database

@author: Tom C.C. HO 
"""
import csv
import sys

def get_accessions(file_path, exclusion_string=None, flatten=True):
    """Return a list (of list) of accessions """
    print('\t# Loading accession list from', file_path + '...')
    accession_lists = list()
    with open(file_path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            if row[9] == 'Replicons':
                print('\t# Header line detected.', file=sys.stderr)
                continue
            if exclusion_string:
                if exclusion_string in row[0]:
                    continue
            accession_lists.append(parse_replicon_field(row[9]))
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