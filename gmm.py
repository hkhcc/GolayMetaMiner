#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GolayMetaMiner: a software for k-mer based identification of clade-specific 
                targets

Created on Tue Jan  7 21:43:46 2020
@author: Tom C.C. HO 
"""
import hashlib
import os
import re
import sys
import time

import scipy.signal
import urllib.parse
import urllib.request

import matplotlib.pyplot as plt


DEBUG = False
VERSION = '1.0 beta build 20200107'
CACHE_DIR = 'gmm-cache'
RESULT_DIR = 'gmm-runs'
KMER_SIZE = 16

def check_create_dir(directory):
    """Create a storage directory if not already present."""
    dir_path = os.path.join(os.path.abspath(os.path.curdir), directory)
    if not os.path.isdir(dir_path):
        try:
            os.mkdir(dir_path)
        except:
            raise PermissionError(dir_path, 'could not be created!')
        print(dir_path, 'successfully created.', file=sys.stderr)
    print(dir_path, '... OK!', file=sys.stderr)

def wget(url):
    """Return plain text response from a URL"""
    d = None
    passed = False
    while not passed:
        try:
            with urllib.request.urlopen(url) as response:
                d = response.read().decode('ascii')
            passed = True
        except:
            print('# Waiting...', file=sys.stderr)
            time.sleep(5)
    return d

def load_fasta(fasta_path):
    """Return the sequence title and data from a single-sequence FASTA file."""
    file_content = ''
    with open(fasta_path) as f:
        file_content = f.read()
    file_content = file_content.replace('\r', '')
    title, sequence = file_content.split(sep='\n', maxsplit=1)
    # clean up the title and sequence
    title = title.replace('>', '')
    sequence = sequence.replace('\n', '')
    print(title, 'with', len(sequence), 'characters loaded.', file=sys.stderr)
    return title, sequence

def load_genome(accession, save=True, start=0, end=0, complementary_strand=False):
    """Load the genome (from cache) and return the parsed FASTA."""
    if end < start:
        start, end = end, start
        if not complementary_strand:
            complementary_strand = True
            print('# Coerced to complementary strand', file=sys.stderr)
    if start == 0 and end == 0 and complementary_strand == False:
        genome_file = hashlib.sha1(accession.encode('ascii')).hexdigest() + '.fasta'
    else:
        tag = str(start) + '-' + str(end)
        if complementary_strand:
            tag = 'c' + tag
        genome_file = hashlib.sha1(accession.encode('ascii')).hexdigest() + '.' + tag + '.fasta'
    genome_file_path = os.path.join(os.path.abspath(os.path.curdir), CACHE_DIR, genome_file)
    if os.path.isfile(genome_file_path):
        pass
    else:
        efetch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text'
        efetch_url += '&id=' + accession
        if start != 0 and end != 0:
            efetch_url += '&seq_start=' + str(start) + '&seq_stop=' + str(end)
            if complementary_strand:
                efetch_url += '&strand=2'
        print(efetch_url, file=sys.stderr)
        fasta_data = wget(efetch_url)
        try:
            with open(genome_file_path, 'w') as f:
                print(fasta_data, file=f)
        except:
            raise IOError(genome_file_path, 'could not be created!')
    return load_fasta(genome_file_path)

def clean_sequence(sequence, coerce_to='A'):
    """Return the cleaned sequence."""
    sequence = str.strip(sequence)
    sequence = str.upper(sequence)
    sequence = re.sub('[^ATCG]', coerce_to, sequence)
    print('# Cleaned sequence length:', len(sequence), file=sys.stderr)
    return sequence

def pad_sequence(sequence, k=KMER_SIZE):
    """Return the padded sequence."""
    sequence = sequence + sequence[0:k-1]
    print('# Padded sequence length for kmer generation:', len(sequence), file=sys.stderr)
    return sequence

def kmerize(sequence, k=KMER_SIZE, coerce_to='A', circular=True, both_strands=True):
    """Return the k-mer dictionary of a genome sequence."""
    # pre-process the sequence
    sequence = clean_sequence(sequence)
    if circular:
        print('# Topology: circular', file=sys.stderr)
    else:
        print('# Topology: linear', file=sys.stderr)
    if both_strands:
        print('# BOTH strands would be processed.', file=sys.stderr)
    else:
        print('# Only the input strand would be processed.', file=sys.stderr)
    # initialize the k-mer dictionary
    kmers = dict()
    if circular:
        # pad the sequence
        orig_length = len(sequence)
        sequence = pad_sequence(sequence)
        for i in range(orig_length):
            kmer = sequence[i:i+k]
            if kmer not in kmers:
                kmers[kmer] = list()
            kmers[kmer].append(i)
        if both_strands:
            sequence = str(sequence).translate(str.maketrans('ATCG', 'TAGC'))[::-1]
            print('# Reversed sequence length for kmer generation:', len(sequence), file=sys.stderr)
            for j in range(orig_length):
                kmer = sequence[j:j+k]
                if kmer not in kmers:
                    kmers[kmer] = list()
                kmers[kmer].append(-j-1)
    else:
        print('# Unpadded sequence length for kmer generation:', len(sequence), file=sys.stderr)
        for i in range(len(sequence)-k+1):
            kmer = sequence[i:i+k]
            if kmer not in kmers:
                kmers[kmer] = list()
            kmers[kmer].append(i)
        if both_strands:
            sequence = str(sequence).translate(str.maketrans('ATCG', 'TAGC'))[::-1]
            print('# Reversed sequence length for kmer generation:', len(sequence), file=sys.stderr)
            for i in range(len(sequence)-k+1):
                kmer = sequence[j:j+k]
                if kmer not in kmers:
                    kmers[kmer] = list()
                kmers[kmer].append(-j-1)
    print('# Number of unique kmers counted:', len(kmers), file=sys.stderr)
    return kmers

# check and create the genome cache and run result directories
check_create_dir(CACHE_DIR)
check_create_dir(RESULT_DIR)

# prepare the non-target genomes
nt_pool = dict()
nt_accessions = ['CP022546.1', 'LR130759.1']

for accession in nt_accessions:
    nt_title, nt_sequence = load_genome(accession)
    nt_kmers = kmerize(nt_sequence)
    nt_pool = dict(nt_pool, **nt_kmers)
    print('# Non-target pool size:', len(nt_pool), file=sys.stderr)

# load the M. kansasii genome
title, sequence = load_genome('NC_022663.1')
cleaned_sequence = clean_sequence(sequence)
padded_sequence = pad_sequence(cleaned_sequence)

uniqueness_array = []
for i in range(len(sequence)):
    this_kmer = padded_sequence[i:i+KMER_SIZE]
    if this_kmer not in nt_pool:
        uniqueness_array.append(1)
    else:
        uniqueness_array.append(0)

# s_array = scipy.signal.savgol_filter(uniqueness_array, 501, 3)
s_array = uniqueness_array

plt.plot(s_array)
plt.hlines(0.5, 1411362, 1411475)
plt.hlines(0.5, 1411546, 1414663)
plt.hlines(0.5, 1414932, 1416471)
plt.show()

