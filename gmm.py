#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GolayMetaMiner: a software for k-mer based identification of clade-specific 
                targets

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
import numpy as np

from multiprocessing import Pool, Manager
from parse_ncbi_table import get_accessions

DEBUG = False
VERSION = '1.0 beta build 20200107'
CACHE_DIR = 'gmm-cache'
RESULT_DIR = 'gmm-runs'
KMER_SIZE = 12
PERCENTILE = 99.9
THREADS = 6


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
    print('\t', title, 'with', len(sequence), 'characters loaded.', file=sys.stderr)
    return title, sequence

def load_genome(accession, save=True, start=0, end=0, complementary_strand=False):
    """Load the genome (from cache) and return the parsed FASTA."""
    if end < start:
        start, end = end, start
        if not complementary_strand:
            complementary_strand = True
            if DEBUG:
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
        if DEBUG:
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
    if DEBUG:
        print('# Cleaned sequence length:', len(sequence), file=sys.stderr)
    return sequence

def pad_sequence(sequence, k=KMER_SIZE):
    """Return the padded sequence."""
    sequence = sequence + sequence[0:k-1]
    if DEBUG:
        print('# Padded sequence length for kmer generation:', len(sequence), file=sys.stderr)
    return sequence

def kmerize(sequence, k=KMER_SIZE, coerce_to='A', circular=True, both_strands=True):
    """Return the k-mer set of a genome sequence."""
    # pre-process the sequence
    sequence = clean_sequence(sequence)
    if DEBUG:
        if circular:
            print('# Topology: circular', file=sys.stderr)
        else:
            print('# Topology: linear', file=sys.stderr)
        if both_strands:
            print('# BOTH strands would be processed.', file=sys.stderr)
        else:
            print('# Only the input strand would be processed.', file=sys.stderr)
    # initialize the (unique) k-mer set
    kmers = dict()
    if circular:
        # pad the sequence
        orig_length = len(sequence)
        sequence = pad_sequence(sequence)
        for i in range(orig_length):
            kmer = sequence[i:i+k]
            if not kmer in kmers:
                kmers[kmer] = 1
        if both_strands:
            sequence = str(sequence).translate(str.maketrans('ATCG', 'TAGC'))[::-1]
            if DEBUG:
                print('# Reversed sequence length for kmer generation:', len(sequence), file=sys.stderr)
            for j in range(orig_length):
                kmer = sequence[j:j+k]
                if not kmer in kmers:
                    kmers[kmer] = 1
    else:
        if DEBUG:
            print('# Unpadded sequence length for kmer generation:', len(sequence), file=sys.stderr)
        for i in range(len(sequence)-k+1):
            kmer = sequence[i:i+k]
            if not kmer in kmers:
                kmers[kmer] = 1
        if both_strands:
            sequence = str(sequence).translate(str.maketrans('ATCG', 'TAGC'))[::-1]
            if DEBUG:
                print('# Reversed sequence length for kmer generation:', len(sequence), file=sys.stderr)
            for i in range(len(sequence)-k+1):
                kmer = sequence[j:j+k]
                if not kmer in kmers:
                    kmers[kmer] = 1
    if DEBUG:
        print('# Number of unique kmers counted:', len(kmers), file=sys.stderr)
    return kmers

def addto_kmer_pool(kmer_pool, accessions, k=KMER_SIZE):
    """Add k-mers from accessions to kmer-pool"""
    for accession in accessions:
        title, sequence = load_genome(accession)
        kmers = kmerize(sequence, k)
        kmer_pool.update(kmers)
        print('\tPool occupancy:', len(kmer_pool), '/', 4**KMER_SIZE, 
              '(' + str(round(100*len(kmer_pool)/4**KMER_SIZE, 2)) + '%)', file=sys.stderr)
        print('\t', sys.getsizeof(kmer_pool), file=sys.stderr)
        del kmers
        if DEBUG:
            print('# k = ', k, file=sys.stderr)
            print('# k-mer pool size:', len(kmer_pool), file=sys.stderr)

# check and create the genome cache and run result directories
check_create_dir(CACHE_DIR)
check_create_dir(RESULT_DIR)
      
# load the M. kansasii genome
print('## Loading target genome...', file=sys.stderr)
title, sequence = load_genome('NC_022663.1')
cleaned_sequence = clean_sequence(sequence)
padded_sequence = pad_sequence(cleaned_sequence)

# prepare the non-target genomes
print('## Loading non-target genome...', file=sys.stderr)
nt_accessions = get_accessions('mycobacterium_complete.csv', exclusion_string='kansasii')
task_list = list()
m = Manager()
nt_pool = m.dict()
for nt_accession in nt_accessions:
    task_list.append((nt_pool, [nt_accession,]))
p = Pool(THREADS)
p.starmap(addto_kmer_pool, task_list)

print('## Non-target k-mer pool generation finished.', file=sys.stderr)


# determine the uniqueness of the genome using a sliding-window approach
print('## Eliminating non-target k-mers...', file=sys.stderr)
uniqueness_array = np.zeros(len(sequence))
for i in range(len(sequence)):
    this_kmer = padded_sequence[i:i+KMER_SIZE]
    if this_kmer not in nt_pool:
        uniqueness_array[i] = 1

# manually release the memory
del nt_pool

# determine the conservedness of the genome (among the target species)
print('## Selecting conserved k-mers...', file=sys.stderr)
conservedness_array = np.zeros(len(sequence), dtype=np.int8)
t_accessions = ['CP019887.1', 'CP019886.1', 'CP019884.1',
                'CP019885.1', 'CP019883.1', 'CP019888.1']
for t_accession in t_accessions:
    t_pool = set()
    t_pool = addto_kmer_pool(t_pool, [t_accession,])
    for i in range(len(sequence)):
        this_kmer = padded_sequence[i:i+KMER_SIZE]
        if this_kmer in t_pool:
            conservedness_array[i] += 1
    del t_pool

u_array = scipy.signal.savgol_filter(uniqueness_array, 501, 3)
print('# Uniqueness (min):', np.min(u_array), file=sys.stderr)
print('# Uniqueness (50th centile):', np.percentile(u_array, 50), file=sys.stderr)
print('# Uniqueness (max):', np.max(u_array), file=sys.stderr)
c_array = scipy.signal.savgol_filter(conservedness_array/len(t_accessions), 501, 3)
print('# Conservedness (min):', np.min(c_array), file=sys.stderr)
print('# Conservedness (50th centile):', np.percentile(c_array, 50), file=sys.stderr)
print('# Conservedness (max):', np.max(c_array), file=sys.stderr)
plt.plot(u_array)
plt.plot(c_array)
plt.show()

