#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GolayMetaMiner: a software for k-mer based identification of clade-specific 
                targets

@author: Tom C.C. HO 
"""
import argparse
import hashlib
import os
import re
import sys
import textwrap
import time

import scipy.signal
import urllib.parse
import urllib.request

import matplotlib.pyplot as plt
import numpy as np

from multiprocessing import Pool, Manager, cpu_count
from threading import RLock
from parse_ncbi_table import get_accessions, parse_simple_acclist

# initialize the timer after major imports
from datetime import datetime
t0 = datetime.now().timestamp()

DEBUG = False
VERSION = '1.0 codebase 20200301'
CACHE_DIR = 'gmm-cache'
RESULT_DIR = 'gmm-runs'
PERCENTILE = 99.99
KMER_SIZE = 12
THREADS = cpu_count()
MIN_LEN = 50

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

def addto_kmer_pool(kmer_pool, accessions, k):
    """Add k-mers from accessions to kmer-pool"""
    lock = RLock()
    for accession in accessions:
        title, sequence = load_genome(accession)
        kmers = kmerize(sequence, k)        
        lock.acquire()
        kmer_pool.update(kmers)
        print('\tPool occupancy:', len(kmer_pool), '/', 4**KMER_SIZE, 
              '(' + str(round(100*len(kmer_pool)/4**KMER_SIZE, 2)) + '%)', file=sys.stderr)
        print('\t', sys.getsizeof(kmer_pool), file=sys.stderr)
        del kmers
        lock.release()
        
def detect_roi(u_array, c_array, u_cutoff, min_length=50):
    """Report u_array regions that are above u_cutoff"""
    roi = list()
    in_region = False
    base_pos = 1
    for u_score, c_score in zip(u_array, c_array):
        if in_region == False and u_score >= u_cutoff:
            in_region = True # turn on recording
            roi.append([base_pos, 0])
        elif in_region == True and u_score >= u_cutoff:
            pass
        elif in_region == True and u_score < u_cutoff:
            in_region = False # turn off recording
            roi[-1][1] = base_pos
        else:
            pass
        base_pos += 1
    len_filtered_roi = list()
    for region in roi:
        if (region[1] - region[0] + 1) >= min_length:
            len_filtered_roi.append(region)
    return len_filtered_roi

def print_time():
    """Print the time since script start"""
    print('[{time:10.4f}s]'.format(time=datetime.now().timestamp()-t0), file=sys.stderr)

if __name__ == "__main__":
    
    # check and create the genome cache and run result directories
    print_time()
    check_create_dir(CACHE_DIR)
    check_create_dir(RESULT_DIR)
    
    # parse the command line arguments
    print_time()
    parser = argparse.ArgumentParser(description='GolayMetaMiner: a software for k-mer based identification of clade-specific targets')
    parser.add_argument('--primary_target', required=True, type=str, help='Accession number of the primary target genome.')
    secondary_target_sources = parser.add_mutually_exclusive_group()
    secondary_target_sources.add_argument('--secondary_targets', type=str, nargs='+', help='Accession number(s) of the secondary target genomes.')
    secondary_target_sources.add_argument('--secondary_target_list', type=str, help='Path to a TXT list of secondary target accession numbers.')
    parser.add_argument('--reporting_centile', type=float, default=PERCENTILE, help='Specific target centile score threshold (default: ' + str(PERCENTILE) + ').')
    parser.add_argument('--kmer_size', type=int, default=KMER_SIZE, help='Genome analysis k-mer size (default: ' + str(KMER_SIZE)+ ').')
    parser.add_argument('--num_threads', type=int, default=THREADS, help='Number of CPU threads to use (default: ' + str(THREADS) + ').')
    parser.add_argument('--min_length', type=int, default=MIN_LEN, help='Minimum ROI length (default: ' + str(MIN_LEN) + ').')
    non_target_sources = parser.add_mutually_exclusive_group()
    non_target_sources.add_argument('--non_targets', type=str, nargs='+', help='Accession number(s) of the non-target genomes.')
    non_target_sources.add_argument('--non_target_ncbi_table', type=str, help='Path to a CSV table of non-target genomes from https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/')
    parser.add_argument('--exclusion_string', type=str, help='Species-name exclusion for table-based input of non-target genomes.')
    args = parser.parse_args()
    
    # Early checking of some essential arguments
    if not args.secondary_targets and not args.secondary_target_list:
        raise ValueError('Please specify either --secondary_targets or --secondary_target_list')
    if not args.non_target_ncbi_table and not args.non_targets:
        raise ValueError('Please specify either --non_targets or --non_target_ncbi_table')
    
    # Change analysis settings (if necessary)
    if args.reporting_centile:
        PERCENTILE = args.reporting_centile
    if args.kmer_size:
        KMER_SIZE = args.kmer_size
    if args.num_threads:
        THREADS = args.num_threads
    if args.min_length:
        MIN_LEN = args.min_length
    
    # Print the current settings
    print('====================', file=sys.stderr)
    print_time()
    print('ANALYSIS SETTINGS:', file=sys.stderr)
    print('Reporting centile cutoff: {}'.format(PERCENTILE), file=sys.stderr)
    print('k-mer size: {}'.format(KMER_SIZE), file=sys.stderr)
    print('Number of CPU threads to use: {}'.format(THREADS), file=sys.stderr)
    print('Minimum target region length to report: {}'.format(MIN_LEN), file=sys.stderr)
    print('====================', file=sys.stderr)

    # load the primary target genome
    print('## Loading target genome...', file=sys.stderr)
    print_time()
    title, sequence = load_genome(args.primary_target)
    cleaned_sequence = clean_sequence(sequence)
    padded_sequence = pad_sequence(cleaned_sequence)
    
    # prepare the non-target genomes
    print('## Loading non-target genomes...', file=sys.stderr)
    print_time()
    if args.non_target_ncbi_table:
        nt_accessions = get_accessions(args.non_target_ncbi_table, exclusion_string=args.exclusion_string)
    else:
        nt_accessions = args.non_targets
    task_list = list()
    m = Manager()
    nt_pool = m.dict()
    for nt_accession in nt_accessions:
        task_list.append((nt_pool, [nt_accession,], KMER_SIZE))
    p = Pool(THREADS)
    p.starmap(addto_kmer_pool, task_list)
    
    print('## Non-target k-mer pool generation finished.', file=sys.stderr)
    print_time()

    # this step serves to bypass a Windows-specific bug
    print('\t', type(nt_pool), len(nt_pool), sys.getsizeof(nt_pool), file=sys.stderr)
    print('\tConverting shared dict() object to set()...', file=sys.stderr)
    nt_pool = set(nt_pool.keys())
    print('\t', type(nt_pool), len(nt_pool), sys.getsizeof(nt_pool), file=sys.stderr)

    # determine the uniqueness of the genome using a sliding-window approach
    print('## Eliminating non-target k-mers...', file=sys.stderr)
    print_time()
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
    if args.secondary_target_list:
        t_accessions = parse_simple_acclist(args.secondary_target_list)
    else:
        t_accessions = args.secondary_targets
    for t_accession in t_accessions:
        print_time()
        t_pool = dict()
        addto_kmer_pool(t_pool, [t_accession,], KMER_SIZE)
        for i in range(len(sequence)):
            this_kmer = padded_sequence[i:i+KMER_SIZE]
            if this_kmer in t_pool:
                conservedness_array[i] += 1
        del t_pool
    
    u_array = scipy.signal.savgol_filter(uniqueness_array, 501, 3)
    u_cutoff = np.percentile(u_array, PERCENTILE)
    print('# Uniqueness (min):', np.min(u_array), file=sys.stderr)
    print('# Uniqueness (50th centile):', np.percentile(u_array, 50), file=sys.stderr)
    print('# Uniqueness cutoff score (i.e. ' + str(PERCENTILE) + 'th centile:', u_cutoff, file=sys.stderr)
    print('# Uniqueness (max):', np.max(u_array), file=sys.stderr)
    c_array = scipy.signal.savgol_filter(conservedness_array/len(t_accessions), 501, 3)
    print('# Conservedness (min):', np.min(c_array), file=sys.stderr)
    print('# Conservedness (50th centile):', np.percentile(c_array, 50), file=sys.stderr)
    print('# Conservedness (max):', np.max(c_array), file=sys.stderr)
    print_time()
    
    # report the regions of interest
    print('## Reporting ROIs...', file=sys.stderr)
    print_time()
    roi_list = detect_roi(u_array, c_array, u_cutoff, min_length=MIN_LEN)
    for roi in roi_list:
        print('# start:', roi[0], 'end:', roi[1], 'length:', roi[1] - roi[0] + 1, file=sys.stderr)
        print('>{title} ({start}-{end}, len {length})'.format(title=title, start=roi[0], end=roi[1], length=roi[1]-roi[0]+1))
        roi_sequence = sequence[roi[0]-1:roi[1]] # roi[0] and roi[1] are 1-based 
        seq_lines = textwrap.wrap(roi_sequence, width=50)
        for line in seq_lines:
            print(line)
    # plot the uniqueness and conservedness 
    plt.plot(u_array)
    plt.axhline(u_cutoff, c='red')
    plt.plot(c_array)
    plt.title('{genome_title} genome plot'.format(genome_title=str.upper(title)), fontweight='bold')
    plt.legend(('Uniqueness score (U-score)', 'U-score cutoff at ' + str(PERCENTILE) + 'th percentile', 'Conservedness score (C-score)'))
    plt.xlabel('Genome position (bp)')
    plt.ylabel('Score')
    print('## Finished plotting. Rendering plot...', file=sys.stderr)
    print_time()
    plt.show()

    
    
