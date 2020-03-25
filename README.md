# GolayMetaMiner: a software for k-mer based identification of clade-specific targets
This is the official repository of GolayMetaMiner. A manuscript on a specific diagnostic assay developed using this software is being submitted.

1.  Installing the software
  * The scripts are written in Python 3 and depend on a number of modules. The easiest way to install the required modules is by installing Anaconda. Anaconda may be downloaded from https://www.anaconda.com/distribution/
  * During installation, the Windows user is suggested to select installation of “Just Me” (as recommended), “Add Anaconda to my PATH environment variable” (despite not recommended) and “Register Anaconda as my default Python 3.x”.
  * The scripts require a working Internet connection to function properly. 
  * To download the software, you can click the download link of this repository, or simply
```
git clone https://github.com/hkhcc/GolayMetaMiner
```
  
2. Checking the installation
  * Check that Python is properly installed. In the Windows command processor (cmd.exe), type
```
python
```
and the Anaconda Python interpreter should be invoked:
```
C:\Users\COMPUTER_USER\Documents\Python Scripts\GolayMetaMiner>python
Python 3.7.4 (default, Aug  9 2019, 18:34:13) [MSC v.1915 64 bit (AMD64)] :: Anaconda, Inc. on win32

Warning:
This Python interpreter is in a conda environment, but the environment has
not been activated.  Libraries may fail to load.  To activate this environment
please see https://conda.io/activation

Type "help", "copyright", "credits" or "license" for more information.
>>>
```
You may then try to run the script gmm.py. Running the script with no arguments provided displays a (hopefully informative) help message, prompting for input. 
Note that, in Linux, ```python``` should probably be replace by ```python3``` or similar.
```
C:\Users\COMPUTER_USER\Documents\Python Scripts\GolayMetaMiner>python gmm.py
[    0.0000s]
C:\Users\COMPUTER_USER\Documents\Python Scripts\GolayMetaMiner\gmm-cache successfully created.
C:\Users\COMPUTER_USER\Documents\Python Scripts\GolayMetaMiner\gmm-cache ... OK!
C:\Users\COMPUTER_USER\Documents\Python Scripts\GolayMetaMiner\gmm-runs successfully created.
C:\Users\COMPUTER_USER\Documents\Python Scripts\GolayMetaMiner\gmm-runs ... OK!
[    0.0000s]
usage: gmm.py [-h] --primary_target PRIMARY_TARGET
              [--secondary_targets SECONDARY_TARGETS [SECONDARY_TARGETS ...] |
              --secondary_target_list SECONDARY_TARGET_LIST]
              [--reporting_centile REPORTING_CENTILE] [--kmer_size KMER_SIZE]
              [--num_threads NUM_THREADS] [--min_length MIN_LENGTH]
              [--non_targets NON_TARGETS [NON_TARGETS ...] |
              --non_target_ncbi_table NON_TARGET_NCBI_TABLE]
              [--exclusion_string EXCLUSION_STRING]
gmm.py: error: the following arguments are required: --primary_target
```
