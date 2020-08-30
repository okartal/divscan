# divscan - Jensen-Shannon divergence for genomics

Divscan combines a set of python functions and a command-line script to perform genomic scans for Jensen-Shannon divergence over a given set of BED-like genome position files (GPFs). 

## How to run divscan in a conda environment

### Requirements

- POSIX-like operating systems
- git (<https://git-scm.com/>)
- miniconda3 (<https://conda.io/miniconda.html>), optional: Anaconda3

### Installation

1. Clone the repository.
2. Ensure that bin/divscan.py is executable.
    ```
    $ chmod +x divscan.py
    ```
3. Put divscan.py into your path of binaries, for example if a folder ``bin`` is already in your path, do
   ```
   $ cd <path-to-bin>
   $ ln -s <path-to-divscan.py> # use sudo depending on <path-to-bin>
   ```
4. Create conda environment:
   ```
   $ conda env create -f environment.yml
   ```

### Test and run

1. Activate environment
   ```
   $ conda activate divscan
   ```
2. Test
   ```
   $ shannon -h
   ```
