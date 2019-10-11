# Shannon - Measuring Genomic Diversity using Information Theory

*NOTE*: The command line app is usable but installation, API, and implementation may change in the future

## Requirements

* git (see <https://git-scm.com/>)
* Miniconda for Python3 (see <https://conda.io/miniconda.html>)

## Install

1. Clone the repository.
2. Go into the folder and ensure that shannon is executable.
    ```
    $ chmod +x shannon
    ```
3. Put shannon in your path. Example:
   ```
   $ cd <path-to-bin>
   $ ln -s <path-to-shannon> # use sudo depending on <path-to-bin>
   ```
4. Create conda environment:
   ```
   $ conda env create -f environment.yml
   ```

## Run

1. Activate environment
   ```
   $ conda activate shannon
   ```
2. Run help
   ```
   $ shannon -h
   ```
