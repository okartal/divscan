# Shannon - Measuring Genomic Diversity using Information Theory

## Requirements

* git (see <https://git-scm.com/>)
* Miniconda for Python3 (see <https://conda.io/miniconda.html>)

## Install

1. Clone the repository.
2. Make shannon and test.sh executable.
    ```
    $ chmod +x shannon
    $ chmod +x tests/test.sh
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
   $ source activate shannon
   ```
2. Run help
   ```
   $ shannon -h
   ```
