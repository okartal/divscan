# shannon - Genome Analysis using Information Theory

## Requirements

* conda
* git

## Installation

1. Clone repository
2. Make shannon and test.sh executable
    ```shell
    $ chmod +x shannon
    $ chmod +x tests/test.sh
    ```
3. Put shannon in your path. Example:
   ```shell
   $ cd <path-to-bin>
   $ ln -s <path-to-shannon> # use sudo depending on <path-to-bin>
   ```
4. Create conda environment
   ```shell
   $ conda create --name <env-name> --file requirements.txt
   ```

## Testing

1. Activate environment
   ```shell
   $ source activate <env-name>
   ```
2. Run test
   ```shell
   $ cd tests
   $ ./test.sh
   ```
