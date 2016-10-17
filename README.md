# shannon - Genome Analysis using Information Theory

## Testing

1. Create conda environment and activate it.

   ```shell
   $ conda create --name <env-name> --file requirements.txt
   ```

2. Make sure shannon and tests/test.sh are executable

   ```shell
   $ chmod u+x shannon
   ```

3. Make sure shannon is in your path. My preferred approach is to create a symlink in /usr/local/bin to the shannon script in the repository.

   ```shell
   $ cd <path-to-bin>
   $ ln -s <path-to-shannon> # use sudo depending on <path-to-bin>
   ```

4. Test

   ```shell
   $ cd tests
   $ ./test.sh
   ```
