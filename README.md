# sierra-local
sierra-local is a local implementation of Stanford University's HIVdb web interface Sierra. It generates drug resistance predictions given HIV-1 sequence data, and outputs this in a standard format.

## Rationale

## Dependencies
- Python 2.7+
- [NucAmino](https://github.com/hivdb/nucamino)

## Installation
1. Clone this repository.
    ```
    git clone https://github.com/PoonLab/sierra-local.git
    ```
2. [Build `NucAmino` binary using Docker](https://github.com/hivdb/nucamino). The NucAmino binary needs to be located in the `sierra-local` root directory.
3. Pull algorithm and essential data using `update_HIVDB.py`. This needs to be done before first use.
    ```cmd
    python ./scripts/update_HIVDB.py
    ```

## Using sierra-local
Example use case:
```
python sierralocal.py SEQUENCES.fasta -o OUTPUT.json
```

## About Us
This project was developed at the Poon Lab under the Department of Pathology and Laboratory Medicine, Schulich School of Medicine and Dentistry, Western University, London, Ontario.
