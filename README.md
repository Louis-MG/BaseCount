# BaseCount

Counts bases proportions in sequences. Parallel execution.

# Installation

```bash
mamba create -n basecount
mamba activate basecount
mamba install pip
pip install argparse biopython multiprocess
```

# Usage

```bash
python3 BaseCount.py -f $(test/*.gz) -t 3 -o output
```

Works with (un)compressed fastq. Parallelised.
