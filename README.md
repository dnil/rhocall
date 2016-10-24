# rhocall
Call regions of homozygosity and make tentative UPD calls

## Usage ##

```Usage: rhocall [OPTIONS] VCF

  Call runs of autozygosity.

Options:
  -m, --max_hets FLOAT            Max heterozygotes per Mb in a homozygous
                                  block
  -f, --max_het_fraction FLOAT    Max heterozygotes over homozygotes fraction
                                  in a homozygous block
  -n, --minimum_homs INTEGER      Minimum absolute number of homozygotes to
                                  report a block
  -s, --shortest_block INTEGER    Shortest block
  -u, --flag_upd_at_fraction FLOAT
                                  Flag UPD if homozygous blocks span this
                                  fraction of total chr size
  -k, --individual INTEGER        Index of individual in vcf/bcf, 0-based.
  -s, --block_constant INTEGER    Should be adapted to type of analysis(exome
                                  or genome)
  -v, --verbose
  --help                          Show this message and exit.
```


## Installation ##

```pip install -e .```


