# rhocall
Call regions of homozygosity and make tentative UPD calls

## Usage ##

```
Usage: rhocall [OPTIONS] VCF

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

Usage: rhocall [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  annotate  Markup VCF file using rho-calls.
  call      Call runs of autozygosity.
  tally     Tally runs of autozygosity from rhofile.

Usage: rhocall call [OPTIONS] VCF

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

Usage: rhocall tally [OPTIONS] ROH

  Tally runs of autozygosity from rhofile. Accepts a bcftools roh style TSV-
  file with CHR,POS,AZ,QUAL.

Options:
  -q, --quality_threshold FLOAT   Minimum quality that counts towards region
                                  totals.
  -u, --flag_upd_at_fraction FLOAT
                                  Flag UPD if this fraction of chr quality
                                  positions called AZ.
  -v, --verbose
  -o, --output FILENAME
  --help                          Show this message and exit.

Usage: rhocall annotate [OPTIONS] VCF

  Markup VCF file using rho-calls. VCF files annotated with GENMOD style
  inheritance patterns are accepted.

Options:
  -r FILENAME                     Bcftools roh style TSV file with
                                  CHR,POS,AZ,QUAL.
  -q, --quality_threshold FLOAT   Minimum quality calls that are imported in
                                  region totals.
  -u, --flag_upd_at_fraction FLOAT
                                  Flag UPD if this fraction of chr quality
                                  positions called AZ.
  -v, --verbose
  -o, --output FILENAME
  --help                          Show this message and exit.

```

## Example ##

```
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' anon-SweGen_STR_NSPHS_1000samples_snp_freq_hg19.vcf.gz | bgzip -c > anon_SweGen_161019_snp_freq_hg19.tab.gz
bcftools roh --AF-file anon_SweGen_161019_snp_freq_hg19.tab.gz -I 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.bcf > 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.roh
rhocall tally 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.roh
rhocall annotate -r 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.roh 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.bcf
```

## Installation ##

```
pip install -e .
```


