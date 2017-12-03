# rhocall
Call regions of homozygosity and make tentative UPD calls

## Usage ##

```
Usage: rhocall [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  aggregate  Aggregate runs of autozygosity from rhofile...
  annotate  Markup VCF file using rho-calls.
  call      Call runs of autozygosity.
  tally     Tally runs of autozygosity from rhofile.
```
### rhocall call ###
```
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
```

### rhocall aggregate ####
```
Usage: rhocall aggregate [OPTIONS] ROH

  Aggregate runs of autozygosity from rhofile into windowed rho BED file.
  Accepts a bcftools roh style TSV-file with CHR,POS,AZ,QUAL.

Options:
  -q, --quality_threshold FLOAT  Minimum quality trusted to start or end ROH-
                                 windows.
  -v, --verbose
  -o, --output FILENAME
  --help                         Show this message and exit.
```

### rhocall tally ###
```
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
```

### rhocall annotate ###
```
Usage: rhocall annotate [OPTIONS] VCF

  Markup VCF file using rho-calls. Use BED file to mark all variants in AZ
  windows. Alternatively, use a bcftools v>=1.4 file with RG entries to mark
  all vars. With the --no-v14 flag, use an older bcftools v<=1.2 style roh
  TSV to mark only selected AZ variants. Roh is broken in bcftools v1.3  -
  do not use.

Options:
  -r FILENAME                     Bcftools roh style TSV file with
                                  CHR,POS,AZ,QUAL.
  --v14 / --no-v14                Bcftools v1.4 or newer roh file including RG
                                  calls.
  -b FILENAME                     BED file with AZ windows.
  -q, --quality_threshold FLOAT   Minimum quality calls that are imported in
                                  region totals.
  -u, --flag_upd_at_fraction FLOAT
                                  Flag UPD if this fraction of chr quality
                                  positions called AZ.
  -v, --verbose
  -o, --output FILENAME
  --help                          Show this message and exit.

```

### rhoviz ###

```
Usage: rhoviz [OPTIONS]

  Visualise the rhocall output file. Genomic regions labeled rho are visualised as red lines. Additionally, the fraction of homozygous snps are visualised as black dots.


Options:
  -r FILENAME                     Input RHO file produced from rhocall

  --help                          Show this message and exit.

  -i                              Input vcf file

  -d                              output directory, the files will be named dir/chr.png,

  -w                              window size(bases)

  -m                              minimum number of snvs for each plotted bin

  -M                              maximum number of snvs for each plotted bin

  --minaf                         minimum allele frequency(this variable must be set to 0 if the allele frequency is not annotated)

  --maxaf                         maximum allele frequency

  --aftag AFTAG                   the allele frequency tag

  -q                              do not add snvs to a bin if there quality is less than this value
  -p                              Size of the points (pixels)

  -n                              include variants, even if they are not labeled PASS


```

## Examples ##

### Suggested workflow ###

#### Preparation ####
```
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' popfreq.vcf.gz | bgzip -c > popfreq.tab.gz
```

#### Call ROH with bcftools ####
Please see the [samtools project](https://samtools.github.io/bcftools/) for installation instructions, and 
please refer to [Narasimhan et al, 2016](http://bioinformatics.oxfordjournals.org/content/early/2016/01/30/bioinformatics.btw044) regarding method details.

```
bcftools roh --AF-file popfreq.tab.gz -I sample.bcf > sample.roh
```

#### Aggregate ROH calls into windows, and mark up variant file (VCF/BCF) ####
```
rhocall aggregate sample.roh -o sample.roh.bed
rhocall annotate -b sample.roh.bed -o sample.rho.vcf sample.bcf
```

#### Obtain per chromosome overview ####
```
rhocall tally sample.roh -o sample.roh.tally.tsv
```

### Additional usage examples ###

```
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' anon-SweGen_STR_NSPHS_1000samples_snp_freq_hg19.vcf.gz | bgzip -c > anon_SweGen_161019_snp_freq_hg19.tab.gz
bcftools roh --AF-file anon_SweGen_161019_snp_freq_hg19.tab.gz -I 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.bcf > 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.roh
# bcftools <=1.2
rhocall tally 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.roh -o 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.roh.tally.tsv
rhocall annotate --no-v14 -r 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.roh 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.bcf -o 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.roh.vcf
rhocall aggregate 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.roh -o 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.roh.bed
rhocall annotate -b 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.roh.bed -o 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.rho.vcf 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.bcf
# bcftools >=1.4
rhocall annotate --v14 -r 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.14.roh -o 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.14.rho.vcf 2016-14676_sorted_md_rreal_brecal_gvcf_vrecal_comb_BOTH.bcf
```

## Test files ##
The test directory contains test files from the [BCFtools/RoH project](https://samtools.github.io/bcftools/howtos/roh-calling.html).

## Installation ##
The cyvcf2 install process appears to be jinxed on certain systems/setups. 
In practice this means that a chained pip install on a naive system may fail. Installation of each requirement for cyvcf2 prior to installing it appears to work unconditionally.
```
pip install numpy; pip install Cython
pip install -r requirements.txt
pip install -e .
```


