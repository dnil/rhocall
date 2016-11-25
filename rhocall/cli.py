import click
import logging

from cyvcf2 import VCF

from rhocall.log import configure_stream, LEVELS
from .run_rho import run_rhocall
from .run_annotate import run_annotate
from .run_tally import run_tally

from .prints import (output_bed_header)

logger = logging.getLogger(__name__)

@click.group()
def cli():
    pass

@click.command()
@click.argument('vcf', type=click.Path(exists=True))
@click.option('--max_hets', '-m', 
    default=2.0,
    help='Max heterozygotes per Mb in a homozygous block'
)
@click.option('--max_het_fraction', '-f', 
    type=float, default=2,
    help='Max heterozygotes over homozygotes fraction in a homozygous block'
)
@click.option('--minimum_homs', '-n',
    default = 5, 
    help='Minimum absolute number of homozygotes to report a block'
)
@click.option('--shortest_block', '-s', 
    default=100000,
    help='Shortest block'
)
@click.option('--flag_upd_at_fraction', '-u',
    default=0.3,
    help='Flag UPD if homozygous blocks span this fraction of total chr size'
)
@click.option('--individual','-k', 
    type=int, 
    help='Index of individual in vcf/bcf, 0-based.'
)
@click.option('--block_constant', '-s', 
    default=100000,
    help='Should be adapted to type of analysis(exome or genome)'
)
@click.option('-v', '--verbose',
    count=True,
    default=2
)
@click.pass_context
def call(ctx,vcf, max_hets, max_het_fraction, minimum_homs, shortest_block, 
        flag_upd_at_fraction, individual, block_constant, verbose):
    """Call runs of autozygosity."""
    loglevel = LEVELS.get(min(verbose, 3))
    configure_stream(level=loglevel)
    
    proband_vcf = VCF(vcf)
    
    if len(proband_vcf.samples) > 1:
        try:
            individual = int(individual)
        except TypeError:
            logger.warning("Please specify which individual that should be checked")
            ctx.abort()
    else:
        individual = 0
    
    output_bed_header()
    
    run_rhocall(
        proband_vcf=proband_vcf, 
        block_constant=block_constant, 
        max_hets=max_hets, 
        max_het_fraction=max_het_fraction, 
        minimum_homs=minimum_homs, 
        shortest_block=shortest_block, 
        flag_UPD_at_fraction=flag_upd_at_fraction, 
        individual=individual
    )


@click.command()
@click.argument('roh', type=click.File('r'))
@click.option('--quality_threshold', '-q', 
    default=10.0,
    help='Minimum quality that counts towards region totals.'
)
@click.option('--flag_upd_at_fraction', '-u',    
    default=0.4,
    help='Flag UPD if this fraction of chr quality positions called AZ.'
)
@click.option('-v', '--verbose',
    count=True,
    default=2
)
@click.option('--output','-o', type=click.File('w'), default='-')
def tally(roh, quality_threshold, flag_upd_at_fraction, output, verbose):
    """Tally runs of autozygosity from rhofile.
    Accepts a bcftools roh style TSV-file with CHR,POS,AZ,QUAL."""
    loglevel = LEVELS.get(min(verbose, 3))
    configure_stream(level=loglevel)

    run_tally(
        roh=roh,
        quality_threshold=quality_threshold,
        flag_upd_at_fraction=flag_upd_at_fraction,
        output=output
    )


@click.command()
@click.argument('vcf', type=click.Path(exists=True))
@click.option('roh', '-r', type=click.File('r'),
              help='Bcftools roh style TSV file with CHR,POS,AZ,QUAL.')
@click.option('--quality_threshold', '-q', 
    default=10.0,
    help='Minimum quality calls that are imported in region totals.'
)
@click.option('--flag_upd_at_fraction', '-u',
    default=0.4,
    help='Flag UPD if this fraction of chr quality positions called AZ.'
)
@click.option('-v', '--verbose',
    count=True,
    default=2
)
@click.option('--output','-o',type=click.File('w'), default='-')
def annotate(vcf, roh, quality_threshold, flag_upd_at_fraction,output,verbose):
    """Markup VCF file using rho-calls. VCF files annotated with GENMOD style
    inheritance patterns are accepted."""
    loglevel = LEVELS.get(min(verbose, 3))
    configure_stream(level=loglevel)

    proband_vcf = VCF(vcf)
    
    # add this command to VCF header
    # add additional tag to VCF header INFO

    #output_vcfheader(output)

    run_annotate(
        proband_vcf=proband_vcf,
        roh=roh,        
        quality_threshold=quality_threshold,
        flag_upd_at_fraction=flag_upd_at_fraction,
        output=output
    )


cli.add_command(call)
cli.add_command(tally)
cli.add_command(annotate)
