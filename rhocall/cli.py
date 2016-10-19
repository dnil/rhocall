import click
import logging

from cyvcf2 import VCF

from rhocall.log import configure_stream, LEVELS
from .run_rho import run_rhocall

from .prints import (output_bed_header)

logger = logging.getLogger(__name__)

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
def cli(ctx, vcf, max_hets, max_het_fraction, minimum_homs, shortest_block, 
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
