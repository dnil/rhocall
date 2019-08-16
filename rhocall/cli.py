import click
import logging
import inspect
import os

from cyvcf2 import VCF

from rhocall.log import configure_stream, LEVELS
from .run_rho import run_rhocall
from .run_annotate import run_annotate
from .run_annotate_bcfroh import run_annotate_rg
from .run_annotate_var import run_annotate_var
from .run_tally import run_tally
from .run_aggregate import run_aggregate
from .run_viz import generate_bins, generate_plots, generate_wig, extract_roh

from .prints import (output_bed_header)

from rhocall import __version__

logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
        '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

@click.group()
@click.version_option(__version__)
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
def call(ctx,vcf, max_hets, max_het_fraction, minimum_homs, shortest_block,
        flag_upd_at_fraction, individual, block_constant, verbose):
    """Call runs of autozygosity. (deprecated: use bcftools roh instead."""
    loglevel = LEVELS.get(min(verbose, 3))
    configure_stream(level=loglevel)

    proband_vcf = VCF(vcf)

    if len(proband_vcf.samples) > 1:
        try:
            individual = int(individual)
        except TypeError:
            logger.warning("Please specify which individual to check.")
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

    logger.info("Running rhocall tally {0}".format(__version__))

    run_tally(
        roh=roh,
        quality_threshold=quality_threshold,
        flag_upd_at_fraction=flag_upd_at_fraction,
        output=output
    )


@click.command()
@click.argument('roh', type=click.File('r'))
@click.option('--quality_threshold', '-q',
    default=10.0,
    help='Minimum quality trusted to start or end ROH-windows.'
)
@click.option('-v', '--verbose',
    count=True,
    default=2
)
@click.option('--output','-o', type=click.File('w'), default='-')
def aggregate(roh, quality_threshold, output, verbose):
    """Aggregate runs of autozygosity from rhofile into windowed rho BED file.
    Accepts a bcftools roh style TSV-file with CHR,POS,AZ,QUAL."""
    loglevel = LEVELS.get(min(verbose, 3))
    configure_stream(level=loglevel)

    logger.info("Running rhocall aggregate {0}".format(__version__))

    run_aggregate(
        roh=roh,
        quality_threshold=quality_threshold,
        output=output
    )


@click.command()
@click.argument('vcf', type=click.Path(exists=True))
@click.option('roh', '-r', type=click.File('r'),
              help='Bcftools roh style TSV file with CHR,POS,AZ,QUAL.')
@click.option('--v14/--no-v14', default=True,
	      help='Bcftools v1.4 or newer roh file including RG calls.')
@click.option('bed', '-b', type=click.File('r'),
              help='BED file with AZ windows.')
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
def annotate(vcf, roh, bed, v14, quality_threshold, flag_upd_at_fraction,
	     output, verbose):
    """Markup VCF file using rho-calls. Use BED file to mark all variants in AZ
    windows. Alternatively, use a bcftools v>=1.4 file with RG entries to mark
    all vars. With the --no-v14 flag, use an older bcftools v<=1.2 style roh TSV
    to mark only selected AZ variants. Roh is broken in bcftools v1.3
    - do not use."""
    loglevel = LEVELS.get(min(verbose, 3))
    configure_stream(level=loglevel)

    proband_vcf = VCF(vcf)

    # add this command to VCF header

    ## This is for logging the command line string ##
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    argument_list = [
        i+'='+str(values[i]) for i in values if values[i] and
        i not in ['frame']
    ]

    logger.info("Running rhocall annotate  {0}".format(__version__))
    logger.debug("Arguments: {0}".format(', '.join(argument_list)))

    ## add additional tags to VCF header
    proband_vcf.add_to_header('##rhocall_version={0}'.format(__version__))
    proband_vcf.add_to_header("##rhocall_arguments={0}".format(', '.join(argument_list)))

    if roh and not bed and not v14:
         run_annotate_var(
		proband_vcf=proband_vcf,
		roh=roh,
		quality_threshold=quality_threshold,
		flag_upd_at_fraction=flag_upd_at_fraction,
		output=output
		)
    elif roh and v14 and not bed:
        run_annotate_rg(
		proband_vcf=proband_vcf,
		bcfroh=roh,
		quality_threshold=quality_threshold,
		flag_upd_at_fraction=flag_upd_at_fraction,
		output=output
		)
    elif bed and not roh:
        run_annotate(
		proband_vcf=proband_vcf,
		bed=bed,
		quality_threshold=quality_threshold,
		flag_upd_at_fraction=flag_upd_at_fraction,
		output=output
		)
    else:
        click.echo("""Cannot use both BED and ROH at once. Please apply
                    them sequentially instead.""")


@click.command()
@click.argument('vcf', type=click.File('r'), required=True) # help='Input (sorted) vcf file'
@click.option('--out_dir', type=click.Path('w'), exists=False,required=True,
                help='Output directory. The files will be named out_dir/chr.png. One picture is drawn per chromosome.')
@click.option('--wig/--no-wig', default=True,
                help='Produce wig file.')
@click.option('--pointsize','-p', type=int,default=8,
                help='Size of the points (pixels)')
@click.option('--rho', '-r', type=click.File('r'), required=True,
                help='Input RHO file produced from rhocall')
@click.option('--minsnv', '-m',type=int,default=2,
                help='Minimum number of snvs for each plotted bin')
@click.option('--maxsnv', '-M', type=int, default=20,
                help='Maximum number of snvs for each plotted bin')
@click.option('--minaf', type=float, default=0.1,
                help='Minimum allele frequency. This variable must be set to 0 if the allele frequency is not annotated.')
@click.option('--maxaf', type=float, default=0.9,
                help='Maximum allele frequency')
@click.option('--aftag', type=str, default="1000GAF",
                help='The allele frequency tag to use.')
@click.option('--minqual','-q', type=int, default=98,
                help='Do not add SNVs to a bin if their quality is less than this value.')
@click.option('--mnv/--no-mnv', default=False, help='Include MNV')
@click.option('--window','-w', type=int, default=10000,
                help='Window size(bases)')
@click.option('--rsid/--no-rsid','-s', default=False,
                help='Skip variants not containing an rsid')
@click.option('--filter/--no-filter','-n', default=True,
                help='include variants, even if they are not labeled PASS')
@click.option('-v', '--verbose', count=True, default=2)
def viz(vcf, out_dir, wig, pointsize, rho, minsnv, maxsnv, minaf, maxaf, aftag,
        minqual, mnv, window, rsid, filter, verbose):
    """Plot binned zygosity and RHO-regions."""

    loglevel = LEVELS.get(min(verbose, 3))
    configure_stream(level=loglevel)

    os.system( "mkdir {}".format(out_dir) )
    binned_zygosity=generate_bins(vcf, window, filter, mnv, minqual,
                                    rsid, minaf, aftag, maxaf, minsnv)
    roh=extract_roh(rho)
    if wig:
        out_file_basename = out_dir + "/output"
        generate_wig(binned_zygosity,roh, window, out_file_basename)
    else:
        generate_plots(binned_zygosity,roh,window, pointsize, out_dir)

cli.add_command(call)
cli.add_command(tally)
cli.add_command(aggregate)
cli.add_command(annotate)
cli.add_command(viz)
