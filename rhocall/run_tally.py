import logging

logger = logging.getLogger(__name__)

def run_tally(roh, quality_threshold, flag_upd_at_fraction, output):
    """Markup VCF file using rho-calls. VCF files annotated with GENMOD style
    inheritance patterns are accepted."""

    chrom_hw = dict()
    chrom_az = dict()
    # output chrom in order of appearance
    chrom_ordered_list = []

    for r in roh:

        if r[0] == '#':
            continue

        col = r.rstrip().split('\t')
        chrom = str(col[0])
        pos = int(col[1])
        az = int(col[2])
        qual = float(col[3])

        if chrom not in chrom_ordered_list:
            chrom_ordered_list = chrom_ordered_list + [chrom]
            if chrom not in chrom_az:
                chrom_az[chrom] = 0
            if chrom not in chrom_hw:
                chrom_hw[chrom] = 0

        if qual >= quality_threshold:               
            if az == 0:
                if chrom in chrom_hw:
                    chrom_hw[chrom] = chrom_hw[chrom] + 1
                else:
                    chrom_hw[chrom] = 1
            else:
                if chrom in chrom_az:
                    chrom_az[chrom] = chrom_az[chrom] + 1
                else:
                    chrom_az[chrom] = 1

        else:
            # ignore low quality
            pass
    
    output.write("chrom\thw\taz\tazf\tflag\n")
    for chrom in chrom_ordered_list:
        az_fraction = chrom_az[chrom] / (chrom_hw[chrom] + chrom_az[chrom])
        flag = ""
        if az_fraction >= flag_upd_at_fraction:
            flag = "HIGH_AZ"
            # alert only for chromX, chromY calls?
        output.write("%s\t%d\t%d\t%.2f\t%s\n" % (chrom, chrom_hw[chrom], chrom_az[chrom], az_fraction, flag))
