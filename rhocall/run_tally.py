import logging

logger = logging.getLogger(__name__)

def run_tally(roh, quality_threshold, flag_upd_at_fraction, output):
    """Markup VCF file using rho-calls. VCF files annotated with GENMOD style
    inheritance patterns are accepted."""

    chr_hw = dict()
    chr_az = dict()
    # output chr in order of appearance
    chr_ordered_list = []

    for r in roh:

        if r[0] == '#':
            continue

        col = r.rstrip().split('\t')
        chr = str(col[0])
        pos = int(col[1])
        az = int(col[2])
        qual = float(col[3])

        if chr not in chr_ordered_list:
            chr_ordered_list = chr_ordered_list + [chr]
            if chr not in chr_az:
                chr_az[chr] = 0
            if chr not in chr_hw:
                chr_hw[chr] = 0

        if qual >= quality_threshold:                      
            if az == 0:
                if chr in chr_hw:
                    chr_hw[chr] = chr_hw[chr] + 1
                else:
                    chr_hw[chr] = 1
            else:
                if chr in chr_az:
                    chr_az[chr] = chr_az[chr] + 1
                else:
                    chr_az[chr] = 1

        else:
            # ignore low quality
            pass
    
    output.write("chr\thw\taz\taz\tflag\n")
    for chr in chr_ordered_list:
        az_fraction = chr_az[chr]/(chr_hw[chr]+chr_az[chr])
        flag = ""
        if az_fraction >= flag_upd_at_fraction:
            flag = "HIGH_AZ"
            # alert only for chrX, chrY calls?
        output.write("%s\t%d\t%d\t%.2f\t%s\n" % (chr, chr_hw[chr], chr_az[chr], az_fraction, flag))
