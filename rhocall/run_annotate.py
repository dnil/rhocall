import logging

logger = logging.getLogger(__name__)

def run_annotate(proband_vcf, bed, quality_threshold, flag_upd_at_fraction, output):
    """Markup VCF file using rho-call BED file."""

    az_info_header={'ID' : 'AZ', 'Number' : 1, 'Type' : 'Flag', 
                    'Source' : 'rhocall', 'Version' : '0.1',
                    'Description' : "Autozygous positon call"}
    proband_vcf.add_info_to_header(az_info_header);

    hw_info_header={'ID' : 'HW', 'Number' : 1, 'Type' : 'Flag', 
                    'Source' : 'rhocall', 'Version' : '0.1',
                    'Description' : "Hardy-Weinberg equilibrium (non-autozyous) positon call"}
    proband_vcf.add_info_to_header(hw_info_header);

    # pyvcf2 does not seem to play with floats yet. Setting type to string for now.
    azqual_info_header={'ID' : 'AZQUAL', 'Number' : 1, 'Type' : 'String', 
                    'Source' : 'rhocall', 'Version' : '0.1',
                    'Description' : 'Autozygous positon call quality'}
    proband_vcf.add_info_to_header(azqual_info_header);

    var = next(proband_vcf)

    for r in bed:

        if r[0] == '#':
            continue

        col = r.rstrip().split('\t')
        chr = str(col[0])
        start = int(col[1])
        end = int(col[2])
        az = int(col[3])
        qual = float(col[4])

#        print("looking for chr %s %d" % (chr, pos))
        passed_win = False
        while not passed_win:
#            print("testing var chr %s %d" % (var.CHROM, var.start))
            if var.CHROM == chr and var.end >= start and var.end <= end:

                if az == 1:
                    var.INFO['AZ'] = True
                else:
                    var.INFO['HW'] = True
                    
                var.INFO['AZQUAL']=str(qual)                
                output.write(str(var))
                var = next(proband_vcf)
            elif var.CHROM == chr and var.start < start:
                # before next win (and not in win) - write and pull new var
                output.write(str(var))
                var = next(proband_vcf)
            elif var.CHROM == chr and var.end > end:
                # var is after last window position, same chr
                # pull new window, but no new variant and don't write var yet
                passed_win = True
            elif var.CHROM != chr:
                # we ask that chromosomes are sorted the same way in both files
                if var.end > start:
                    # likely, we have pulled a new win on a new chr. 
                    # dump var without label, and pull new
                    output.write(str(var))
                    var = next(proband_vcf)
                elif var.end <= start:
                    # if we have pulled a new var on a new chr 
                    # say, skipping one/first chr entirely, or last win extended
                    # to chromosome end
                    # pull new window, but no new variant and don't write var yet
                    passed_win = True
            else:
                # not found, but passed the due position?!
                # var = next(proband_vcf)
                click.echo("Oops? Dunno what to do with this win and var sequence!")



    
