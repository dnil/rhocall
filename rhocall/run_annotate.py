import logging

from rhocall import __version__ 

logger = logging.getLogger(__name__)

def run_annotate(proband_vcf, bed, quality_threshold, flag_upd_at_fraction, output):
    """Markup VCF file using rho-call BED file."""

    az_info_header={'ID' : 'AZ', 'Number' : 1, 'Type' : 'Flag', 
                    'Source' : 'rhocall', 'Version' : __version__,
                    'Description' : "Autozygous positon call"}
    proband_vcf.add_info_to_header(az_info_header);

    hw_info_header={'ID' : 'HW', 'Number' : 1, 'Type' : 'Flag', 
                    'Source' : 'rhocall', 'Version' : __version__,
                    'Description' : "Hardy-Weinberg equilibrium (non-autozyous) positon call"}
    proband_vcf.add_info_to_header(hw_info_header);

    # pyvcf2 does not seem to play with floats yet. Setting type to string for now.
    azqual_info_header={'ID' : 'AZQUAL', 'Number' : 1, 'Type' : 'String', 
                    'Source' : 'rhocall', 'Version' : __version__,
                    'Description' : 'Autozygous positon call quality'}
    proband_vcf.add_info_to_header(azqual_info_header);

    # pyvcf2 does not seem to play with floats yet. Setting type to string for now.
    azlength_info_header={'ID' : 'AZLENGTH', 'Number' : 1, 'Type' : 'String', 
                    'Source' : 'rhocall', 'Version' : __version__,
                    'Description' : 'Autozygous region length'}
    proband_vcf.add_info_to_header(azlength_info_header);

    aztype_info_header={'ID' : 'AZTYPE', 'Number' : 1, 'Type' : 'String', 
                    'Source' : 'rhocall', 'Version' : __version__,
                    'Description' : 'Autozygous region type'}
    proband_vcf.add_info_to_header(aztype_info_header);

    output.write(proband_vcf.raw_header)

    var = next(proband_vcf)

    for r in bed:

        if r[0] == '#':
            continue

        col = r.rstrip().split('\t')
        chrom = str(col[0])
        start = int(col[1])
        end = int(col[2])
        az = int(col[3])
        qual = float(col[4])
        # placeholder for future development: classify into UPD,SEX,DEL,IBD
        aztype = 'ND'
        azlength = end - start + 1

        logger.debug("Looking at bed window {0} {1}={2}".format(chrom,start,end))

#        print("looking for chrom %s %d" % (chrom, pos))
        passed_win = False
        found_var = False
        
        while not passed_win:
#            print("testing var chrom %s %d" % (var.CHROM, var.start))
            if var.CHROM == chrom and var.end >= start and var.end <= end:
                # variant in window - print, and go to next var
                if az == 1:
                    var.INFO['AZ'] = True
                else:
                    var.INFO['HW'] = True

                var.INFO['AZQUAL']=str(qual)
                var.INFO['AZLENGTH']=str(azlength)
                var.INFO['AZTYPE']=str(aztype)
                
                output.write(str(var))
                var = next(proband_vcf)
                found_var = True
            elif var.CHROM == chrom and var.start < start:
                # before next win (and on same chr) - write and pull new var
                output.write(str(var))
                var = next(proband_vcf)
            elif var.CHROM == chrom and var.end > end:
                # var is after last window position, same chr
                # pull new window, but no new variant and don't write var yet
                passed_win = True
            elif var.CHROM != chrom:
                # first time around, assume there are many variants, at least one
                # per bed interval
                
                # then we either just exited win by drawing var from new chr
                if found_var:
                    # pull next win, without deciding the fate of this var yet
                    passed_win = True
                else:
                    # or window is on new chr, and we need to draw new vars to
                    # get there
                    var = next(proband_vcf)
                    output.write(str(var))
            else:
                # not found, but passed the due position?!
                # var = next(proband_vcf)
                logger.error("Oops? Unexpected window/variant set!")



    
