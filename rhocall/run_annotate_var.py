import logging

logger = logging.getLogger(__name__)

def run_annotate_var(proband_vcf, roh, quality_threshold, flag_upd_at_fraction, output):
    """Markup VCF file using rho-calls. VCF files annotated with GENMOD style
    inheritance patterns are accepted."""

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

    for r in roh:
        
        if r[0] == '#':
            continue

        col = r.rstrip().split('\t')
        chr = str(col[0])
        pos = int(col[1])
        az = int(col[2])
        qual = float(col[3])

#        print("looking for chr %s %d" % (chr, pos))
        found_var = False
        while not found_var:
#            print("testing var chr %s %d" % (var.CHROM, var.start))
            if var.CHROM == chr and var.end == pos:
                found_var = True
                if az == 1:
                    var.INFO['AZ'] = True
                else:
                    var.INFO['HW'] = True
                    
                var.INFO['AZQUAL']=str(qual)
                output.write(str(var))
            elif var.CHROM == chr and var.start < pos:   
                output.write(str(var))

            else:
                # not found, but passed the due position?!
                pass

            var = next(proband_vcf)
    
