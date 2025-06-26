import logging
from rhocall import __version__

logger = logging.getLogger(__name__)


def run_annotate_var(proband_vcf, roh, quality_threshold, flag_upd_at_fraction, output):
    """Markup VCF file using rho-calls. VCF files annotated with GENMOD style
    inheritance patterns are accepted."""

    az_info_header = {
        "ID": "AZ",
        "Number": 1,
        "Type": "Flag",
        "Source": "rhocall",
        "Version": __version__,
        "Description": "Autozygous positon call",
    }
    proband_vcf.add_info_to_header(az_info_header)

    hw_info_header = {
        "ID": "HW",
        "Number": 1,
        "Type": "Flag",
        "Source": "rhocall",
        "Version": __version__,
        "Description": "Hardy-Weinberg equilibrium (non-autozyous) positon call",
    }
    proband_vcf.add_info_to_header(hw_info_header)

    # pyvcf2 does not seem to play with floats yet. Setting type to string for now.
    azqual_info_header = {
        "ID": "AZQUAL",
        "Number": 1,
        "Type": "String",
        "Source": "rhocall",
        "Version": __version__,
        "Description": "Autozygous positon call quality",
    }
    proband_vcf.add_info_to_header(azqual_info_header)

    aztype_info_header = {
        "ID": "AZTYPE",
        "Number": 1,
        "Type": "String",
        "Source": "rhocall",
        "Version": __version__,
        "Description": "Autozygous region type",
    }
    proband_vcf.add_info_to_header(aztype_info_header)

    output.write(proband_vcf.raw_header)
    var = next(proband_vcf)

    for r in roh:

        if r[0] == "#":
            continue

        col = r.rstrip().split("\t")
        chrom = str(col[0])
        pos = int(col[1])
        az = int(col[2])
        qual = float(col[3])
        # placeholder for future development: classify into UPD,SEX,DEL,IBD
        aztype = "ND"

        #        print("looking for chrom %s %d" % (chrom, pos))
        found_var = False
        new_var = True

        while not found_var:
            #            print("testing var chrom %s %d" % (var.CHROM, var.start))
            if var.CHROM == chrom and var.end == pos:
                found_var = True
                if az == 1:
                    var.INFO["AZ"] = True
                else:
                    var.INFO["HW"] = True

                var.INFO["AZQUAL"] = str(qual)
                var.INFO["AZTYPE"] = str(aztype)

                output.write(str(var))
                new_var = True
            elif var.CHROM == chrom and var.start < pos:
                output.write(str(var))
                new_var = True
            else:
                # not found, but passed the due position?!
                # fake finding variant, and pull new tsv (but not new var)
                logger.warning(
                    "Variant given in rho-file was not found in variant file. "
                    "Are you sure this is the sample you want to tag?"
                )
                found_var = True
                new_var = False

            if new_var:
                var = next(proband_vcf)
