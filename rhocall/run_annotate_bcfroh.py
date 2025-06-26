import logging

from rhocall import __version__

logger = logging.getLogger(__name__)


def run_annotate_rg(proband_vcf, bcfroh, quality_threshold, flag_upd_at_fraction, output):
    """Markup VCF file using rho-call BED file."""

    az_info_header = {
        "ID": "AZ",
        "Number": 0,
        "Type": "Flag",
        "Source": "rhocall",
        "Version": __version__,
        "Description": "Autozygous positon call",
    }
    proband_vcf.add_info_to_header(az_info_header)

    hw_info_header = {
        "ID": "HW",
        "Number": 0,
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

    # pyvcf2 does not seem to play with floats yet. Setting type to string for now.
    azlength_info_header = {
        "ID": "AZLENGTH",
        "Number": 1,
        "Type": "String",
        "Source": "rhocall",
        "Version": __version__,
        "Description": "Autozygous region length",
    }
    proband_vcf.add_info_to_header(azlength_info_header)

    # pyvcf2 does not seem to play with floats yet. Setting type to string for now.
    azlength_info_header = {
        "ID": "AZMARKERS",
        "Number": 1,
        "Type": "String",
        "Source": "rhocall",
        "Version": __version__,
        "Description": "Autozygous region length",
    }
    proband_vcf.add_info_to_header(azlength_info_header)

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

    var = next(proband_vcf, False)

    found_any_var = False

    for r in bcfroh:

        if r[0] == "#":
            continue

        col = r.rstrip().split("\t")
        if col[0] != "RG":
            continue

        # trust RG entries
        if col[0] == "RG":
            sample = str(col[1])
            chrom = str(col[2])
            start = int(col[3])
            end = int(col[4])
            azlength = int(col[5])
            nmarkers = int(col[6])
            qual = float(col[7])

        found_var = False

        # placeholder for future development: classify into UPD,SEX,DEL,IBD
        aztype = "ND"
        azlength = end - start + 1

        logger.debug(
            "looking for roh window %s %d-%d nmarkers %d qual %f"
            % (chrom, start, end, nmarkers, qual)
        )

        passed_win = False

        while var and not passed_win:
            #            print("testing var chrom %s %d" % (var.CHROM, var.start))
            if var.CHROM == chrom and var.end >= start and var.end <= end:
                # variant in window - print, and go to next var
                var.INFO["AZ"] = True
                var.INFO["AZQUAL"] = str(qual)
                var.INFO["AZLENGTH"] = str(azlength)
                var.INFO["AZMARKERS"] = str(nmarkers)
                var.INFO["AZTYPE"] = str(aztype)
                found_var = True
                found_any_var = True

                output.write(str(var))
                var = next(proband_vcf, False)

            elif var.CHROM == chrom and var.start < start:
                # before next win (and on same chr) - write and pull new var
                output.write(str(var))
                var = next(proband_vcf, False)

            elif var.CHROM == chrom and var.end > end:
                # var is after last window position, same chr
                # pull new window, but no new variant and don't write var yet
                passed_win = True
                # when written?
            elif var.CHROM != chrom:
                # first time around, assume there are many variants, at least one
                # per bed interval.

                if found_var:
                    # then we either just exited win by drawing var from new chr
                    # so pull next win, without deciding the fate of this var yet
                    passed_win = True
                else:
                    # or window is on new chr, and we need to draw new vars to
                    # get there - essentially "before next win"
                    output.write(str(var))
                    var = next(proband_vcf, False)
                    logger.debug(
                        "Win chr %s not same as var chr %s: keep drawing new vars (end %s)."
                        % (chrom, var.CHROM, var.end)
                    )

            else:
                # not found, but passed the due position?!
                # var = next(proband_vcf)
                logger.error("Oops? Unexpected window/variant set!")

    # out of windows. print any remaining unflagged vars.
    while var:
        output.write(str(var))
        var = next(proband_vcf, False)

    if not found_any_var:
        logger.warning(
            "No variants found for at least one window."
            " - were correctly matched files used for annotation?"
        )
