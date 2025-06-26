from __future__ import division

from rhocall import __version__


def output_bed_header():
    print("#rhocall version {0}".format(__version__))


def end_block(block_chr, block_start, block_end, block_hets, block_homs):
    print(
        "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(
            block_chr, block_start, block_end, block_end - block_start + 1, block_hets, block_homs
        )
    )


def end_chr(current_chr, chr_size_spanned, chr_blocked, chr_hets, chr_homs, flag_UPD_at_fraction):
    flag = "Normal"

    if (chr_blocked / chr_size_spanned) > flag_UPD_at_fraction:
        if current_chr == "X" or current_chr == "chrX":
            flag = "HemiX"
        elif current_chr == "Y" or current_chr == "chrY":
            flag = "HemiY"
        elif (
            current_chr == "MT"
            or current_chr == "chrMT"
            or current_chr == "M"
            or current_chr == "chrM"
        ):
            flag = "HemiMT"
        else:
            flag = "UPD"
    else:
        if current_chr == "X" or current_chr == "chrX":
            flag = "HetX"
        elif current_chr == "Y" or current_chr == "chrY":
            if chr_hets + chr_homs < 1000:
                flag = "NoY"

            flag = "HetY"
        elif (
            current_chr == "MT"
            or current_chr == "chrMT"
            or current_chr == "M"
            or current_chr == "chrM"
        ):
            flag = "HetMT"

        print(
            "{0}\t1\t{1}\tCHROMOSOME_TOTAL\t{2}\t{3}\t{4}\t{5}".format(
                current_chr, chr_size_spanned, chr_blocked, chr_hets, chr_homs, flag
            )
        )
