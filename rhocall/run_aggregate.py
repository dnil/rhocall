import logging

logger = logging.getLogger(__name__)

from .win import Win


def run_aggregate(roh, quality_threshold, output):
    """Aggregate bcftools TSV style per variant roh calls into windowed BED file.
    Trust calls at or above given quality threshold."""

    # current window object..
    win = Win()

    in_win = False

    for r in roh:

        if r[0] == "#":
            continue

        col = r.rstrip().split("\t")
        chr = str(col[0])
        pos = int(col[1])
        az = int(col[2])
        qual = float(col[3])

        if chr != win.chr:
            # new chr! leave last win and dump win (if printable)
            # dont start new win yet - await first 1..
            in_win = False
            win.dump(output)
            # clear the dumped window, setting chr to current and output false.
            win.reset(chr=chr, start=pos, end=pos, sum=0, count=0, print_me=False)

        if qual < quality_threshold:
            # low qual cannot start win, cannot end win - only extend
            if in_win and az == 1:
                # init new tentative window instead..
                win.extend(pos, qual)
            continue

        if in_win and az == 1:
            # extend win
            win.extend(pos, qual)

        elif in_win and az == 0:
            # dump win to output, drop.
            win.dump(output)
            win.reset(chr=chr, start=pos, end=pos, sum=0, count=0, print_me=False)
            in_win = False

        elif not in_win and az == 1:
            in_win = True
            win.reset(chr=chr, start=pos, end=pos, sum=qual, count=1, print_me=True)
        # if not in win and not az - do nothing..

    # just in case the last part of the last chr is in an az block
    if in_win:
        win.dump(output)
