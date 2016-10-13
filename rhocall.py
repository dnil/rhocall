from cyvcf2 import *
import argparse
import numpy as np

version = "0.1"

### Argument parsing

parser = argparse.ArgumentParser(description='Call runs of autozygosity.')
parser.add_argument('--input_vcf','-i', type=str,
                    required=True, help='Input (sorted) vcf file')
parser.add_argument('--max_hets', '-m', type=float, default=2,
                    help='Max heterozygotes per Mb in a homozygous block')
parser.add_argument('--max_het_fraction', '-f', type=float, default=2,
                    help='Max heterozygotes over homozygotes fraction in a homozygous block')
# 4.8 might be a normal fraction

parser.add_argument('--minimum_homs', '-n', type=int, default = 5, 
                    help='Minimum absolute number of homozygotes to report a block')
parser.add_argument('--shortest_block', '-s', type=int, default=100000,
                    help='Shortest block')
parser.add_argument('--flag_UPD_at_fraction', '-u', type=float, default=0.3,
                    help='Flag UPD if homozygous blocks span this fraction of total chr size')
parser.add_argument('--individual','-k', type=int, default=0,
                    help='Index of individual in vcf/bcf, 0-based.')
parser.add_argument('--DEBUG', help='Enable debug output.', action='store_true')

args = parser.parse_args()

### Output format functions

def output_bed_header():
    print "#rhocall version %s" % (version)

def end_block(block_chr, block_start, block_end, block_hets, block_homs):
    print "%s\t%d\t%d\t%d\t%d\t%d" % ( block_chr, block_start, block_end, block_end-block_start+1, block_hets, block_homs )

def end_chr(current_chr, chr_size_spanned, chr_blocked, chr_hets, chr_homs):
    flag = "Normal"

    if chr_blocked / chr_size_spanned > args.flag_UPD_at_fraction:
        if current_chr == "X" or current_chr == "chrX":
            flag = "HemiX"
        elif current_chr == "Y" or current_chr == "chrY":
            flag = "HemiY"
        elif current_chr == "MT" or current_chr == "chrMT" or current_chr == "M" or current_chr == "chrM":
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
        elif current_chr == "MT" or current_chr == "chrMT" or current_chr == "M" or current_chr == "chrM":
            flag = "HetMT"
        
        print "%s\t1\t%d\tCHROMOSOME_TOTAL\t%d\t%d\t%d\t%s" % (current_chr, chr_size_spanned, chr_blocked, chr_hets, chr_homs, flag)

### Open vcf file 
# assuming a sorted vcf - check needed?
proband_vcf = VCF(args.input_vcf)

# start off output stream
output_bed_header()

# block parameters
hom_block_start = None
current_chr = None
hom_block_end = None
current_block_hets = 0
current_block_homs = 0

#chromosome parameters
chr_blocked = 0
chr_homs = 0 
chr_hets = 0
chr_last_pos_seen = None

# state variables for parser
seen_hom = False
block_ends = False
start_new_block = False
chr_ends = False

for var in proband_vcf:

    # ignore low qual and indels for now
    if var.FILTER:
        if args.DEBUG:
            print "Filter %s" % var.FILTER
        continue
    
    if var.end - var.start > 1:
        if args.DEBUG:
            print "Indel size %d" % (var.end - var.start)
        continue

    # new block if new chr
    if current_chr != None and var.CHROM != current_chr:
        block_ends = True
        chr_ends = True
        start_new_block = True
        if args.DEBUG:
            print "Leaving chromosome %s for %s." % ( current_chr, var.CHROM )
    elif seen_hom == False:
        if args.DEBUG:
            print type(var.gt_types)
            print len(var.gt_types)
            #        cpy = np.array(var.gt_bases)
            #        print "%d" % cpy.ndim
            print "%s %d" % (var.CHROM, var.start)

        if var.gt_types[args.individual] == proband_vcf.HOM_REF or var.gt_types[args.individual] == proband_vcf.HOM_ALT:
            # found homozygous variant
            start_new_block = True
            if args.DEBUG:
                print "Start new block."

    elif seen_hom == True:
        if var.gt_types[args.individual] == proband_vcf.HET:
            if var.gt_types[args.individual] == proband_vcf.HET:
                # update nr of hets. Tempting to leave the het count/block end at last hom?
                current_block_hets += 1
                hom_block_end = var.end
                if var.end > hom_block_start:
                    hom_block_end = var.end
                elif chr_last_pos_seen > hom_block_start:
                    hom_block_end = chr_last_pos_seen
                else:
                    print "ERROR - block size negative!"

                block_size = hom_block_end - hom_block_start + 1

                if args.DEBUG:
                    print "hets: %d homs: %d block_size: %d" % (current_block_hets, current_block_homs, block_size)
                if ( (current_block_hets) * 1000000 / block_size > args.max_hets ) or (current_block_hets) /current_block_homs > args.max_het_fraction :
                    block_ends = True
                    if args.DEBUG:
                        print "Leaving block with %d hets and block size %d." % ( current_block_hets, block_size )

        elif var.gt_types[args.individual] == proband_vcf.HOM_REF or var.gt_types[args.individual] == proband_vcf.HOM_ALT:
            current_block_homs += 1    
            current_block_end = var.end

    if block_ends:
        if hom_block_start != None:

            if var.end > hom_block_start:                
                hom_block_end = var.end
            elif chr_last_pos_seen > hom_block_start:
                hom_block_end = chr_last_pos_seen
            else:
                print "ERROR - block size negative!"
            
            # Nb - chr end/start!
            block_size = hom_block_end - hom_block_start + 1

            if current_block_hets * 1000000 / block_size > args.max_hets or current_block_hets / current_block_homs > args.max_het_fraction :
                if block_size > args.shortest_block and current_block_homs > args.minimum_homs:
                    end_block(current_chr, hom_block_start, hom_block_end, current_block_hets, current_block_homs)

            chr_blocked += block_size

            current_block_homs = 0
            current_block_hets = 0
            hom_block_start = None
            hom_block_end = None
            seen_hom = False

        else:
            if args.DEBUG:
                print "Block ends before one was started at %s %d." % (var.CHROM, var.start)

        block_ends = False
    
    if chr_ends:
        chr_size_spanned = chr_last_pos_seen
        end_chr(current_chr, chr_size_spanned, chr_blocked, chr_hets, chr_homs)
        
        chr_homs = 0
        chr_hets = 0
        current_chr = var.CHROM
        chr_blocked = 0
        chr_ends = False 

    if start_new_block:
        # new block; double check homoz for new chr.. just ignore if het.
        if var.gt_types[args.individual] == 0 or var.gt_types[args.individual] == 3:
            seen_hom = True
            current_block_homs = 1
            hom_block_start = var.start
            hom_block_end = var.end
            current_chr = var.CHROM

        start_new_block = False

    # finally, update chr totals and then return to beginning of loop to draw new var
    if var.gt_types[args.individual] == proband_vcf.HOM_REF or var.gt_types[args.individual] == proband_vcf.HOM_ALT:
        chr_homs += 1
    elif var.gt_types[args.individual] == proband_vcf.HET:
        chr_hets += 1

    chr_last_pos_seen = var.end

if (hom_block_end and hom_block_start):
# and check for block once more after seeing last var
    block_size = hom_block_end - hom_block_start + 1
    if current_block_hets / block_size * 1000000 > args.max_hets or current_block_hets /current_block_homs > args.max_het_fraction :
        if block_size > args.shortest_block and current_block_homs > args.minimum_homs:
            end_block(current_chr, hom_block_start, hom_block_end, current_block_hets, block_homs)



