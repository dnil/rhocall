import logging

from .prints import (end_block, end_chr)
from .utils import (skip_variant, check_homozygote, check_heterozygote,
                    check_block_end, check_valid_block)

logger = logging.getLogger(__name__)

def run_rhocall(ctx,proband_vcf, block_constant, max_hets, max_het_fraction, 
                 minimum_homs, shortest_block, flag_UPD_at_fraction, individual):
    """docstring for run_rho_call"""
    block_constant = 1000000
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
        
        if not skip_variant(var):
            # new block if new chr
            if (current_chr != None and var.CHROM != current_chr):
                block_ends = True
                chr_ends = True
                start_new_block = True
                logger.debug("Leaving chromosome %s for %s." % 
                              (current_chr, var.CHROM))

            elif seen_hom == False:
                logger.debug(type(var.gt_types))
                logger.debug(len(var.gt_types))
                logger.debug("%s %d" % (var.CHROM, var.start))

                if check_homozygote(var, proband_vcf, individual):
                    # found homozygous variant
                    start_new_block = True
                    logger.debug("Start new block.")
            
            elif seen_hom == True:
                if check_heterozygote(var, proband_vcf, individual):
                    # update nr of hets. Tempting to leave the het count/block end at last hom?
                    current_block_hets += 1
                    hom_block_end = var.end
                    if var.end > hom_block_start:
                        hom_block_end = var.end
                    elif chr_last_pos_seen > hom_block_start:
                        hom_block_end = chr_last_pos_seen
                    else:
                        logger.error("Block size negative!")
                        #Raise exception??

                    block_size = hom_block_end - hom_block_start + 1

                    logger.debug("hets: %d homs: %d block_size: %d" % 
                                 (current_block_hets, current_block_homs, block_size))
                    
                    block_ends = check_block_end(
                        nr_hets = current_block_hets,
                        nr_homs = current_block_homs,
                        constant = block_constant,
                        block_size = block_size,
                        max_hets = max_hets,
                        max_het_fraction = max_het_fraction
                    )
                    if block_ends:
                        logger.debug("Leaving block with %d hets and block"\
                                     " size %d." % ( current_block_hets, block_size ))
            
                elif check_homozygote(var, proband_vcf, individual):
                    current_block_homs += 1    
                    current_block_end = var.end
            
            if block_ends:
                if hom_block_start != None:

                    if var.end > hom_block_start:
                        hom_block_end = var.end
                    elif chr_last_pos_seen > hom_block_start:
                        hom_block_end = chr_last_pos_seen
                    else:
                        logger.error("Block size negative!")
                        # Raise exception??
                
                    # Nb - chr end/start!
                    block_size = hom_block_end - hom_block_start + 1
                    
                    block_ends = check_block_end(
                        nr_hets = current_block_hets,
                        nr_homs = current_block_homs,
                        constant = block_constant,
                        block_size = block_size,
                        max_hets = max_hets,
                        max_het_fraction = max_het_fraction
                    )
                    
            
                    if block_ends:
                        valid_block = check_valid_block(
                            block_size=block_size, 
                            shortest_block_tres=shortest_block, 
                            nr_homs=current_block_homs, 
                            minimum_homs_tres=minimum_homs
                        )
                        if valid_block:
                            end_block(current_chr, hom_block_start, hom_block_end, current_block_hets, current_block_homs)
            
                    chr_blocked += block_size
            
                    current_block_homs = 0
                    current_block_hets = 0
                    hom_block_start = None
                    hom_block_end = None
                    seen_hom = False
            
                else:
                    logger.info("Block ends before one was started at %s %d." %
                                (var.CHROM, var.start))
            
                block_ends = False
            
            if chr_ends:
                chr_size_spanned = chr_last_pos_seen
                end_chr(current_chr, chr_size_spanned, chr_blocked, chr_hets, chr_homs,
                        flag_UPD_at_fraction)
            
                chr_homs = 0
                chr_hets = 0
                current_chr = var.CHROM
                chr_blocked = 0
                chr_ends = False 
            
            if start_new_block:
                # new block; double check homoz for new chr.. just ignore if het.
                if check_homozygote(var, proband_vcf, individual):
                    seen_hom = True
                    current_block_homs = 1
                    hom_block_start = var.start
                    hom_block_end = var.end
                    current_chr = var.CHROM
            
                start_new_block = False
            
            # finally, update chr totals and then return to beginning of loop to draw new var
            if check_homozygote(var, proband_vcf, individual):
                chr_homs += 1
            elif check_heterozygote(var, proband_vcf, individual):
                chr_hets += 1
            
            chr_last_pos_seen = var.end

    if (hom_block_end and hom_block_start):
    # and check for block once more after seeing last var
        block_size = hom_block_end - hom_block_start + 1
        block_ends = check_block_end(
            nr_hets = current_block_hets,
            nr_homs = current_block_homs,
            constant = block_constant,
            block_size = block_size,
            max_hets = max_hets,
            max_het_fraction = max_het_fraction
        )
        
        if block_ends:
            valid_block = check_valid_block(
                block_size=block_size, 
                shortest_block_tresh=shortest_block, 
                nr_homs=current_block_homs, 
                minimum_homs_tres=minimum_homs
            )
            
            if valid_block:
                end_block(current_chr, hom_block_start, hom_block_end, current_block_hets, block_homs)
    
