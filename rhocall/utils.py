from __future__ import division

import logging

logger = logging.getLogger(__name__)

def skip_variant(variant):
    """Check if a variant should be skipped
    
        Args:
            variant(cyvcf2.Variant)
        
        Returns:
            bool: if variant should be skipped
    """
    # ignore low qual and indels for now
    if variant.FILTER:
        logger.debug("Filter %s" % variant.FILTER)
        return True

    if variant.end - variant.start > 1:
        logger.debug("Indel size %d" % (variant.end - variant.start))
        return True
    
    return False
    

def check_homozygote(variant, vcf, ind_index):
    """Check if a variant is homozygote
    
        Args:
            variant(cyvcf2.Variant)
            vcf(cyvcf2.VCF)
            ind_index(int)
    
        Returns:
            bool
    """
    
    if variant.gt_types[ind_index] == vcf.HOM_REF:
        return True

    if variant.gt_types[ind_index] == vcf.HOM_ALT:
        return True
    
    return False

def check_heterozygote(variant, vcf, ind_index):
    """Check if a variant is heterozygote
    
        Args:
            variant(cyvcf2.Variant)
            vcf(cyvcf2.VCF)
            ind_index(int)
    
        Returns:
            bool
    """
    
    if variant.gt_types[ind_index] == vcf.HET:
        return True

    return False

def check_block_end(nr_hets, nr_homs, constant, block_size, max_hets, max_het_fraction):
    """Check if the fraction of hets aborts the block"""
    
    if (nr_hets * constant / block_size) > max_hets:
        return True
    
    if (nr_hets / nr_homs) > max_het_fraction:
        return True
    
    return False

def check_valid_block(block_size, shortest_block_tres, nr_homs, minimum_homs_tres):
    """Check if a block is of valid size
    
        Args:
            block_size(int): Size of block
            shortest_block_tres(int): User defined minimum block size
            nr_homs(int): Nr observed homozygotes in block
            minimum_homs_tres(int): User defined minimum homs for one block
        
        Returns:
            bool: If block is valid
    """
    if block_size > shortest_block_tres:
        if nr_homs > minimum_homs_tres:
            return True
    
    return False



