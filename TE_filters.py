"""
TE_filters.py

Various filters for TE sites generated by popoolationTE from MA data

Filters on:
 - Distance from masked
 - Depth
 - Chromosome
 - Absence range

 - Singleton
    - Binomial
    - TE depth
    - window edge

 - Global



IMPLEMENTED in TE_filters.py (May 4th 2016)

 - DISTANCE FROM MASKED SEQUENCE
 - HARD WINDOW EDGE CORRECTION FILTER
 - FILTERING BY # OF SAMPLES IN WINDOW - global and single
 - WINDOW QUALITY FILTER
    - depth
    - binomial
  
"""
import sys
from polyLine import polyLine
from pyBinom import pbinom

### BEGIN - DISTANCE FROM MASKED SEQUENCE ###
def dist_masked_filter(site, masked_ranges, min_dist):
    """
    Returns true if given site is at least min_dist bp away from any masked_ranges.
    """
    chrom = site.chrom
    pos = int(site.pos)
    right_pos = pos + int(min_dist)
    left_pos = pos - int(min_dist)
    far_from_masked = True
   
    #For all possible points within min_dist of pos 
    for point in range(left_pos,right_pos):

        #Check if any points fall within a masked range
        for (low,high) in masked_ranges[chrom]:
            if point <= int(high) and point >= int(low):
                far_from_masked = False
                return far_from_masked

    return far_from_masked 
### END - DISTANCE FROM MASKED SEQUENCE ###


### BEGIN - FILTER EDGE CORRECTION ###
def get_samples_in_range(TE_pos,target_site,range_size):
    """
    Given a target_site, find any TEs present within 500bp (over window edges) in any other sample.
    """ 
    n = range_size*2
    chrom = target_site.chrom
    pos = target_site.pos
    left_pos, right_pos = None, None 
    new_window_sites = []
    
    #Get pos_range
    pos_idx = int(pos/n)

    #Start position of window
    window_start = pos_idx*n
    
    #Get position within window for target_site
    window_pos = int(pos) - int(window_start)
   
    # Target site in right side of window, more than 500bp from the left edge 
    #    ==> know there are no TEs in other samples within this window 
    # Only have  to check the right window
    if window_pos > n/2:
        #Want to check range pos:pos+500
        pos_right_border = pos + range_size
    
        #Get new pos_range idx 
        pos_idx = int(pos_right_border/n)
        window = TE_pos[chrom].get(pos_idx,None)

        #Find sites within neighbouring window which are within 500bp of target_site
        if window != None:
            for site in window:
               if site.pos < pos_right_border:
                    new_window_sites.append(site)
        
    #Only check the left window -- same methods as right window
    if window_pos < n/2:
        pos_left_border = pos - range_size
        pos_idx = int(pos_left_border/n)
        window = TE_pos[chrom].get(pos_idx,None)

        if window != None:
            for site in window:
                if site.pos > pos_left_border:
                    new_window_sites.append(site)

    return new_window_sites

def singleton_edge_filter(all_TE_set, site, window_size):
    """
    Filter for sites with no trace of a TE within window_size/2 bp (Across window edges)
    """
    site_pass = True
    half_window = window_size/2
    
    #Get list of TEs within half_window bp of given site
    new_window_TEs = get_samples_in_range(all_TE_set, site, half_window)
   
    #Look for any new TEs not from the same sample as focal site 
    for TE in new_window_TEs:

        #Insert present in other sample within half_window bp of site --> site fails to pass filter
        if TE.sample != site.sample:
            site_pass = False
            return site_pass
    
    return site_pass

### END FILTER EDGE CORRECTION ###


### BEGIN - FILTER FOR TE-WINDOWS BY # SAMPLES PRESENT ###
def filter_window_by_sample_count(window,sample_count):
    """
    Return True for windows which contain sample_count samples

    sample_count == 1 --> Returns singletons
    sample_count == #Samples in analysis --> global TEs
    """
    window_pass = False

    if num_sample_in_window(window) == sample_count:
        window_pass = True

    return True

def num_sample_in_window(window):
    """
    window -> list of pLine objects in TE_pos[chrom][window_idx] -->  window
    """
    samples_present = set() 
    for pline in window:
        samples_present.add(pline.sample) 
    
    return len(samples_present)
### END - FILTER FOR TE-WINDOWS BY #Samples Present ###


### BEGIN - WINDOW QUALITY FILTER ###
def passTEDepth(TE_site, min_depth, max_depth):
    """
    Requires that TE read count > min_depth and < max_depth
    """
    altReads = TE_site.TE_reads_forward + TE_site.TE_reads_reverse
    return (altReads >= min_depth) and (altReads < max_depth)

def getPval(altReads, refReads):
    """
    Get Pval of TE site assuming heterozygous for TE (p=0.5)
    """
    return pbinom(altReads,refReads,0.5)

def low_passBinom(TE_site, min_prob):
    """ 
    Return True if TE_site passes binomial test (prob > min_prob (0.02)).
    """
    altReads = TE_site.TE_reads_forward + TE_site.TE_reads_reverse
    refReads = TE_site.nonTE_reads_forward + TE_site.nonTE_reads_reverse
    binomResult = (getPval(altReads, refReads) > min_prob)
    return (binomResult)


def filter_window_quality(window, min_TE_depth, max_TE_depth, binom_min_prob):
    """
    Given a TE window, return True if at least one site in the window passes binomial and depth filters.
    """
    min_TE_depth = 3
    max_TE_depth = 30
    binom_min_prob = 0.02

    TE_depth_pass = False
    binom_pass = False

    for site in window:
         if low_passBinom(site, binom_min_prob) and passTEDepth(site,min_TE_depth,max_TE_depth):
            TE_depth_pass = True
            binom_pass = True

    return (binom_pass and TE_depth)

### END - WINDOW QUALITY FILTER ###



























