"""
excision.py

Goal:
    Looks for Cases where insertion has be excised in one sample.
    Find TE insertion sites homozygous in n-1, where one sample has a frequency != 1
- possible duplicate methods with polyFilter?
"""

import sys
from polyLine import polyLine
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import time
from scipy.misc import imread
from scipy.misc import imresize
import matplotlib.image as mpimg
import subprocess
import re
from pyBinom import pbinom
from TE_plots import *
from polyFilter import *

CC = ["CC_A","CC_B","CC_C","CC_D","CC_E","CC_F","CC_G","CC_H","CC_I","CC_J","CC_K","CC_L","CC_M","CC_N","CC_O","CC_P"]
GP = ["GP2-3_A","GP2-3_B","GP2-3_C","GP2-3_D","GP2-3_E","GP2-3_F","GP2-3_G","GP2-3_H","GP2-3_I","GP2-3_J","GP2-3_K","GP2-3_L","GP2-3_M","GP2-3_N","GP2-3_O","GP2-3_P"]
#Excludes samples G and M from CC genotype
CC_mod = ["CC_A","CC_B","CC_C","CC_D","CC_E","CC_F","CC_H","CC_I","CC_J","CC_K","CC_L","CC_N","CC_O","CC_P"]

both_genotypes = CC_mod + GP
target_genotype = both_genotypes
n = 1000

def outlier_freq(te_sites):
    """
    Returns the # of samples that do not have a freq of 1
        Excision represented by n - 1 samples @ frequency of 1 & one sample at freq 0.5 or 0.0
    """
    #Dict of {"CC_A":[TE_site, TE_site ...],"CC_B":[...],...}
    samples = {}

    for site in te_sites:
        samples[site.sample] = samples.get(site.sample,[])
        samples[site.sample].append(site)

    #Samples at the site
    num_samples = len(samples.keys())

    #num_samples MUST be equal to len(target_genotype)
    if num_samples != len(target_genotype):
        print "Error num sample not global"

    freq_1_samples = 0
    not_freq_1_samples = 0

    for tar_sample in samples.keys():
        sample_freq = max([site.freq for site in samples[tar_sample]])

        if sample_freq == 1.0:
            freq_1_samples += 1
        else:
            not_freq_1_samples += 1
   
    #excision case is where not_freq_1_samples = 1 and freq_1_samples = len(target_genotype) - 1
    return not_freq_1_samples 


def sample_depth_filter(sample_sites):
    """
    For all sites present for one sample in a single 1000bp window
    If at least one of the sites passes depth filter ==> sample in this window passes
    """
    sample_depth = True

    min_depth = 3
    sample_depth = True

    for site in sample_sites:
        if site.total_TE_depth < min_depth:
            sample_depth = False    

    return sample_depth 


def window_depth_filter(window):
    """
    For a given 1000bp window,
    Return True if every sample present passes depth filters

    """    
    depth_pass = True

    samples = {}
    
    for site in window:
        samples[site.sample] = samples.get(site.sample,[])
        samples[site.sample].append(site)
    
    for sample in samples.keys():
        if sample_depth_filter(samples[sample]) != True:
            depth_pass = False
    
    return depth_pass


def filter_global_for_outlier(TE_dict):
    """
    Filters global TE sites, TE present across all samples in any freq, for sites
        where n - 1 samples have a insertion frequency of 1.
    """              
    filt_TE_dict = {}

    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom].keys():
            if window_depth_filter(TE_dict[chrom][pos_range]):
                #Window passed depth filters

                if outlier_freq(TE_dict[chrom][pos_range]) == 1:
                    #Window has one outlier of low freq 
                    filt_TE_dict[chrom] = filt_TE_dict.get(chrom,{})
                    filt_TE_dict[chrom][pos_range] = TE_dict[chrom][pos_range]

    return filt_TE_dict


if __name__ == "__main__":
    TE_pos = get_TE_dict(target_genotype,n) 
    TE_pos = filter_abs_ranges(TE_pos)
    global_sites = get_global_ranges(TE_pos)    

    filt_global = filter_global_for_outlier(global_sites)
    print "# global_sites: ",  count_windows(global_sites)
    print "# filt_global_sites: ",  count_windows(filt_global)

    write_TE_file(filt_global,"excision_global_sites.txt")









