"""
excl_global.py
    - partial clone of polyFilter
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
import operator

CC = ["CC_A","CC_B","CC_C","CC_D","CC_E","CC_F","CC_G","CC_H","CC_I","CC_J","CC_K","CC_L","CC_M","CC_N","CC_O","CC_P"]
GP = ["GP2-3_A","GP2-3_B","GP2-3_C","GP2-3_D","GP2-3_E","GP2-3_F","GP2-3_G","GP2-3_H","GP2-3_I","GP2-3_J","GP2-3_K","GP2-3_L","GP2-3_M","GP2-3_N","GP2-3_O","GP2-3_P"]

#Excludes samples G and M from CC genotype
CC_mod = ["CC_A","CC_B","CC_C","CC_D","CC_E","CC_F","CC_H","CC_I","CC_J","CC_K","CC_L","CC_N","CC_O","CC_P"]

both_genotypes = CC_mod + GP
target_genotype = both_genotypes

#Remember calculated Pvals to speed finding Pvals
n = 1000


def sample_depth_filter(sample_sites):
    """
    Given list of TE sites for one sample in a TE window,
        Return True if at least one site passes minimum depth filters
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
    For TE insertions in a given 1000bp window,
    Return True if every sample present passes depth filters
    """    
    depth_pass = True
    samples = {}
   
    #Load all window sites into dict of samples: sites 
    for site in window:
        samples[site.sample] = samples.get(site.sample,[])
        samples[site.sample].append(site)
   
    #Make sure each sample has atleast one site which passes depth filter 
    for sample in samples.keys():
        if sample_depth_filter(samples[sample]) != True:
            depth_pass = False
    
    return depth_pass


def filter_global(TE_dict):
    """
    Given TE_dict set of global TE sites,
        Return sites where all samples pass minimum depth filter.
    """              
    filt_TE_dict = {}

    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom].keys():
            if window_depth_filter(TE_dict[chrom][pos_range]):
                #Window passed depth filters
                filt_TE_dict[chrom] = filt_TE_dict.get(chrom,{})
                filt_TE_dict[chrom][pos_range] = TE_dict[chrom][pos_range]

    return filt_TE_dict



def get_excl_sites(target_dict,rem_dict):
    """
    Returns TE_dict = target_dict - rem_dict
    """
    target_dict_excl = {}

    #Check for sites in target_dict and not in rem_dict
    for chrom in target_dict.keys():

        #Chromosome missing from rem_dict
        if rem_dict.get(chrom,None) == None:
            target_dict_excl[chrom] = target_dict_excl.get(chrom,{})
            target_dict_excl[chrom] = CC[chrom]

        else:
            for pos_range in target_dict[chrom]:
                if rem_dict[chrom].get(pos_range,None) == None:
                    #pos window not in rem_dict
                    target_dict_excl[chrom] = target_dict_excl.get(chrom,{})
                    target_dict_excl[chrom][pos_range] = target_dict[chrom][pos_range]
                     
    return target_dict_excl 


def hetz_filter(window):
    """
    
    """

    freq_pass = True
    samples = {}
   
    #Load all window sites into dict of samples: sites 
    for site in window:
        samples[site.sample] = samples.get(site.sample,[])
        samples[site.sample].append(site)
   
    #Make sure each sample has atleast one site which passes freq filter 
    for sample in samples.keys():
        
        #Need at least one site for this sample with freq != 1.0
        site_pass = False
        for site in samples[sample]:
            if site.freq != 1.0:
                site_pass = True 

        #One sample failed - no sites for this sample that are not homz.
        if site_pass == False:
            freq_pass = False           
  
    return freq_pass


def get_hetz_sites(TE_dict):
    """
    Duplicate?
    """
    hetz_TE_dict = {}

    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom].keys():
            if hetz_filter(TE_dict[chrom][pos_range]):
                #Window passed depth filters
                hetz_TE_dict[chrom] = hetz_TE_dict.get(chrom,{})
                hetz_TE_dict[chrom][pos_range] = TE_dict[chrom][pos_range]

    return hetz_TE_dict

def write_TE_file(TE_dict,filename):
    """
    DUPLICATE
    writes TE_dict out to filename
    """
    file_out = open(filename,"w")
    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom].keys():
            for pLine in TE_dict[chrom][pos_range]:
                out_line = pLine.sample + "\t" + pLine.raw_line
                #out_line = pLine.raw_line
                file_out.write(out_line)    
            file_out.write('\n')

    file_out.close()

def get_family_counts(TE_dict):
    """
    Returns family count dict {LTR:count,unknown:count ... } given TE_dict.
    """
    families = {}
    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom].keys():
            family_site_count = {}
            for site in TE_dict[chrom][pos_range]:
                family_site_count[site.family] = families.get(site.family,0)
                family_site_count[site.family] += 1

            pos_family  = max(family_site_count.iteritems(), key=operator.itemgetter(1))[0]             
            families[pos_family] = families.get(pos_family,0)
            families[pos_family] += 1

    return families

if __name__ == "__main__":

   #        FOR GLOBAL ACROSS BOTH GENOTYPES
    all_pos = get_TE_dict(target_genotype,n)
    all_pos = get_ranges_by_sample_count(all_pos,len(target_genotype) - 1)
#    all_pos = get_hetz_sites(all_pos)
    all_pos = filter_global(all_pos)
    write_TE_file(all_pos,"./genotype_global/glob_hom_allpos.txt")

#    plot_freq_hist(all_pos,"./genotype_global/noglob_hom_all_pos_freq.png")
#    plot_depth_hist(all_pos,"./genotype_global/noglob_hom_all_pos_depth.png")
 
    print "hetz_all_pos",count_windows(all_pos)


#
#
#

