"""
polyFilter.py


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
import operator

CC = ["CC_A","CC_B","CC_C","CC_D","CC_E","CC_F","CC_G","CC_H","CC_I","CC_J","CC_K","CC_L","CC_M","CC_N","CC_O","CC_P"]
GP = ["GP2-3_A","GP2-3_B","GP2-3_C","GP2-3_D","GP2-3_E","GP2-3_F","GP2-3_G","GP2-3_H","GP2-3_I","GP2-3_J","GP2-3_K","GP2-3_L","GP2-3_M","GP2-3_N","GP2-3_O","GP2-3_P"]
#Excludes samples G and M from CC genotype
CC_mod = ["CC_A","CC_B","CC_C","CC_D","CC_E","CC_F","CC_H","CC_I","CC_J","CC_K","CC_L","CC_N","CC_O","CC_P"]
wo_AP_CC = ["CC_B","CC_C","CC_D","CC_E","CC_F","CC_H","CC_I","CC_J","CC_K","CC_L","CC_N","CC_O"]
wo_AP_GP = ["GP2-3_B","GP2-3_C","GP2-3_D","GP2-3_E","GP2-3_F","GP2-3_G","GP2-3_H","GP2-3_I","GP2-3_J","GP2-3_K","GP2-3_L","GP2-3_M","GP2-3_N","GP2-3_O"]

both_genotypes = CC_mod + GP
wo_AP_both = wo_AP_CC + wo_AP_GP
target_genotype = both_genotypes
#Remember calculated Pvals to speed finding Pvals
Pval_dict = {}
n = 1000
poly_dir = "./poly_files/"

def loadTEranges(TE_file_loc):
    """
    Given location of TE ranges returns dictionary of TE ranges.
    Load TE ranges from a text file of the format CHROM START STOP.
    """
    TE_ranges = {}
    CHROM,START,STOP = 0,1,2

    with open(TE_file_loc) as TE_file:
        for line in TE_file:
            line_col = str.split(line)
            TE_ranges[line_col[CHROM]] = TE_ranges.get(line_col[CHROM],[])
            TE_ranges[line_col[CHROM]].append((line_col[START],line_col[STOP]))

    TE_file.close()
    return TE_ranges


def get_valid_ranges(TE_dict,TE_ranges):
    """
    Given TE_dict and masked ranges TE_ranges, 
        return sites from TE_dict that pass masked distance filter from TE_ranges sites.
    """
    valid_TE = {}
    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom].keys():
            for site in TE_dict[chrom][pos_range]:
                if far_from_masked(site,TE_ranges):
                    valid_TE[chrom] = valid_TE.get(chrom,{})
                    valid_TE[chrom][pos_range] = valid_TE[chrom].get(pos_range,[])
                    valid_TE[chrom][pos_range].append(site) 
    return valid_TE


def far_from_masked(site,TE_ranges):
    """
    Return True if given site is at least dist bp away from any masked range in TE_ranges.
    """
    dist = 150
    chrom  = site.chrom
    pos = site.pos
    right_pos = pos + dist
    left_pos = pos - dist
    near = False

    #TE range = left_pos --> high_pos
    for point in range(int(left_pos),int(right_pos)):
        for (low,high) in TE_ranges[chrom]:
            if point <= int(high) and point >= int(low):
                near = True

    return not near



def get_all_TE_depth(TE_dict):
    """
    Returns list of depth of all TEs in TE_dict    
    """
    depths = []
    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom]:
            for site in TE_dict[chrom][pos_range]:
                depths.append(site.total_TE_depth)
    return depths


def get_samples_in_range(TE_pos,target_site,range_size):
    """
    Given a target_site, find any TEs present within 500bp (over window edges) in any other sample.
    """ 
    global n
    chrom = target_site.chrom
    pos = target_site.pos
    left_pos, right_pos = None, None 

    new_window_sites = []
    pos_idx = int(pos/n)
    window_start = pos_idx*n
    window_pos = int(pos) - int(window_start)
    
    # Only have  to check the right window
    if window_pos > n/2:
        #Want to check range pos:pos+500
        pos_right_border = pos + range_size
        pos_idx = int(pos_right_border/n)
        window = TE_pos[chrom].get(pos_idx,None)

        if window != None:
            for site in window:
               if site.pos < pos_right_border:
                    new_window_sites.append(site)
        
    #Only check the left window  
    if window_pos < n/2:
        pos_left_border = pos - range_size
        pos_idx = int(pos_left_border/n)
        window = TE_pos[chrom].get(pos_idx,None)

        if window != None:
            for site in window:
                if site.pos > pos_left_border:
                    new_window_sites.append(site)

    return new_window_sites


def filt_singleton_edge(TE_pos, singletons):
    """
    Fix problem with windows causing false positive singletons
    """
    global n
    window_width = n/2
    new_singletons = {}

    #Loop through every singleton TE read
    for chrom in singletons.keys():
        for site_range in singletons[chrom].keys():
            for singleton_site in singletons[chrom][site_range]:
            
                #Get samples with TEs present in the new window 500 bp on either side of TE
                new_window_samples = get_samples_in_range(TE_pos,singleton_site,window_width)
                true_singleton = True

                for TE in new_window_samples:
                    if TE.sample != singleton_site.sample:
                        #TE of other sample in window
                        true_singleton = False
               
                #Singleton passed filters 
                if true_singleton:
                    new_singletons[chrom] = new_singletons.get(chrom,{})                     
                    new_singletons[chrom][site_range] = new_singletons[chrom].get(site_range,[])
                    new_singletons[chrom][site_range].append(singleton_site) 

    return new_singletons


def num_sample_in_window(window):
    """
    window -> list of pLine objects in TE_pos[chrom][window_idx] -->  window
    """
    samples_present = set() 
    for pline in window:
        samples_present.add(pline.sample) 
    
    return len(samples_present)


def get_singleton_ranges(TE_pos):
    """
    Find TE windows unique to one sample
    """
    TE_singletons = {}
    
    for chrom in TE_pos.keys():
        for pos_range in TE_pos[chrom].keys():
            if num_sample_in_window(TE_pos[chrom][pos_range]) == 1:
                TE_singletons[chrom] = TE_singletons.get(chrom,{})
                TE_singletons[chrom][pos_range] = TE_pos[chrom][pos_range]

    return TE_singletons


def low_passBinom(TE_site, min_prob):
    """ 
    Return True if TE_site passes binomial test (prob > min_prob (0.02)).
    """
    altReads = TE_site.TE_reads_forward + TE_site.TE_reads_reverse
    refReads = TE_site.nonTE_reads_forward + TE_site.nonTE_reads_reverse
    binomResult = (getPval(TE_site) > min_prob)
    return (binomResult)


def high_passBinom(TE_site, min_prob):
    """ 
    Return True if TE site passes high binomial test (1 - prob) > 0.02
        - Possibly to filter for true heterozygotes    
    """
    altReads = TE_site.TE_reads_forward + TE_site.TE_reads_reverse
    refReads = TE_site.nonTE_reads_forward + TE_site.nonTE_reads_reverse
    binomResult = (1 - getPval(TE_site) > min_prob)
    return (binomResult)


def passDepth(TE_site, min_depth):
    """
    Requires min_depth TE reads
    """
    altReads = TE_site.TE_reads_forward + TE_site.TE_reads_reverse
    refReads = TE_site.nonTE_reads_forward + TE_site.nonTE_reads_reverse
    return (altReads + refReads) > min_depth


def passTEDepth(TE_site, min_depth, max_depth):
    """
    Requires that TE read count > min_depth and < max_depth
    """
    altReads = TE_site.TE_reads_forward + TE_site.TE_reads_reverse
    return (altReads >= min_depth) and (altReads < max_depth)


def getPval(TE_site):
    """
    Get Pval of TE site assuming heterozygous for TE (p=0.5)
    """
    altReads = TE_site.TE_reads_forward + TE_site.TE_reads_reverse
    refReads = TE_site.nonTE_reads_forward + TE_site.nonTE_reads_reverse
    return pbinom(altReads,refReads,0.5)


def getAllPvals(TE_pos):
    """
    Return a list of P values from cumulative binomial test of TE vs. ref reads
    """
    Pvals = []
    for chrom in TE_pos.keys():
        for pos in TE_pos[chrom].keys():
            for site in TE_pos[chrom][pos]:
                altReads = site.TE_reads_forward + site.TE_reads_reverse
                refReads = site.nonTE_reads_forward + site.nonTE_reads_reverse
                read_key = (altReads, refReads)
    
                #If pval not been found before calc it
                if Pval_dict.get(read_key,None) == None:
                    Pval_dict[read_key] = getPval(site)

                Pvals.append(Pval_dict[read_key])

    return Pvals


def get_global_ranges(TE_pos):
    """
    Generates a dict of sites where the TE is present in every sample
    """
    global_TE = {}
    for chrom in TE_pos.keys():
        for pos_range in TE_pos[chrom].keys():
            if num_sample_in_window(TE_pos[chrom][pos_range]) == len(target_genotype):
                global_TE[chrom] = global_TE.get(chrom,{})
                global_TE[chrom][pos_range] = TE_pos[chrom][pos_range]

    return global_TE 


def get_ranges_by_sample_count(TE_pos,sample_count):
    """
    Generates a dict of sites where a TE is present in sample_count samples.

    sample_count == 1 --> Returns singletons
    sample_count == #Samples in analysis --> global TEs
    """
    TEs_present = {}

    for chrom in TE_pos.keys():
            for pos_range in TE_pos[chrom].keys():
                if num_sample_in_window(TE_pos[chrom][pos_range]) == sample_count:
                    TEs_present[chrom] = TEs_present.get(chrom,{})
                    TEs_present[chrom][pos_range] = TE_pos[chrom][pos_range]

    return TEs_present

def window_filter(window):
    """
    Given a TE window, return True if atleast one site in the window passes binomial and depth filters.
    """
    min_TE_depth = 3
    max_TE_depth = 30
    min_prob = 0.02

    forward = False
    Reverse = False
    Coverage = 0
    family = None

    #IF ONE site in the window passes the binomial test --> window passes 
    # IF One site (not necessarily the same) passes the binomial test --> window passes
    #       What would a site be like with the sites being passed from exclusive sites?
    #      --  Make it the same site to pass these tests?
    TE_depth_pass = False
    binom_pass = False
    for site in window:
         if low_passBinom(site, min_prob) and passTEDepth(site,min_TE_depth,max_TE_depth):
            TE_depth_pass = True
            binom_pass = True

    return binom_pass  


def get_filtered_ranges(TE_pos):
    """
    Given TE_dict return new TE dict of TE windows which pass window_filter().

    """
    filtered_TEs = {}
    
    for chrom in TE_pos.keys():
        for pos_range in TE_pos[chrom].keys():
            if window_filter(TE_pos[chrom][pos_range]):
                #Window passed filters
                filtered_TEs[chrom] = filtered_TEs.get(chrom,{})
                filtered_TEs[chrom][pos_range] = TE_pos[chrom][pos_range]

    return filtered_TEs


def count_windows(TE_pos):
    """
    Return the number of TE windows present in a TE dict.
    """
    count = 0
    for chrom in TE_pos.keys():
        count += len(TE_pos[chrom].keys())    

    return count

def count_sites(TE_pos):
    """
    Return the number of sites present in a TE dict.
        Multiple sites per TE window.
    """
    count = 0
    for chrom in TE_pos.keys():
        for pos_range in TE_pos[chrom]:
            count += len(TE_pos[chrom][pos_range])    

    return count
   
 
def write_TE_file(TE_dict,filename):
    """
    Write out TE dict to given filename.
    """
    file_out = open(filename,"w")
    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom].keys():
            for pLine in TE_dict[chrom][pos_range]:
                out_line = pLine.sample + "\t" + pLine.raw_line
                #out_line = pLine.raw_line
                file_out.write(out_line)    

    file_out.close()


def init_TE_dict():
    """
    Initializes a TE dict with chromosome keys.
    """
    TE_pos = {"mitochondrion":{},"chloroplast":{}} 
    for i in range(0,33):
        chrom_str = "pseudo" + str(i)
        TE_pos[chrom_str] = TE_pos.get(chrom_str,{})

    return TE_pos


def get_TE_dict(file_name_list,window_size):
    """
    Given list of filenames, generates a TE dict of the format:
        {chrom:{pos_range:[polyLine,polyLine, ...],pos_range:[polyLine,polyLine,...], ...}, ...}

    Removes 'bad' chromosomes: pseudo0, mitochondrion, chloroplast
    """
    global poly_dir
    TE_pos = init_TE_dict()

    for file_name in  (file_name_list):
            with open(poly_dir + file_name) as file:
                for line in file:
                    #Init polyLine object representing one line in the polymorphism file
                    pLine = polyLine(line,file_name)
                    #Get the postion key
                    pos_key = int(float(pLine.pos)) / window_size
    
                    #print pos_key
                    TE_pos[pLine.chrom][pos_key] = TE_pos[pLine.chrom].get(pos_key,[])
                    TE_pos[pLine.chrom][pos_key].append(pLine)

    bad_chroms = ["pseudo0","mitochondrion","chloroplast"]
    for key in bad_chroms:
        TE_pos.pop(key,None)

    return TE_pos


def get_all_freq(TE_pos):
    """
    Return list of all insertion frequencies from provided TE dict.
    """
    freq = []
    for chrom in TE_pos.keys():
        for pos in TE_pos[chrom].keys():
            for site in TE_pos[chrom][pos]:
                freq.append(site.freq)
    return freq


def get_directional_info(TE_dict):
    """
    Returns directional info of insertions in provided TE dict.
    """
    dir_counts = {"F":0,"R":0,"FR":0,"RF":0}
    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom].keys():
            for site in TE_dict[chrom][pos_range]:
                dir_counts[site.direction] = dir_counts.get(site.direction,0)
                dir_counts[site.direction] = dir_counts[site.direction] + 1

    return dir_counts


def get_high_freq(TE_dict):
    """
    Returns TE dict of TE windows where at least one site has a frequency of 1.0.
    """
    high_freq = {}
    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom].keys():
            for site in TE_dict[chrom][pos_range]:
                if site.freq == 1.0:
                    passed = True
                    high_freq[chrom] = high_freq.get(chrom,{})
                    high_freq[chrom][pos_range] = high_freq[chrom].get(pos_range,[])
                    high_freq[chrom][pos_range].append(site)

    return high_freq

  
def get_low_freq(TE_dict):
    """
    Given TE dict returns windows where all sites do not have a frequency of 1.0.
    """
    high_freq = {}
    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom].keys():
            for site in TE_dict[chrom][pos_range]:
                if site.freq != 1.0:
                    high_freq[chrom] = high_freq.get(chrom,{})
                    high_freq[chrom][pos_range] = high_freq[chrom].get(pos_range,[])
                    high_freq[chrom][pos_range].append(site)

    return high_freq

 
def get_range_dists(TE_dict):
    """
    Returns a list of the bp difference from start and end range
    """
    range_sizes = []
    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom].keys():
            diff = None
            for site in TE_dict[chrom][pos_range]:

                if site.direction == "R":
                    diff = int(site.start_range_reverse) - int(site.end_range_reverse)

                if site.direction == "F":
                    diff = int(site.start_range_forward) - int(site.end_range_forward)

                if site.direction == "FR":
                    diff_F = int(site.start_range_forward) - int(site.end_range_forward)
                    diff_R = int(site.start_range_reverse) - int(site.end_range_reverse)
                    diff = (diff_R + diff_F)/2.0

                range_sizes.append(diff)
    return range_sizes


def range_filter(site, min_range):
    """
    Returns True if given site has a TE-absence range of at least min_range.
    """
    diff = None

    if site.direction == "R":
        diff = int(site.start_range_reverse) - int(site.end_range_reverse)

    if site.direction == "F":
        diff = int(site.start_range_forward) - int(site.end_range_forward)

    if site.direction == "FR":
        diff_F = int(site.start_range_forward) - int(site.end_range_forward)
        diff_R = int(site.start_range_reverse) - int(site.end_range_reverse)
        diff = max(diff_F,diff_R)

    return abs(diff) >= min_range


def filter_abs_ranges(TE_dict):
    """
    Returns TE sites with valid TE-absence ranges. 
    """
    passed_ranges = {}
    min_range = 73
    for chrom in TE_dict.keys():
       for pos_range in TE_dict[chrom].keys():
           for site in TE_dict[chrom][pos_range]:
                if range_filter(site,min_range):
                    passed_ranges[chrom] = passed_ranges.get(chrom,{})
                    passed_ranges[chrom][pos_range] = passed_ranges[chrom].get(pos_range,[])
                    passed_ranges[chrom][pos_range].append(site)

    return passed_ranges

 
def get_family_counts(TE_dict):
    """ 
    Returns family count information of TE insertions.
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
    #!Base pair window size
    
    # ----------------------------------------------
    #CC_pos = get_TE_dict(CC_mod,n)
    #GP_pos = get_TE_dict(GP,n)
    #print "CC_pos count", count_windows(CC_pos) 
    #print "GP_pos count", count_windows(GP_pos) 
    #CC_range = get_ranges_by_sample_range(CC_pos,range(0,8))
    #GP_range = get_ranges_by_sample_range(GP_pos,range(0,8))
    #print "CC_global count", count_windows(CC_range)
    #print "GP_global count", count_windows(GP_range) 
    #mut_excl_global = get_excl_sites(CC_range,GP_range)

    #mut_excl_freq = get_all_freq(mut_excl_global)
    #plot_hist(mut_excl_freq,"Estimated TE frequency at site","mut_excl_freq_global_te_sites.png")
    ## -----------------------------------------------
    
    TE_pos = get_TE_dict(target_genotype,n)
    print "TE_pos",count_windows(TE_pos)
    print "TE_pos",count_sites(TE_pos)
    TE_pos = filter_abs_ranges(TE_pos)
    print "TE_pos_after",count_windows(TE_pos)
    print "TE_pos_after",count_sites(TE_pos)


    #REMOVE SITES WITHIN 50BP OF MASKED SEQUENCE 
#    masked_ranges = remTE_range.loadTEranges("masked_ranges/incl_small_50bp+_range_parsed.txt")      
#    print masked_ranges.keys()
#    TE_pos_valid = remTE_range.get_valid_ranges(TE_pos,masked_ranges)
#    write_TE_file(TE_pos_valid,"TE_pos_valid.txt")
#    print "valid", count_windows(TE_pos_valid)




    plot_count_by_samples_in(TE_pos,target_genotype,"num_TE_sample_count.png")

    single_sites = get_singleton_ranges(TE_pos)
    global_sites = get_global_ranges(TE_pos)

    plot_depth_hist(single_sites,"no_AP_final_unfilt_singletons_depth.png")
    plot_freq_hist(single_sites,"no_APfinal_unfilt_singletons_freq.png")

    print "single", count_windows(single_sites)
    print "singleton dir", get_directional_info(single_sites)
    print "global ranges: " + str(count_windows(global_sites))
    print "global dir", get_directional_info(global_sites)

    filter_singletons = get_filtered_ranges(single_sites)
    print "filtered with low freq and by min & max depth"
    print "filtered singletons: " + str(count_windows(filter_singletons))
    print "filtered dir", get_directional_info(filter_singletons)

    #New filter edge cases
    window_filt_single = filt_singleton_edge(TE_pos,filter_singletons)
    
    #Apply absence range cutoff
    window_filt_single = filter_abs_ranges(window_filt_single)

    window_filt_single = get_filtered_ranges(window_filt_single)


#    masked_ranges = remTE_range.loadTEranges("masked_ranges/incl_small_50bp+_range_parsed.txt")      
#    window_filt_single = remTE_range.get_valid_ranges(window_filt_single,masked_ranges)
  
#HERE
    print "size before",count_windows(window_filt_single) 
    masked_ranges = loadTEranges("masked_ranges/no_small_gtf_parsed.txt")    
    window_filt_single = get_valid_ranges(window_filt_single,masked_ranges)    
    print "size after",count_windows(window_filt_single) 

    print "window filter singletons ", count_windows(window_filt_single)
    print "window singletons: " + str(count_windows(window_filt_single))
    print "window dir", get_directional_info(window_filt_single)

#    plot_depth_hist(window_filt_single,"all_final_singletons_depth.png")
#    plot_freq_hist(window_filt_single,"all_final_singletons_freq.png")
 

    window_filt_single_low = get_low_freq(window_filt_single)
    print "window filter singletons low ", count_windows(window_filt_single_low)
    print "window singletons low: " + str(count_windows(window_filt_single_low))
    print "window dir low", get_directional_info(window_filt_single_low)

    plot_depth_hist(window_filt_single_low,"final_singletons_depth.png")
    plot_freq_hist(window_filt_single_low,"final_singletons_freq.png")
 
    window_filt_single_homz = get_high_freq(window_filt_single)
    print "window filter singletons homz ", count_windows(window_filt_single_homz)
    print "window singletons homz: " + str(count_windows(window_filt_single_homz))
    print "window dir homz", get_directional_info(window_filt_single_homz)

    plot_depth_hist(window_filt_single_homz,"final_Hom_singletons_depth.png")
    plot_freq_hist(window_filt_single_homz,"final_Hom_singletons_freq.png")
 
    write_TE_file(single_sites,"temp_unfiltered_singletons.txt")
    write_TE_file(window_filt_single,"temp_singletons.txt")
    write_TE_file(window_filt_single_low,"templow_singletons.txt")
    write_TE_file(window_filt_single_homz,"temp_homz_singletons.txt")

   # print(get_family_counts(window_filt_single_low))

    sys.exit()

    global_sites = filter_ranges(global_sites)
    filter_singletons = filter_ranges(filter_singletons)
    window_filt_single = filter_ranges(window_filt_single)
    window_filt_single_low = filter_ranges(window_filt_single_low)
    window_filt_single_homz = filter_ranges(window_filt_single_homz)
    print "AFTER"
    print "global",count_windows(global_sites)
    print "filter1",count_windows(filter_singletons)
    print "filt2",count_windows(window_filt_single)
    print "filt2_low",count_windows(window_filt_single_low)
    print "filt2_homz",count_windows(window_filt_single_homz)
    



    #Plot range sizes
    TE_pos_ranges = filter_ranges(TE_pos)
    plot_freq_hist(TE_pos_ranges,"ra_all_TE_ranges.png")
    plot_freq_hist(global_sites,"ra_global_Freq.png")
    plot_freq_hist(window_filt_single,"ra_wind_freq.png")
    plot_freq_hist(global_sites,"ra_global_Freq.png")
    plot_freq_hist(window_filt_single,"ra_wind_freq.png")



    ref_TE_loc = "refs_TE.txt"


#    high_freq = get_high_freq(window_filt_single)
    write_TE_file(window_filt_single_homz,"high_freq_singletons.txt")
    write_TE_file(window_filt_single,"all_singletons.txt")
#    write_TE_file(global_sites,"global_sites.txt")
    sys.exit()
#    plot_depth_hist(single_sites,"unfilt_singleton_TE_depth.png")
#    plot_depth_hist(window_filt_single,"singletons_TE_depth.png") 
#    plot_depth_hist(global_sites,"global_TE_depth.png")


#    plot_freq_hist(single_sites,"unfilt_singleton_freq.png")
    plot_freq_hist(window_filt_single,"no_max_singletons_freq.png") 
#    plot_freq_hist(global_sites,"global_freq.png")

    
    
#    write_singleton_file(window_filt_single,"window_filt_singletons.txt")

#    plot_TEDepth_hist(window_filt_single,"singleton_depths.png")
#    plot_TEDepth_hist(global_sites,"global_depths.png")

 
    sys.exit() 

   
    #Directional support data
    print "all"
    print get_directional_info(TE_pos)

    print "global"
    print get_directional_info(global_sites)

    print "singleton"
    print_directional_info(single_sites)

    print "filt_singleton"
    print_directional_info(window_filt_single)



    sys.exit() 
    #Plot frequency of TE by distance to window edge 







