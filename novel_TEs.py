from TE_set import TE_set
from polyLine import polyLine
import TE_filters as TE_filters

window_size = 1000

CC = ["CC_A","CC_B","CC_C","CC_D","CC_E","CC_F","CC_G","CC_H","CC_I","CC_J","CC_K","CC_L","CC_M","CC_N","CC_O","CC_P"]
GP = ["GP2-3_A","GP2-3_B","GP2-3_C","GP2-3_D","GP2-3_E","GP2-3_F","GP2-3_G","GP2-3_H","GP2-3_I","GP2-3_J","GP2-3_K","GP2-3_L","GP2-3_M","GP2-3_N","GP2-3_O","GP2-3_P"]

#Excludes samples G and M from CC genotype
CC_mod = ["CC_A","CC_B","CC_C","CC_D","CC_E","CC_F","CC_H","CC_I","CC_J","CC_K","CC_L","CC_N","CC_O","CC_P"]
wo_AP_CC = ["CC_B","CC_C","CC_D","CC_E","CC_F","CC_H","CC_I","CC_J","CC_K","CC_L","CC_N","CC_O"]
wo_AP_GP = ["GP2-3_B","GP2-3_C","GP2-3_D","GP2-3_E","GP2-3_F","GP2-3_G","GP2-3_H","GP2-3_I","GP2-3_J","GP2-3_K","GP2-3_L","GP2-3_M","GP2-3_N","GP2-3_O"]

both_genotypes = CC_mod + GP
wo_AP_both = wo_AP_CC + wo_AP_GP
target_genotype = both_genotypes
poly_dir = "./poly_files/"

target_genotype_paths = []
for sample in target_genotype:
    target_genotype_paths.append(poly_dir + sample)    


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




if __name__ == "__main__":
    #Initialize empty TE_set
    all_TE = TE_set(window_size)
        
    #Load sites from .txt files into TE_set
    all_TE.load_sites_txt(target_genotype_paths)
    print all_TE
    
    #Filter based on TE-absence range
    min_TE_abs_range = 73
    abs_filt = all_TE.get_site_filtered_set(TE_filters.TE_absence_range_filter,min_TE_abs_range)
    print abs_filt
    
    #Filter for singletons and global 
    singleton_TE = abs_filt.get_window_filtered_set(TE_filters.filter_window_by_sample_count,1)
    global_TE = abs_filt.get_window_filtered_set(TE_filters.filter_window_by_sample_count,len(target_genotype))
    print singleton_TE
    print global_TE
    
    #Filter singletons based on quality
    min_TE_depth = 3
    max_TE_depth = 30
    min_binom_prob = 0.02
    filt_singleton_TE = singleton_TE.get_window_filtered_set(TE_filters.filter_singleton_window, min_TE_depth, max_TE_depth, min_binom_prob)
    print filt_singleton_TE
    
    #Filter out false-positive singletons based on hard cutoff window

#    edge_filt_singleton_TE = filt_singleton_TE.get_site_filtered_set(TE_filters.singleton_edge_filter, all_TE, window_size) --> 28 sites
    edge_filt_singleton_TE = filt_singleton_TE.get_site_filtered_set(TE_filters.singleton_edge_filter, abs_filt, window_size) #--> 29 sites
    print edge_filt_singleton_TE
    
    #Filter based on distance to masked sequence
    masked_ranges = loadTEranges("masked_ranges/no_small_gtf_parsed.txt")
    min_masked_dist = 75
    masked_edge_filt_singleton_TE = edge_filt_singleton_TE.get_site_filtered_set(TE_filters.dist_masked_filter,masked_ranges, min_masked_dist)
    print masked_edge_filt_singleton_TE
    
    
    
