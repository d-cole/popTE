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

#!!! - SAMPLES G and M


CC = ["CC_A","CC_B","CC_C","CC_D","CC_E","CC_F","CC_G","CC_H","CC_I","CC_J","CC_K","CC_L","CC_M","CC_N","CC_O","CC_P"]
GP = ["GP2-3_A","GP2-3_B","GP2-3_C","GP2-3_D","GP2-3_E","GP2-3_F","GP2-3_G","GP2-3_H","GP2-3_I","GP2-3_J","GP2-3_K","GP2-3_L","GP2-3_M","GP2-3_N","GP2-3_O","GP2-3_P"]

#Excludes samples G and M from CC genotype
CC_mod = ["CC_A","CC_B","CC_C","CC_D","CC_E","CC_F","CC_H","CC_I","CC_J","CC_K","CC_L","CC_N","CC_O","CC_P"]

both_genotypes = CC_mod + GP
target_genotype = both_genotypes

#Remember calculated Pvals to speed finding Pvals
Pval_dict = {}
n = 1000


def get_samples_in_range(TE_pos,target_site,range_size):
    """
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
    if window_pos > 500:
        #Want to check range pos:pos+500
        pos_right_border = pos + range_size
        print "right",pos_idx
        pos_idx = int(pos_right_border/n)
        window = TE_pos[chrom].get(pos_idx,None)

        if window != None:
            for site in window:
               if site.pos < pos_right_border:
                    new_window_sites.append(site)
        
    #Only check the left window  
    if window_pos < 500:
        pos_left_border = pos - range_size
        pos_idx = int(pos_left_border/n)
        print "left",pos_idx
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
            
                print "single range: ", site_range 
                #Get samples with TEs present in the new window 500 bp on either side of TE
                new_window_samples = get_samples_in_range(TE_pos,singleton_site,window_width)
                false_singleton = False 

                for TE in new_window_samples:
                    if TE.sample != singleton_site.sample:
                        false_singleton = True
                
                if false_singleton is not True:
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
    """
    TE_singletons = {}
    
    for chrom in TE_pos.keys():
        for pos_range in TE_pos[chrom].keys():
            if num_sample_in_window(TE_pos[chrom][pos_range]) == 1:
                TE_singletons[chrom] = TE_singletons.get(chrom,{})
                TE_singletons[chrom][pos_range] = TE_pos[chrom][pos_range]

    return TE_singletons


def passBinom(TE_site, min_prob):
    """ 
    
    """
    altReads = TE_site.TE_reads_forward + TE_site.TE_reads_reverse
    refReads = TE_site.nonTE_reads_forward + TE_site.nonTE_reads_reverse
    binomResult = (getPval(TE_site) > min_prob)
    return (binomResult and (altReads + refReads) > 3)


def passDepth(TE_site, min_depth):
    """
    """
    altReads = TE_site.TE_reads_forward + TE_site.TE_reads_reverse
    refReads = TE_site.nonTE_reads_forward + TE_site.nonTE_reads_reverse
    return (altReads + refReads) > min_depth


def passTEDepth(TE_site, min_depth, max_depth):
    """
    """
    altReads = TE_site.TE_reads_forward + TE_site.TE_reads_reverse
    return (altReads > min_depth) and (altReads < max_depth)


def getPval(TE_site):
    """
    """
    altReads = TE_site.TE_reads_forward + TE_site.TE_reads_reverse
    refReads = TE_site.nonTE_reads_forward + TE_site.nonTE_reads_reverse
    return pbinom(altReads,refReads,0.5)


def getAllPvals(TE_pos):
    """
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


def get_ranges_by_sample_range(TE_pos, sample_range):
    """
    Generates a dict of sites where a TE is present in sample_count or more samples
    """
    TEs_present = {}

    for chrom in TE_pos.keys():
            for pos_range in TE_pos[chrom].keys():
                if num_sample_in_window(TE_pos[chrom][pos_range]) in sample_range:
                    TEs_present[chrom] = TEs_present.get(chrom,{})
                    TEs_present[chrom][pos_range] = TE_pos[chrom][pos_range]

    return TEs_present


def window_filter(window):
    """
    Filters:
        - Must have same family?
        - Must have forward and reverse insertions present
        - Frequency around 0.5?
        - Atleast converage of ?
    """
    min_TE_depth = 3
    max_TE_depth = 100
    min_prob = 0.02

    forward = False
    Reverse = False
    Coverage = 0
    family = None

    # At least one forward and one reverse of the same family
    #pass_dir_test = True
    #family_directions={"F":[],"R":[],"FR":[]}
    #for pLine in window:
    #    family_directions[pLine.direction].append(pLine.family)
    #
    #if (len(family_directions["FR"]) != 0) or (len(family_directions["F"] \
    #    + family_directions["R"]) != len(set(family_directions["F"] + family_directions["R"]))):
    #    
    #    pass_dir_test = True
    #           
  
    #IF ONE site in the window passes the binomial test --> window passes 
    # IF One site (not necessarily the same) passes the binomial test --> window passes
    #       What would a site be like with the sites being passed from exclusive sites?
    #      --  Make it the same site to pass these tests?
    TE_depth_pass = False
    binom_pass = False
    for site in window:
        if passBinom(site, min_prob):
            binom_pass = True
        if passTEDepth(site,min_TE_depth,max_TE_depth):
            TE_depth_pass = True
 
    # Minimum average of coverage across sites in window
#    pass_coverage = False
#    num_sites = 0.0
#    coverage_total = 0.0 
#    for pLine in window:
#        num_sites += 1.0 
#        coverage_total += (pLine.TE_reads_forward + pLine.TE_reads_reverse)
#    
#    if coverage_total/num_sites > min_avg_TE_coverage:
#        pass_coverage = True
#
    return binom_pass  

def freq_by_TE_depth(TE_pos,depth):
    depth_freq = []
    for chrom in TE_pos.keys():
        for site_range in TE_pos[chrom].keys():
            for pos in TE_pos[chrom][site_range]:
                if passTEDepth(pos, depth):
                    depth_freq.append(pos.freq)
    return depth_freq
                     

def get_filtered_ranges(TE_pos):
    """
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
    """
    count = 0
    for chrom in TE_pos.keys():
        count += len(TE_pos[chrom].keys())    

    return count
   
 
def write_singleton_file(singletons,filename):
    """
    """
    file_out = open(filename,"w")
    for chrom in singletons.keys():
        for pos_range in singletons[chrom].keys():
            for pLine in singletons[chrom][pos_range]:
                out_line = pLine.sample + "\t" + pLine.raw_line
                file_out.write(out_line)    

    file_out.close()



def sites_per_bp_range():
    """
    """
    pass


def global_sites_per_n():
    """
    """
    pass


def init_TE_dict():
    """
    """
    TE_pos = {"mitochondrion":{},"chloroplast":{}} 
    for i in range(0,33):
        chrom_str = "pseudo" + str(i)
        TE_pos[chrom_str] = TE_pos.get(chrom_str,{})

    return TE_pos


def get_TE_dict(file_name_list,window_size):
    """
    """
    TE_pos = init_TE_dict()

    for file_name in  (file_name_list):
            with open("./poly_files/" + file_name) as file:
                for line in file:
                    #Init polyLine object representing one line in the polymorphism file
                    pLine = polyLine(line,file_name)
                    #Get the postion key
                    pos_key = int(float(pLine.pos)) / window_size
    
                    #print pos_key
                    TE_pos[pLine.chrom][pos_key] = TE_pos[pLine.chrom].get(pos_key,[])
                    TE_pos[pLine.chrom][pos_key].append(pLine)
    return TE_pos


def make_plot(fig_num,x_vals,y_vals,x_label,y_label,filename,in_label=None):
    """ 
    Plot the given data with axis labels. Legend/data labels optional.
    Saves the image as the file specified by filename

    *From A1 code
    """
    f,axarr = plt.subplots(1,1,sharex=False, sharey=False)
    axarr.plot(x_vals,y_vals,label=in_label)
    axarr.set_xlabel(x_label)
    axarr.set_ylabel(y_label)
    plt.savefig(filename,bbox_inches='tight')

def plot_freq_TE_depth(TE_pos,filename):
    freq = []
    TE_depth = []
    for chrom in TE_pos.keys():
        for pos_range in TE_pos[chrom].keys():
            for site in TE_pos[chrom][pos_range]:
                freq.append(site.freq)                
                TE_depth.append(site.total_TE_depth)

    f,axarr = plt.subplots(1,1,sharex=False, sharey=False)
    axarr.scatter(TE_depth,freq)
    axarr.set_xlabel("TE depth")
    axarr.set_ylabel("frequency")
    plt.savefig(filename,bbox_inches='tight')

    return

def plot_hist(x_vals,x_label,filename):
    f,axarr = plt.subplots(1,1,sharex=True,sharey=True)
    axarr.hist(x_vals)
    axarr.set_ylabel("Frequency")
    axarr.set_xlabel(x_label)
    plt.savefig(filename,bbox_inches='tight')


def plot_all_num_with_TE_by_window_size(file_name_list):
    for num_sample in range(1,len(file_name_list) + 1):
        plot_num_samples_with_TE_by_window_size(file_name_list,num_sample)
    return


def get_all_freq(TE_pos):
    freq = []
    for chrom in TE_pos.keys():
        for pos in TE_pos[chrom].keys():
            for site in TE_pos[chrom][pos]:
                freq.append(site.freq)
    return freq

     
def get_avg_freq(TE_pos):
    freq = get_all_freq(TE_pos)
    total = sum(freq) 
    return total/float(len(freq))


def plot_num_samples_with_TE_by_window_size(file_name_list,num_samples):
    """
    """
    counts = []
    window_sizes = []

    for window_size in range(100,5000,100):

        TE_pos = get_TE_dict(file_name_list,window_size)
        TE_pos_only_num_samples_sites = get_ranges_by_sample_count(TE_pos,num_samples)
        windows_with_num_samples = count_windows(TE_pos_only_num_samples_sites)
        
        window_sizes.append(window_size) 
        counts.append(windows_with_num_samples)

    filename = "num_sites_" + str(num_samples) + "_vs.Window_size.png"
    make_plot(0,window_sizes,counts,"Window size (bp)","Number of sites with TE present in " + str(num_samples) + " samples",filename)
    
    return


def get_excl_sites(CC_sites, GP_sites):
    """
    """
   #s.symmetric_difference(t) 
    mut_excl_sites = {}

    for chrom in CC_sites:
        for site_key in CC_sites[chrom].keys():
            if GP_sites[chrom].get(site_key, None) == None:
                mut_excl_sites[chrom] = mut_excl_sites.get(chrom, {})
                mut_excl_sites[chrom][site_key] = CC_sites[chrom][site_key]
            
    for chrom in GP_sites:
        for site_key in GP_sites[chrom].keys():
            if CC_sites[chrom].get(site_key, None) == None:
                mut_excl_sites[chrom] = mut_excl_sites.get(chrom, {})
                mut_excl_sites[chrom][site_key] = GP_sites[chrom][site_key]
    
    return mut_excl_sites


def plot_freq_window_pos(TE_pos,n,file_prefix):
    """
    """
    freq = []
    pos_diffs = []
    pos_freq_1 = []

    for chrom in TE_pos.keys():
        for pos_range in TE_pos[chrom].keys():
            for site in TE_pos[chrom][pos_range]:
                freq.append(site.freq)
                #Min diff from dist from pos to end of range or dist from pos to start of range
                
                pos_range_start = pos_range*n
                pos_range_end = pos_range*n + n
            
                # GETs difference for either end
                #pos_diff = min(abs((site.pos) - (pos_range_end)), abs((site.pos) - (pos_range_start)))
            
                pos_diff = site.pos - pos_range_start 
                pos_diffs.append(pos_diff)
            
                if site.freq == 1.0:
                    pos_freq_1.append(pos_diff)

    f,axarr = plt.subplots(1,1,sharex=False, sharey=False)
    #axarr.scatter(pos_diffs,freq)
    #axarr.hist(pos_diffs)
    axarr.hist(pos_freq_1)
    axarr.set_xlabel("pos_diff ")
    axarr.set_ylabel("frequency")
    plt.savefig(file_prefix + "freq_dist_window_edge.png",bbox_inches='tight')
    return

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
 
    single_sites = get_singleton_ranges(TE_pos)
    global_sites = get_global_ranges(TE_pos)

    print "singleton ranges: " + str(count_windows(single_sites))
    print "global ranges: " + str(count_windows(global_sites))

    filter_singletons = get_filtered_ranges(single_sites)
    print "filtered singletons: " + str(count_windows(filter_singletons))
    

    #New filter edge cases
    window_filt_single = filt_singleton_edge(TE_pos,filter_singletons)
    print "window filter singletons ", count_windows(window_filt_single)

    
    plot_freq_window_pos(window_filt_single,n,"window_")
#    plot_freq_TE_depth(window_filt_single,"window_filt_single.png")

#    sys.exit()
    #    PLOT FREQUENCY BY TE READ DEPTH
    #plot_freq_TE_depth(global_sites,"global.png")
    #plot_freq_TE_depth(single_sites,"single.png")
    #plot_freq_TE_depth(filter_singletons,"filt_single.png")
    #plot_freq_TE_depth(TE_pos,"all.png")

    #    PLOT FREQUENCY BY RELATIVE POS IN WINDOW 
    #       - ARE THE HOMOZYGOUS SINGLETONS ON WINDOW EDGES?  
#    plot_freq_window_pos(filter_singletons, n,"filt_single")
#    plot_freq_window_pos(TE_pos,n,"all")
#    plot_freq_window_pos(global_sites,n,"global_")
#    plot_freq_window_pos(single_sites,n,"single_")

    sys.exit()
    #write_singleton_file(filter_singletons,"filtered_singletons.txt")
#    plot_num_samples_with_TE_by_window_size(target_genotype,1)
#    plot_num_samples_with_TE_by_window_size(target_genotype,len(target_genotype))
#
#    plot_all_num_with_TE_by_window_size(target_genotype)


    #Plots number of TEs found by number of samples in
    num_samples = []
    counts = [] 
    avg_freq = []
    for i in range(1,len(target_genotype)+1):
        num_samples.append(i)
        target_ranges = get_ranges_by_sample_count(TE_pos,i)
        tar_avg_freq = get_avg_freq(target_ranges)
        avg_freq.append(tar_avg_freq)
        counts.append(count_windows(target_ranges))

    #Plot Number of samples with TE vs. number of TEs found
#    make_plot(0,num_samples,counts,"TE in all x samples","Number of TEs","GP_TEcount_v_num_samples.png")
#    make_plot(0,num_samples,avg_freq,"number of samples TEs in","Avg. TE frequency","TE_freq_by_in_all_samples.png")
#
#    #Get frequencies of TEs at sites
    all_freq = get_all_freq(get_ranges_by_sample_count(TE_pos,len(target_genotype)))
    plot_hist(all_freq,"Estimated TE frequency at site","freq_global_te_sites.png")
#    
#    #Get singleton and filtered singleton frequencies & PLOT
    singleton_freq = get_all_freq(single_sites)
    filt_singleton_freq = get_all_freq(filter_singletons)
    plot_hist(singleton_freq,"Singleton EST TE freq","singleton_freq.png")
    plot_hist(filt_singleton_freq,"Filtered (Depth 3 filter) Singleton Est TE freq","filt_singleton_Freq.png")    


    #Get Pval distributions
#    all_Pvals = getAllPvals(TE_pos)
#    filtered_Pvals = getAllPvals(filter_singletons)
#    singleton_Pvals = getAllPvals(single_sites)
#    global_Pvals = getAllPvals(global_sites)
#    
#    plot_hist(all_Pvals,"Pvalue of TEread vs. nonTE reads","1kbp_all_Pval.png")
#    plot_hist(singleton_Pvals,"Pvalue of TEread vs. nonTE reads","1kbp_t2_singleton_Pvals.png")
#    plot_hist(global_Pvals,"Pvalue of TEreads vs. nonTE reads","1kbp_global_Pvals.png")
#    plot_hist(filtered_Pvals,"Pvalue of TEreads vs. nonTE reads","1kbp_filt_singletons_Pvals.png")








