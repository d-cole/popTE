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

#!!! - SAMPLES G and M


CC = ["CC_A","CC_B","CC_C","CC_D","CC_E","CC_F","CC_G","CC_H","CC_I","CC_J","CC_K","CC_L","CC_M","CC_N","CC_O","CC_P"]
GP = ["GP2-3_A","GP2-3_B","GP2-3_C","GP2-3_D","GP2-3_E","GP2-3_F","GP2-3_G","GP2-3_H","GP2-3_I","GP2-3_J","GP2-3_K","GP2-3_L","GP2-3_M","GP2-3_N","GP2-3_O","GP2-3_P"]

#Excludes samples G and M from CC genotype
CC_mod = ["CC_A","CC_B","CC_C","CC_D","CC_E","CC_F","CC_H","CC_I","CC_J","CC_K","CC_L","CC_N","CC_O","CC_P"]

both_genotypes = CC_mod + GP
target_genotype = both_genotypes

#Remember calculated Pvals to speed finding Pvals
Pval_dict = {}


def num_sample_in_window(window):
    """
    window -> list of pLine objects in TE_pos[chrom][window_idx] -->  window
    """
    samples_present = set() 
    for pline in window:
        samples_present.add(pline.sample) 
    
    return len(samples_present)


def range_filters(singleton_site):
    """
    """
    pass        



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


def passBinom(TE_site):
    """ 
    
    """
    altReads = TE_site.TE_reads_forward + TE_site.TE_reads_reverse
    refReads = TE_site.nonTE_reads_forward + TE_site.nonTE_reads_reverse

    binomResult = "TRUE" in str((subprocess.Popen("~/Documents/spirodela/r/binom.r " + str(altReads)\
             +" "+str(refReads),shell=True,stdout=subprocess.PIPE)).communicate()[0])
    #print binomResult
    return (binomResult and (altReads + refReads) > 3)


def getPval(TE_site):
    """
    """
    altReads = TE_site.TE_reads_forward + TE_site.TE_reads_reverse
    refReads = TE_site.nonTE_reads_forward + TE_site.nonTE_reads_reverse

    binomResult = str((subprocess.Popen("~/Documents/spirodela/r/calcP.r " + str(altReads)\
             +" "+str(refReads),shell=True,stdout=subprocess.PIPE)).communicate()[0])
    binomResult = re.findall(r'"([^"]*)"',binomResult)
    binomResult = float(binomResult[0])
    return binomResult


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


def window_filter(window):
    """
    Filters:
        - Must have same family?
        - Must have forward and reverse insertions present
        - Frequency around 0.5?
        - Atleast converage of ?
    """
    min_avg_TE_coverage = 3

    forward = False
    Reverse = False
    Coverage = 0
    family = None

    # At least one forward and one reverse of the same family
    pass_dir_test = True
    #family_directions={"F":[],"R":[],"FR":[]}

    #for pLine in window:
    #    family_directions[pLine.direction].append(pLine.family)

    ##Temp
    #print family_directions    


    #if (len(family_directions["FR"]) != 0) or (len(family_directions["F"] \
    #    + family_directions["R"]) != len(set(family_directions["F"] + family_directions["R"]))):
    #    
    #    pass_dir_test = True
    #           
  
    #IF ONE site in the window passes the binomial test --> window passes 
    binom_pass = False
    for site in window:
        if passBinom(site):
            binom_pass = True
    
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
    return (binom_pass ) 


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

    for file_name in file_name_list:
            with open(file_name) as file:
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
    f,axarr = plt.subplots(1,1,sharex=True, sharey=True)
    axarr.plot(x_vals,y_vals,label=in_label)
    axarr.set_xlabel(x_label)
    axarr.set_ylabel(y_label)
    plt.savefig(filename,bbox_inches='tight')


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

def get_avg_depths(TE_pos):
    """
    Get array of all 
    """
    pass


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


if __name__ == "__main__":
    #!Base pair window size
    n = 1000

    TE_pos = get_TE_dict(target_genotype,n)

    single_sites = get_singleton_ranges(TE_pos)
    global_sites = get_global_ranges(TE_pos)

    print "singleton ranges: " + str(count_windows(single_sites))
    print "global ranges: " + str(count_windows(global_sites))

    filter_singletons = get_filtered_ranges(single_sites)
    print "filtered singletons: " + str(count_windows(filter_singletons))

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
#    all_freq = get_all_freq(get_ranges_by_sample_count(TE_pos,len(target_genotype)))
#    plot_hist(all_freq,"Estimated TE frequency at site","freq_global_te_sites.png")
#    
#    #Get singleton and filtered singleton frequencies & PLOT
#    singleton_freq = get_all_freq(single_sites)
#    filt_singleton_freq = get_all_freq(filter_singletons)
#    plot_hist(singleton_freq,"Singleton EST TE freq","singleton_freq.png")
#    plot_hist(filt_singleton_freq,"Filtered (Depth 3 filter) Singleton Est TE freq","filt_singleton_Freq.png")    


    #Get Pval distributions
    all_Pvals = getAllPvals(TE_pos)
    filtered_Pvals = getAllPvals(filter_singletons)
    singleton_Pvals = getAllPvals(single_sites)
    global_Pvals = getAllPvals(global_sites)
    
    plot_hist(all_Pvals,"Pvalue of TEread vs. nonTE reads","1kbp_all_Pval.png")
    plot_hist(singleton_Pvals,"Pvalue of TEread vs. nonTE reads","1kbp_t2_singleton_Pvals.png")
    plot_hist(global_Pvals,"Pvalue of TEreads vs. nonTE reads","1kbp_global_Pvals.png")
    plot_hist(filtered_Pvals,"Pvalue of TEreads vs. nonTE reads","1kbp_filt_singletons_Pvals.png")




