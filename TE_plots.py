"""
TE_plots.py

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
from polyFilter import *


def make_plot(fig_num,x_vals,y_vals,x_label,y_label,filename,in_label=None):
    """ 
    Plot the given data with axis labels. Legend/data labels optional.
    Saves the image as the file specified by filename
    *Uses subplot
    *From A1 code
    """
    f,axarr = plt.subplots(1,1,sharex=False, sharey=False)
    axarr.plot(x_vals,y_vals,label=in_label)
    axarr.set_xlabel(x_label)
    axarr.set_ylabel(y_label)
    plt.savefig(filename,bbox_inches='tight')


def plot_freq_window_pos(TE_pos,n,file_prefix):
    """
    Plots the frequency by window position.
        Are hard windows causing false positives? 
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
    axarr.hist(freq) 
    axarr.set_xlabel("pos_diff ")
    axarr.set_ylabel("frequency")
    plt.savefig(file_prefix + "freq_dist_window_edge.png",bbox_inches='tight')
    return


def plot_freq_TE_depth(TE_pos,filename):
    """
    Plot TE frequency by depth
    """
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
    """
    Plots histogram of x_vals, save plot as filename.
    """
    f,axarr = plt.subplots(1,1,sharex=True,sharey=True)
    axarr.hist(x_vals)
    axarr.set_ylabel("Frequency")
    axarr.set_xlabel(x_label)
    plt.savefig(filename,bbox_inches='tight')


def plot_scatter(x_vals,x_label,y_vals,y_label,filename):
    """
    Plot y_vals vs. x_vals, save plot as filename.
    """
    f,axarr = plt.subplots(1,1,sharex=True,sharey=True)
    axarr.scatter(x_vals,y_vals)
    axarr.set_ylabel(y_label)
    axarr.set_xlabel(x_label)
    plt.savefig(filename,bbox_inches='tight')


def plot_freq_hist(TE_dict,filename):
    """
    Plots histogram of TE frequencies, saves as filename.png
        Hides spines and ticks
    """

    freq = []
    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom].keys():
            freq_sum = 0
            site_num = 0
            for site in TE_dict[chrom][pos_range]:
                freq_sum += site.freq
                site_num += 1
            freq.append(float(freq_sum)/float(site_num))

    f,axarr = plt.subplots(1,1,sharex=True,sharey=True)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    axarr.hist(freq)
    axarr.spines['right'].set_visible(False)
    axarr.spines['top'].set_visible(False) 
    axarr.xaxis.set_ticks_position('none')
    axarr.yaxis.set_ticks_position('none')
    axarr.set_ylabel("Number of TEs",fontsize=20)
    axarr.set_xlabel("Frequency of TE insertion",fontsize=20)
    plt.savefig(filename,bbox_inches='tight',dpi=500)


def plot_depth_hist(TE_dict,filename):
    """
    Plot histogram of TE depths, save .png as given filename.
    Hides spines and ticks
    """
    freq = []
    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom].keys():
            freq_sum = 0
            site_num = 0
            for site in TE_dict[chrom][pos_range]:
                if site.total_TE_depth < 30:
                    freq_sum += site.total_TE_depth
                    site_num += 1
            if site_num != 0: 
                freq.append(float(freq_sum)/float(site_num))
    
    f,axarr = plt.subplots(1,1,sharex=True,sharey=True)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    axarr.hist(freq)
    axarr.spines['right'].set_visible(False)
    axarr.spines['top'].set_visible(False) 
    axarr.xaxis.set_ticks_position('none')
    axarr.yaxis.set_ticks_position('none')
    axarr.set_ylabel("Number of TEs",fontsize=20)
    axarr.set_xlabel("Depth of TE insertion",fontsize=20)
    plt.savefig(filename,bbox_inches='tight',dpi=500)


def plot_TEDepth_hist(TE_dict,filename):
    """
    Plots a histogram of the TE depth for all sites in TE_dict.
    """
    TE_depths = get_all_TE_depth(TE_dict)
   
    #Restrict depth below 50 to omit weird 500depth global TEs from skewing scale
    TE_depths = [x for x in TE_depths if x <= 20]
    plot_hist(TE_depths,"TE depths",filename)


def plot_count_by_samples_in(TE_dict,target_genotype,filename):
    """
    Plots number of TEs found by the number of samples that share the TE.
    """
    num_sample_TE = []
    num_sample_count = []

    for num_samples in range(1,len(target_genotype) + 1):
        TE_in_num_samples = count_windows(get_ranges_by_sample_count(TE_dict,num_samples))
    
        print TE_in_num_samples
        print num_samples
        print
        num_sample_TE.append(TE_in_num_samples)
        num_sample_count.append(num_samples)

    make_plot(0,num_sample_count,num_sample_TE,"Number of samples","Number of TEs",filename)
 

def plot_TEDepth_by_samples_in(TE_dict):
    """
    Plots the depth of insertions by the number of samples with the TE present.
    """
    depths = []
    num_sample_TE = []

    for num_samples in range(1,len(target_genotype) + 1):
        TE_in_num_samples = get_ranges_by_sample_count(TE_dict,num_samples)
        TE_depths = get_all_TE_depth(TE_in_num_samples)
   
        #Plots all points     
        depths = depths + TE_depths
        num_sample_TE = num_sample_TE + [num_samples]*len(TE_depths)

        #Plot averages
        #avgDepth = float(sum(TE_depths))/float(len(TE_depths))
        #depths.append(avgDepth)
        #num_sample_TE.append(num_samples)

#    plot_scatter(num_sample_TE,"Num samples TE in", depths, "Total TE depth","TEDepth_by_samples.png") 

    plot_hist(depths,"TE depths","TEDepth_by_sample_count:.png")
    return







