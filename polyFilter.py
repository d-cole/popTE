import sys
from polyLine import polyLine


#!!! - SAMPLES G and M


CC = ["CC_A","CC_B","CC_C","CC_D","CC_E","CC_F","CC_G","CC_H","CC_I","CC_J","CC_K","CC_L","CC_M","CC_N","CC_O","CC_P"]
GP = ["GP2-3_A","GP2-3_B","GP2-3_C","GP2-3_D","GP2-3_E","GP2-3_F","GP2-3_G","GP2-3_H","GP2-3_I","GP2-3_J","GP2-3_K","GP2-3_L","GP2-3_M","GP2-3_N","GP2-3_O","GP2-3_P"]

target_genotype = CC


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


def filter_sites(singleton_TEs):
    """
    """
    pass
    
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



if __name__ == "__main__":
    #!Base pair window size
    n = 200


    #Initialize chromosomes in TE_pos dict
    TE_pos = {"mitochondrion":{},"chloroplast":{}} 
    for i in range(0,33):
        chrom_str = "pseudo" + str(i)
        TE_pos[chrom_str] = TE_pos.get(chrom_str,{})


    #Load all sites from all files into TE_pos
    for file_name in target_genotype:
        
        #Load in all sites from one file
        with open(file_name) as file:
            for line in file:
                #Init polyLine object representing one line in the polymorphism file
                pLine = polyLine(line,file_name)
        #        print  "pos: " + str(int(float(pLine.pos)))

                #Get the postion key
                pos_key = int(float(pLine.pos)) / n 

        #        print pos_key
                TE_pos[pLine.chrom][pos_key] = TE_pos[pLine.chrom].get(pos_key,[])
                TE_pos[pLine.chrom][pos_key].append(pLine)
    


    #Find singleton cases
    i = 0
    for chrom in TE_pos.keys():
        for pos_range in TE_pos[chrom].keys():
            if num_sample_in_window(TE_pos[chrom][pos_range]) == 1:
                i+=1
    print "singleton ranges: " + str(i)


    single = get_singleton_ranges(TE_pos)

    write_singleton_file(single,target_genotype[0][0:2])


