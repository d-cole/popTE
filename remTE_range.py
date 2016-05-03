"""
remTE_range.py
"""

import sys
#
CHROM,START,STOP = 0,1,2

def isDataLine(line):
    """
    Determines in line contains site data
    """
    if len(line) > 1:
        return line[0] != "#"
    return False


def loadTEranges(TE_file_loc):
    """
    Load TE ranges from a text file of the format CHROM START STOP
    """
    TE_ranges = {}

    with open(TE_file_loc) as TE_file:
        for line in TE_file:
            line_col = str.split(line)

            TE_ranges[line_col[CHROM]] = TE_ranges.get(line_col[CHROM],[])
            TE_ranges[line_col[CHROM]].append((line_col[START],line_col[STOP]))


    TE_file.close()

    return TE_ranges

def get_valid_ranges(TE_dict,TE_ranges):
    """
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
    """
    n = 500
    chrom  = site.chrom
    pos = site.pos

    print "here"
    if site.chrom == "pseudo13" and site.sample == "CC_A":
        right_pos = pos + n
        left_pos = pos - n
        print "right", right_pos
        print "left", left_pos

        near = False
        #TE range = left_pos --> high_pos

        for point in range(left_pos,right_pos):
            print "point",point
            for (low,high) in TE_ranges[chrom]:
                if point <= high and point >= low:
                    print "inside"
                    near = True
    
    
    right_pos = pos + n
    left_pos = pos - n
    near = False
    #TE range = left_pos --> high_pos
    for point in range(left_pos,right_pos):
        for (low,high) in TE_ranges[chrom]:
            if point <= high and point >= low:
                near = True

    return not near





def get_valid_ranges(TE_dict,TE_ranges):
    """
    """
    valid_TE = {}
    for chrom in TE_dict.keys():
        for pos_range in TE_dict[chrom].keys():
            for site in TE_dict[chrom][pos_range]:
                if validRange(site,TE_ranges):
                    valid_TE[chrom] = valid_TE.get(chrom,{})
                    valid_TE[chrom][pos_range] = valid_TE[chrom].get(pos_range,[])
                    valid_TE[chrom][pos_range].append(site) 
    return valid_TE


def validRange(site,TE_ranges):
    """
    Determine if the given site is withen any TE range
    """
    chrom  = site.chrom
    pos = site.pos


    if any(float(low) <= float(pos) <= float(high) for (low,high) in TE_ranges[chrom]):
        return False

    return True


if __name__ == "__main__":
    TE_file_loc,vcf_file_loc = sys.argv[1],sys.argv[2]
    loadTEranges(TE_file_loc)

    trimmed_vcf = open(vcf_file_loc[:-4] + "_TEremoved.vcf", "w")
    removed_sites = open(vcf_file_loc[:-4] + "_TE_Sites.vcf","w")

    with open(vcf_file_loc) as vcf_file:
        for line in vcf_file:
            if isDataLine(line):
                if validRange(line):
                    trimmed_vcf.write(line)
                else:
                    removed_sites.write(line)
            else:
                #Write header info 
                trimmed_vcf.write(line)

    removed_sites.close()
    vcf_file.close()
    trimmed_vcf.close()
