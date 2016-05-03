"""
rem_ref_TE.py
    - Given TE polymorphism file writes new polymorphism file sites with TEs within reference TE ranges removed.
"""
import sys
from polyLine import polyLine
#TE = {Pseudo0:[(x,y),(z,l)...],Pseudo1:[(m,x)...]...}
TE_ranges = {}
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
    with open(TE_file_loc) as TE_file:
        for line in TE_file:
            line_col = str.split(line)
            TE_ranges.setdefault(line_col[CHROM],[]).append((line_col[START],line_col[STOP]))

    TE_file.close()
    return

def validRange(line):
    """
    Determine if the given site is withen any TE range
    """
    chrom = line.chrom
    pos = line.pos

    if any( (float(low) - 500) <= float(pos) <= (float(high)+500) for (low,high) in TE_ranges[chrom]):
        return False

    return True


if __name__ == "__main__":
    TE_file_loc,vcf_file_loc = sys.argv[1],sys.argv[2]
    loadTEranges(TE_file_loc)

    trimmed_vcf = open("ref_range_filt.txt", "w")

    with open(vcf_file_loc) as vcf_file:
        for line in vcf_file:
                if validRange(polyLine(line,"CC_A")):
                    trimmed_vcf.write(line)

    trimmed_vcf.close()
