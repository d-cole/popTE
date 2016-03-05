#Col1: the reference sequence ID
#Col2: position in the reference sequence
#Col3: is the TE insertion supported by a forward (F), by a reverse (R) or by both (FR) insertions
#Col4: family of the TE insertion
#Col5: population frequency (1..fixed)
#Col6: order of the TE insertion
#Col7: ID if the TE insertion is present in the reference genome (e.g.: FlyBase ID)
#Col8: comment
#Col9: start of the range of the forward insertion
#Col10: end of the range of the forward insertion
#Col11: population frequency estimated by the forward insertion
#Col12: coverage of the forward insertion
#Col13: TE-presence reads of the forward insertion
#Col14: TE-absence reads of the reverse insertion
#Col15: is the range of the forward insertion overlapping with a forward-range of another TE insertion (0..no; 1..yes)
#Col16: start of the range of the reverse insertion
#Col17: end of the range of the reverse insertion
#Col18: coverage of the reverse insertion
#Col19: TE-presence reads of the reverse insertion
#Col20: TE-absence reads of the reverse insertion
#Col21: is the range of the reverse insertion overlapping with the reverese-range of another TE insertion (0..no; 1..yes) 

class polyLine:
    """
    A class for describing a 'site' in a popTE polymorphism output file
    """

#['mitochondrion', '1660', 'R', 'Unknown', '0.0127388535031847', 'Unknown', '-', '-', '-', '-', '-', '-', '-', '-', '-', '1760', '1834', '0.0127388535031847', '314', '4', '310', '0\n']

    def __init__(self,raw_line,sample_id):
        self.raw_line = raw_line
        self.sample = sample_id
        self.sline = raw_line.split("\t")
        self.chrom = self.sline[0]
        self.pos = self.sline[1]
        self.dirSupport = self.sline[2]
        self.family = self.sline[3]
        self.freq = self.sline[4]
        self.order = self.sline[5]
        self.id = self.sline[6]
        self.comment = self.sline[7]

        self.startRangeForward = self.sline[8]
        self.endRangeForward = self.sline[9]
        self.freqEstForward = self.sline[10]
        self.depthForward = self.sline[11]
        self.TEreadsFwd = self.sline[12]
        self.nonTE_reads_Fwd = self.sline[13]
        self.forward_range_overlap = self.sline[14]

        self.start_range_reverse = self.sline[15]
        self.end_range_reverse = self.sline[16]
        self.depth_reverse = self.sline[17]
        self.TE_reads_reverse = self.sline[18]
        self.nonTE_reads_reverse = self.sline[19]
        self.reverse_range_overlap = self.sline[20] 
        









