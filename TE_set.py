"""
TE_set.py
Represents a set of TE insertion sites found by popoolationTE.
"""

class TE_set:
    """
   
    Attributes:
        - # of TE windows
        - # of TE sites
    """

    def __init__(self,window_size):
        """
        Initialize empty TE_set
        """
        self.num_windows = 0
        self.num_sites = 0
        self.window_size = window_size

        self.sites = {} 


    def load_insert_sites_txt(self,file_locations):
        """
        Given list of file_locations load all specified files adding all sites to TE_set
        """
        for file_name in file_locations:
            with open(file_name) as poly_file:
                for line in poly_file:
                    sample_name = file_name.split("/")[-1] #Should be last item in path 'CC_A etc.'
                    pLine = polyLine(line, sample_name)
                    pos_key = int(float(pLine.pos)) / self.window_size  
    
                    self.sites[pLine.chrom] = self.sites.get(chrom, {})
                    self.sites[pos_key] = self.sites.get(pos_key,[])
                    self.sites[pLine.chrom][pos_key].append(pLine)

        #Remove unwanted chromosomes 
        bad_chroms = ["pseudo0","mitochondrion","chloroplast"]
        for key in bad_chroms:
            self.sites.pop(key,None)

        return
          
    def iter_all_sites(self, per_site_function, *args):
        """
        Apply per_site_function on all sites in self.sites
        """
        for chrom in self.sites.keys():
            for pos_range in self.sites[chrom].keys():
                for pLine in self.sites[chrom][pos_range]:
                    per_site_function(pLine,*args)


    def __write_site__(self,site,*args):
        """
        Given site and output stream in *args writes site to output stream
        """
        out_file = args[0]
        out_file.write(site.sample + "\t" + site.raw_line)


    def write_txt(self,output_loc):
        """
        Write out TE set in modified popTE polymorphism output
        Modified polyLines include sample name at column 0
        """
        output = open(output_loc, "w")
        iter_all_sites(__write_site__,output)
         





    def get_set(self):
        """
        """
        pass
       
 
    def get_filtered_set(self,filter_method):
        """
        """
        pass
    
    def get_sites_info(self, info_method):
        """
        """
        pass


    def remove_site(self,chrom,pos):
        """
        Remove specified TE site from set.
        Remove TE window if now empty
        """
        pass

    
    def remove_TE_window(self,chrom,pos_range):
        """
        Remove specified TE window
        """ 
        pass


    def get_family_info(self):
        """
        Returns TE family information
        """
        pass


    def get_dir_info(self):
        """
        Returns direction of support counts
        """
        pass


    def __repr__(self):
        """
        Modified repr to include information about TE set

        Number of windows
        Number of sites
        Number of unique sample types
        """
        pass


    def __equals__(self,other):
        """
        Determine if two TE sets are equal
        - is this needed?
        """
        pass


    def subtract(self,other):
        """
        Returns new TE_Set = self - other
        
            - Use for removing sites within TE masked range?
            - Can TE masked ranges be translated to TE_set?
        """
        pass











