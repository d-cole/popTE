"""
TE_set.py
Represents a set of TE insertion sites found by popoolationTE.
"""
from polyLine import polyLine

class window_exists_error(Exception):
    """
    """
    def __init__(self,err_str):
        self.err_str = err_str
    def __str__(self):
        return repr(self.err_str)


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
        self._window_size = window_size
        self.sites = {} 


    def add_site(self,site):
        """
        Adds given set to instance of TE_set
        """
        pos_key = int(float(site.pos)) / self._window_size

        self.sites[site.chrom] = self.sites.get(site.chrom, {})
        self.sites[site.chrom][pos_key] = self.sites[site.chrom].get(pos_key,[])
        self.sites[site.chrom][pos_key].append(site)
        return


    def add_window(self,window_contents,chrom,pos_idx):
        """
        Adds given TE window to instance of TE_set
        """
        self.sites[chrom] = self.sites.get(chrom, {})

        if self.sites[chrom].get(pos_idx,None) != None:
            raise window_exists_error("Error adding TE window, window of specified position already exists.")

        self.sites[chrom][pos_idx] = window_contents
        return


    def load_sites_txt(self,file_locations):
        """
        Given list of file_locations load all specified files adding all sites to TE_set
        """
        for file_name in file_locations:
            with open(file_name) as poly_file:
                for line in poly_file:
                    sample_name = file_name.split("/")[-1] #Should be last item in path 'CC_A etc.'
                    pLine = polyLine(line, sample_name)
    
                    self.add_site(pLine) 

        #Remove unwanted chromosomes 
        bad_chroms = ["pseudo0","mitochondrion","chloroplast"]
        for key in bad_chroms:
            self.sites.pop(key,None)

        return

#    def write_txt(self,output_loc):
#        """
#        Write out TE set in modified popTE polymorphism output
#        Modified polyLines include sample name at column 0
#        """
#        output = open(output_loc, "w")
#        iter_all_sites(__write_site__,output)
#         

 
    def get_site_filtered_set(self,filter_method,*args):
        """
        Return a new TE set based on a filter applied on TE sites
        """
        new_set = TE_set(self._window_size)
        for chrom in self.sites.keys():
            for pos_range in self.sites[chrom].keys():
                for site in self.sites[chrom][pos_range]:
                    if filter_method(site,*args):
                        new_set.add_site(site)
        return new_set
   
 
    def get_window_filtered_set(self,filter_method,*args):
        """
        Return a new TE set based on a filter applied on TE windows 
        """
        new_set = TE_set(self._window_size)
        for chrom in self.sites.keys():
            for pos_range in self.sites[chrom].keys():
                if filter_method(self.sites[chrom][pos_range],*args):
                        new_set.add_window(self.sites[chrom][pos_range],chrom,pos_range)

        return new_set
        

    def get_num_sites(self):
        """
        """
        pass         


    def get_num_windows(self):
        """
        Returns the total number of windows in TE set
        """
        num_windows = 0
        for chrom in self.sites.keys():
            num_windows += len(self.sites[chrom])
        return num_windows


    def get_num_sites(self):
        """
        Returns the total number of sites in TE set
        """
        num_sites = 0
        for chrom in self.sites.keys():
            for pos_range in self.sites[chrom].keys():
                num_sites += len(self.sites[chrom][pos_range])

        return num_sites


    def _get_depth(self,site):
        """
        """
        return site.depth


    def get_all_depth(self):
        """
        """
        return self.get_sites_info(self._get_depth)


    def _get_freq(self,site):
        """
        Returns frequency of given site
        """
        return site.freq


    def get_all_freq(self):
        """
        """
        return self.get_sites_info(self._get_freq) 

 
    def get_sites_info(self, info_method):
        """
        """
        info = []
        for chrom in self.sites.keys():
            for pos_range in self.sites[chrom].keys():
                for site in self.sites[chrom][pos_range]:
                    info.append(info_method(site))                   
        return info


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
        num_sites = self.get_num_sites()
        num_windows = self.get_num_windows()

        out_str = "Number of sites: " + str(num_sites) + '\n'
        out_str = out_str + "Number of windows: " + str(num_windows) + '\n'
           
        return out_str 


    def subtract(self,other):
        """
        Returns new TE_Set = self - other
        
            - Use for removing sites within TE masked range?
            - Can TE masked ranges be translated to TE_set?
        """
        pass











