import sys
import matplotlib.pyplot as plt

class element_family:
    
    def __init__(self,head_type):
        self.head_type = head_type
        self.sequences = []

    def add_sequence(self,sequence):
        self.sequences.append(sequence) 

    def get_avg_size(self):
        total_len = 0.0
        num_seq = float(len(self.sequences))
        for sequence in self.sequences:
            total_len += len(sequence)
        
        return total_len/num_seq 

    def get_len_array(self):
        len_array = []
        for seq in self.sequences:
            len_array.append(len(seq)) 

        return len_array



class element:
    
    def __init__(self,header_line,content_lines):
        self.head = header_line
        self.head_type = header_line[self.head.find("#") + 1:self.head.find("(")]
        self.content=""

        for line in content_lines:
            self.content = self.content + (line)








with open("consensi.fa.classified") as file:

    elements = []

    #init params for each element
    header = ""
    content_lines = []

    for line in file:
        if line[0] == ">":
            elements.append(element(header,content_lines))
            
            #Reset for next element
            header = line
            content_lines = []

        else:
            content_lines.append(line)

    #Add the last sequence   
    elements.append(element(header,content_lines))    
 
    #Empty class first item
    elements.remove(elements[0])
            
# Check that elements are being read in correctly
#    for element in elements:
#        sys.stdout.write(element.head)
#        sys.stdout.write(element.content)


    #Ground elements into element_families
    element_groups = []
    for element in elements:
        found = False 
        #Check if element type is already in list, if so add sequence to it
        for element_group in element_groups:
            if element_group.head_type == element.head_type:
                found = True
                element_group.add_sequence(element.content)
        
        if not found:
            #Element type is not already in list
            new_ele_fam = element_family(element.head_type)
            new_ele_fam.add_sequence(element.content)
            element_groups.append(new_ele_fam)




    #Check that the families look right
    for element_family in element_groups:
        print element_family.head_type
        print element_family.get_avg_size()
        print element_family.get_len_array()
       
        plt.clf() 
        len_array = element_family.get_len_array()
        title = element_family.head_type + " n=" + str(len(len_array))
        plt.title(title)
        plt.hist(len_array)
        filename = str.rstrip(element_family.head_type)
        filename = filename.replace("/","_") + "_seq_len_dist.png"
        plt.savefig(filename)








