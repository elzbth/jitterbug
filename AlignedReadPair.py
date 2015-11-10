import sys



class AlignedReadPair:


    def __init__(self, read1, read2):
        """takes two pysam AlignedRead objects"""
        if read1.qname != read2.qname:
            print("error making aligned read pair: read names not the same")
            print(read1)
            print(read2)
            sys.exit(2)

        #self.read1 = read1
        #self.read2 = read2
        self.read1_str = str(read1)
        self.read2_str = str(read2)
        self.read1_is_TE = False
        self.read2_is_TE = False
        self.TE_annot_attr_list = []
        self.TE_map_gff_list = []
        self.interval_chr = ""
        self.interval_start = 0
        self.interval_end = 0
        self.interval_direction = None
        self.anchor_softclipped_pos = None

    def calculate_inside_interval(self, size, read1, read2):
        """calculate the interval of putative insertion sites corresponding to this read pair. \
        This interval is of length size and in the direction to which points the read that maps uniquely to a non-TE location
        a reasonable setting of size is the predicted INSIDE insert size (+  some sdev)"""
        ## INSIDE INTERVAL
        if not self.read1_is_TE:
            #this means that the uniquely mapping read to a non-TE location is read1
            if self.read1.is_reverse:
                self.interval_start = int(read1.pos) - size
                self.interval_end = int(read1.pos)
                self.interval_direction = "rev"
            else:
                self.interval_start = int(read1.aend)
                self.interval_end = int(read1.aend) + size
                self.interval_direction = "fwd"
        elif not self.read2_is_TE:
            #this means that the uniquely mapping read to a non-TE location is read2
            if self.read2.is_reverse:
                self.interval_start = int(read2.pos) - size
                self.interval_end = read2.pos
                self.interval_direction = "rev"
            else:
                self.interval_start = int(read2.aend)
                self.interval_end = int(read2.aend) + size
                self.interval_direction = "fwd"
        else:
            print("error in read pair object - neither read is a TE")


    def calculate_outside_interval(self, size, read1, read2):
        """calculate the interval of putative insertion sites corresponding to this read pair. \
        This interval is of length size and in the direction to which points the read that maps uniquely to a non-TE location, INCLUDING THE MAPPING LOCATION OF THE READS. This accounts for split reads, or truncated reads that would be shorter
        a reasonable setting of size is the predicted OUTSIDE insert(fragment) size (+  some sdev)"""
        ## OUTSIDE INTERVAL
        if not self.read1_is_TE:
            #this means that the uniquely mapping read to a non-TE location is read1
            if read1.is_reverse:
                self.interval_start = int(read1.aend) - size
                self.interval_end = int(read1.aend)
                self.interval_direction = "rev"
            else:
                self.interval_start = int(read1.pos)
                self.interval_end = int(read1.pos) + size
                self.interval_direction = "fwd"
        elif not self.read2_is_TE:
            #this means that the uniquely mapping read to a non-TE location is read2
            if read2.is_reverse:
                self.interval_start = int(read2.aend) - size
                self.interval_end = read2.aend
                self.interval_direction = "rev"
            else:
                self.interval_start = int(read2.pos)
                self.interval_end = int(read2.pos) + size
                self.interval_direction = "fwd"
        else:
            print("error in read pair object - neither read is a TE")
            sys.exit(2)

    def str(self):
#        #string = "read1: %s \n\
#                    read2: %s \n\
#                    read1_is_TE %s \n\
#                    read1_is_rev %s \n\
#                    read2_is_TE %s \n\
#                    read2_is_rev %s \n\
#                    TE_annot_list %s \n\
#                    interval_chr %s\n\
#                    interval_start %d \n\
#                    interval_end %d \n\
#                    interval_direction %s\n" % (self.read1, self.read2, self.read1_is_TE, self.read1.is_reverse, self.read2_is_TE, self.read2.is_reverse, "; ".join([attr['Name'] for attr in self.TE_annot_attr_list]), self.interval_chr, self.interval_start, self.interval_end, self.interval_direction)


        return string


    def str_int(self):
        return "[%s, %d, %d]" % (self.interval_chr, self.interval_start, self.interval_end)

    def str_TE_annot_list(self):
        try:
            return "; ".join([attr['Name'] for attr in self.TE_annot_attr_list])
        except KeyError:
           return "no_name"

    def TE_annot_tag_list(self, tag):
        name_list = []
        for attributes in self.TE_annot_attr_list:
            try:
                name = attributes[tag]
                if name not in name_list:
                    name_list.append(name)
            except KeyError:
                continue

        if len(name_list) == 0:
            name_list = ["undefined"]
        return name_list




    def calc_anchor_is_softclipped(self, read1, read2):
        #if the read is softclipped, return position at which it is clipped. If not, return None
        if self.read1_is_TE: #then read2 is anchor
            if read2.is_reverse and read2.cigar[0][0] == 4: #if the read is rev and softclipped to the left, ie into the interval
                self.anchor_softclipped_pos = read2.pos
            elif read2.cigar[-1][0] == 4: # if the read is fwd an softclipped at the right, ie into the interval
                self.anchor_softclipped_pos = read2.aend
        if self.read2_is_TE: # then read2 is anchor
            if read1.is_reverse and read1.cigar[0][0] == 4:
                self.anchor_softclipped_pos = read1.pos
            elif read1.cigar[-1][0] == 4:
                self.anchor_softclipped_pos = read1.aend
        return None



    def anchor_is_softclipped(self):
        return self.anchor_softclipped_pos

    def to_table(self, cluster_ID, library_name):
        if self.interval_direction == "fwd":
            mate_direction = "rev"
        else:
            mate_direction = "fwd"


        if self.read1_is_TE:
            read1_line = "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("R", cluster_ID, library_name, mate_direction, ".", ".", ".", "mate", self.read1_str)
            read2_line = "%s\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s" % ("R", cluster_ID, library_name, self.interval_direction, self.interval_start, self.interval_end, self.interval_chr, "anchor", self.read2_str)
        elif self.read2_is_TE:
            read2_line = "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("R", cluster_ID, library_name, mate_direction, ".", ".", ".", "mate", self.read2_str)
            read1_line = "%s\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s" % ("R", cluster_ID, library_name, self.interval_direction, self.interval_start, self.interval_end, self.interval_chr, "anchor", self.read1_str)
        return read1_line + "\n" + read2_line








