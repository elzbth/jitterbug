class Cluster:


    def __init__(self, readpair_list):
        self.readpair_list = readpair_list
        self.direction = readpair_list[0].interval_direction
        self.cluster_start = int(min(read.interval_start for read in readpair_list))
        self.cluster_end = int(max(read.interval_end for read in readpair_list))
        self.chr = readpair_list[0].interval_chr
        self.num_reads = len(readpair_list)
        self.span = self.readpair_list[-1].interval_start - self.readpair_list[0].interval_start
        self.intersection_start = int(max(readpair.interval_start for readpair in readpair_list))
        self.intersection_end = int(min(readpair.interval_end for readpair in readpair_list))
        self.softclipped_pos = []
        for read_pair in self.readpair_list:
            pos = read_pair.anchor_is_softclipped()
            if pos:
                self.softclipped_pos.append(pos)

    def get_softclipped_pos(self):

        return self.softclipped_pos

    def get_TE_tags(self, TE_annot_tag):
        tags = []
        for read_pair in self.readpair_list:
            for value in read_pair.TE_annot_tag_list(TE_annot_tag):
                if value not in tags:
                    tags.append(value)
        return tags

    def to_gff(self, cluster_ID, library_name, TE_annot_tag):

        gff_tags = "supporting_reads=%d; cluster_ID=%d; lib=%s; Inserted_TE_tags=%s; direction=%s; span=%d" % (self.num_reads, cluster_ID, library_name, ", ".join(self.get_TE_tags(TE_annot_tag)), self.direction, self.span)

        coords = "%d\t%d" % (self.cluster_start, self.cluster_end)

        gff_line = "\t".join([self.chr, "TE_insertion_single_cluster\tjitterbug", coords, ".", ".", ".", gff_tags ])

        return gff_line

    def to_table(self, cluster_ID, library_name):
        read_lines = ""
        for read_pair in self.readpair_list:
            read_lines += read_pair.to_table(cluster_ID, library_name) + "\n"
        return self.cluster_line(cluster_ID, library_name) + "\n" + read_lines

    def cluster_line(self, cluster_ID, library_name):
        return "%s\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\n" % ("C", cluster_ID, library_name, self.direction, self.cluster_start, self.cluster_end, self.chr, self.num_reads, self.span)

    def is_overlapping(self, other):
        if self.chr != other.chr:
            return False

        if self.cluster_start in range(other.cluster_start, other.cluster_end + 1):
            return True
        if self.cluster_end in range(other.cluster_start, other.cluster_end + 1):
            return True
        if other.cluster_start in range(self.cluster_start, self.cluster_end + 1):
            return True
        if other.cluster_end in range(self.cluster_start, self.cluster_end + 1):
            return True
        return False

    def is_overlapping_strict(self, other):
        if self.chr != other.chr:
            return False

        ## new version ###
        if other.intersection_start <= self.intersection_start and self.intersection_start <= other.intersection_end:
            return True
        if other.intersection_start <= self.intersection_end and self.intersection_end <= other.intersection_end:
            return True
        if self.intersection_start <= other.intersection_start and other.intersection_start <= self.intersection_end:
            return True
        if self.intersection_start <= other.intersection_end and other.intersection_end <= self.intersection_end:
            return True
        ### end new version ###

        #### old version #####
        # if self.intersection_start in range(other.intersection_start, other.intersection_end + 1):
        #     return True
        # if self.intersection_end in range(other.intersection_start, other.intersection_end + 1):
        #     return True
        # if other.intersection_start in range(self.intersection_start, self.intersection_end + 1):
        #     return True
        # if other.intersection_end in range(self.intersection_start, self.intersection_end + 1):
        #     return True
        #### end old version ###
        return False

