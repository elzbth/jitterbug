from Cluster import *
from AlignedReadPair import *

class ClusterPair:

    def __init__(self, fwd_c, rev_c):
        self.fwd_cluster = fwd_c
        self.rev_cluster = rev_c
        self.chr = self.fwd_cluster.chr

        (insertion_interval_start, insertion_interval_end) = self.minimal_interval()
        self.insertion_int_start = insertion_interval_start
        self.insertion_int_end = insertion_interval_end
        self.fwd_softclipped_pos = fwd_c.get_softclipped_pos()
        self.rev_softclipped_pos = rev_c.get_softclipped_pos()

        #list of positions where a read that overlaps the interval has a clipped tail. This list is separate from the list of softclipped reads that
        #are in the clusters, since the former are properly mapped and the latter are discordant
        self.softclipped_pos = []

        #list of tuples (start, end, numreads) where position is a softclipped position interval +/- 4bp and numreads is the number of reads that share that same softclipped position
        #sorted by decreasing support value, so the first one is the best one
        self.softclipped_support = None
        self.num_core_reads = -1

        #ratio of softclipped reads to softclipped reads + core reads. a ratio of 1 means fully homozygous for the TE insertion, a ratio near 0.5 means heterozygous.
        self.zygosity = -1.0

    def get_chr(self):
        return self.chr

    def get_insertion_int_start(self):
        return self.insertion_int_start

    def get_insertion_int_end(self):
        return self.insertion_int_end

    def set_zygosity(self, z):
        self.zygosity = z
        return 1

    def get_my_softclipped_pos(self):
        return self.softclipped_pos

    def get_fwd_cluster(self):
        return self.fwd_cluster

    def get_rev_cluster(self):
        return self.rev_cluster

    def set_my_softclipped_support(self, s):
        self.softclipped_support = s
        return 1

    def set_num_core_reads(self, n):
        self.num_core_reads = n
        return 1

    def append_softclipped_pos(self,pos):
        self.softclipped_pos.append(pos)
        return 1

    def get_clusters_softclipped_pos(self):
        return self.fwd_softclipped_pos + self.rev_softclipped_pos




    def to_table(self, cluster_ID, library_name):
        (insertion_interval_start, insertion_interval_end) = self.calc_insertion_interval()

        if self.softclipped_support:
#            pos_str_list = []
#            for (st, end, support) in self.softclipped_support:
#                pos_str_list.append( "(%d, %d, %d)" % (st, end, support))
#            pos_str = "[" + ", ".join(pos_str_list) + "]"
            (start, stop, support) = self.softclipped_support[0]
            pos_str = "%d\t%d\t%d" % (start, stop, support)

        else:
            pos_str = "-1\t-1\t0"
        insertion_line = "%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%.3f\n" % ("I", cluster_ID, library_name, self.chr, insertion_interval_start, insertion_interval_end, self.fwd_cluster.num_reads, self.rev_cluster.num_reads, self.fwd_cluster.span, self.rev_cluster.span, pos_str, self.num_core_reads, self.zygosity)
#        print "%d %d" % (self.fwd_cluster.span, self.rev_cluster.span)
        read_lines = ""
        for read_pair in self.fwd_cluster.readpair_list + self.rev_cluster.readpair_list:
            read_lines += read_pair.to_table(cluster_ID, library_name) + "\n"

        return insertion_line + "\n" + self.fwd_cluster.cluster_line(cluster_ID, library_name) + "\n" + self.rev_cluster.cluster_line(cluster_ID, library_name) + "\n" + read_lines


    def to_gff(self, cluster_pair_ID, library_name, TE_annot_tag):
        fwd_TE_tags = self.fwd_cluster.get_TE_tags(TE_annot_tag)
        rev_TE_tags = self.rev_cluster.get_TE_tags(TE_annot_tag)

        gff_tags = "supporting_fwd_reads=%d; supporting_rev_reads=%d; cluster_pair_ID=%d; lib=%s; Inserted_TE_tags_fwd=%s; Inserted_TE_tags_rev=%s; fwd_cluster_span=%d; rev_cluster_span=%d" \
                    % (self.fwd_cluster.num_reads, self.rev_cluster.num_reads, cluster_pair_ID, library_name, ", ".join(fwd_TE_tags), ", ".join(rev_TE_tags), self.fwd_cluster.span, self.rev_cluster.span)


        if self.softclipped_support:
            (st, end, num_support) = self.softclipped_support[0]
        else:
            (st, end, num_support) = (-1, -1, 0)
        gff_tags += "; softclipped_pos=(%d, %d); softclipped_support=%d; het_core_reads=%d; zygosity=%.3f" % (st, end, num_support, self.num_core_reads, self.zygosity)

        cluster_pair_chrom = self.fwd_cluster.chr

    #    #DEPRECATED: the coordinates of the region that the two clusters span, ie the UNION of the two clusters' prediction intervals
    #    #now use the insertion_interval coords calculated as the intersection
    #    cluster_pair_start = fwd_cluster[0].interval_start
    #    cluster_pair_end = rev_cluster[-1].interval_end


        (insertion_interval_start, insertion_interval_end) = self.calc_insertion_interval()

        coords = "%d\t%d" % (insertion_interval_start, insertion_interval_end)

        gff_line = "\t".join([cluster_pair_chrom, "jitterbug\tTE_insertion", coords, ".", ".", ".", gff_tags ])

        return gff_line

    def calc_insertion_interval(self):
        interval = self.minimal_interval()
        if interval == None:
            interval = self.intersection_interval_loose()
            print "using loose def of interval!"

        return interval

    def calc_softclipped_support(self, bam_file):

        #get softclipped reads that overlap my interval from the sofclipped bam file. these are concordant reads that have their tail softclipped
        reads = bam_file.fetch(self.chr, self.insertion_int_start, self.insertion_int_end)
        for read in reads:
            pos = softclipped_tail(read)
            if pos:
                self.softclipped_pos.append(pos)

        # to these, add the softclipped reads that were from discordant pairs actually forming the cluster
        softclipped_intervals = self.fwd_cluster.get_softclipped_pos() + self.rev_cluster.get_softclipped_pos() + self.softclipped_pos
        if len(softclipped_intervals) == 0:
            return None

        softclipped_intervals.sort()
        #for each softclipped position, count how many times it falls in the range within 4 bp more or less of another softclipped position
        pos_list = []
        #list of lists, each list containing tuples defining the st and end positions around the softclipped position.
        #each list is a set of mutually overlapping intervals
        pos_int_list = []

        current_int_list = []
        for position in softclipped_intervals:
            (st, end) = (position - 4, position + 5)

            if current_int_list == [] or st <= current_int_list[0][1]:
                current_int_list.append((st, end))
            else:
                pos_int_list.append(current_int_list)
                current_int_list = [(st, end)]
        pos_int_list.append(current_int_list)

        for int_list in pos_int_list:
            start = max(st for (st, end) in int_list)
            end = min(end for (st, end) in int_list)
            support = len(int_list)
            pos_list.append((start, end, support))
        #return list of (position, support) tuples, sorted by decreasing support
        pos_list.sort(key = lambda position : position[2], reverse=True)
        self.softclipped_support = pos_list


    def best_softclipped_support(self):
        return self.softclipped_support[0]
#        return (self.softclipped_support[0][0], self.softclipped_support[0][1])

    def calc_core_read_overlap(proper_pair_bam):
        if self.softclipped_support:
            (start, end, support) = self.best_softclipped_support()
            #get softclipped reads that overlap the softclipped interval from the proper pair bam file
            reads = bam_file.fetch(self.chr, start, end)
            for read in reads:
                if is_core(read, start, end):
                    self.num_core_reads += 1

            self.zygosity = float(support) / float(support + self.num_core_reads)

    def calc_zygosity(self,reads):



#        print "inside calc zygosity"

#        print len(list(reads))

        read_list = list(reads)
        #get softclipped reads that overlap my interval from the sofclipped bam file. these are concordant reads that have their tail softclipped
        #reads = bam_file.fetch(self.chr, self.insertion_int_start, self.insertion_int_end)

        for read in read_list:

            #don't want to consider reads that are not proper pairs
            if not read.is_proper_pair:
                continue
#            print "read"
            pos = softclipped_tail(read)
#            print "read2"
            if pos:
#                print "read3"
                self.append_softclipped_pos(pos)
#            print "read4"



        # to these, add the softclipped reads that were from discordant pairs actually forming the cluster
        softclipped_intervals =  self.get_clusters_softclipped_pos() + self.get_my_softclipped_pos()



        if len(softclipped_intervals) == 0:
            return 1

        softclipped_intervals.sort()
        #for each softclipped position, count how many times it falls in the range within 4 bp more or less of another softclipped position
        pos_list = []
        #list of lists, each list containing tuples defining the st and end positions around the softclipped position.
        #each list is a set of mutually overlapping intervals
        pos_int_list = []

        current_int_list = []
        for position in softclipped_intervals:
            (st, end) = (position - 4, position + 5)

            if current_int_list == [] or st <= current_int_list[0][1]:
                current_int_list.append((st, end))
            else:
                pos_int_list.append(current_int_list)
                current_int_list = [(st, end)]
        pos_int_list.append(current_int_list)

        for int_list in pos_int_list:
            start = max(st for (st, end) in int_list)
            end = min(end for (st, end) in int_list)
            support = len(int_list)
            pos_list.append((start, end, support))
        #return list of (start, end, support) tuples, sorted by decreasing support
        pos_list.sort(key = lambda position : position[2], reverse=True)
        self.set_my_softclipped_support(pos_list)

        (start, end, support) = pos_list[0]

#        print "calc core reads..."
        num_core_reads = 0
        for read in read_list:
            #don't want to consider reads that are not proper pairs
            if not read.is_proper_pair:
                continue

            if is_core(read, start, end):
                num_core_reads += 1
        self.set_num_core_reads(num_core_reads)

        self.set_zygosity( float(support) / float(support + num_core_reads))

#        print "done calc zygosity"

        return 1

    def intersection_interval_loose(self):
        #the coordinates of the INTERSECTION of the two cluster's total prediction intervals

        insertion_interval_start = self.rev_cluster.cluster_start
        insertion_interval_end = self.fwd_cluster.cluster_end
        return (insertion_interval_start, insertion_interval_end)

    def union_interval_tight(self):
        # the intersection of the two minimal intervals, of fwd cluster and rev cluster
        fwd_intersection_range = list(range(self.fwd_cluster.intersection_start, self.fwd_cluster.intersection_end + 1))
        rev_intersection_range = list(range(self.rev_cluster.intersection_start, self.rev_cluster.intersection_end + 1))

        union = fwd_intersection_range + rev_intersection_range

        if union == None:
            return None



        union.sort()
        print union


        if len(union) == 1:
            return [union[0], union[0] + 1]
        else:
            return [union[0], union[-1]]

    def minimal_interval(self):
        start = int(max(readpair.interval_start for readpair in self.fwd_cluster.readpair_list))
        end = int(min(readpair.interval_end for readpair in self.rev_cluster.readpair_list))

        if start > end:
            return None

        else:
            return [start, end]

def softclipped_tail(read):
    """returns the psition at which this read is softclipped at its tail, or None otherwise"""
    #if read is reverse and the last cigar character is S
    if read.is_reverse and read.cigar[-1][0] == 4:
        #return the mapped end position
        return read.aend
    #else if read is fwd and the first cigar character is S
    elif read.cigar[0][0] == 4:
        return read.pos
    else:
        return None

def is_core(read, start, end):
    """returns true if this read overlaps the interval defined by start and end in its "core", ie not in the first or last 7bp, which can be mismapped with a bunch of SNPs"""
    if read.pos < (start - 15) and read.aend > (end + 15):
        return True
    else:
        return False
