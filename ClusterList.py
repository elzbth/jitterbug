from collections import deque
from Cluster import *
from ClusterPair import *
#import pp
# from pp import *
import pysam
#from pysam import csamtools

from multiprocessing import Pool
import multiprocessing as mp 


class ClusterList:

    def __init__(self, read_pair_list):
        #list of AlignedReadPair objects
        self.read_pair_list = read_pair_list

        #list of Cluster Objects
        self.cluster_list = []

        #list of ClusterPair objects



    #cluster the read pairs according to the interval defined by the non-TE mapped read
    def generate_clusters_parallel(self, verbose, num_CPUs, bin_size, psorted_bamfile_name, bed_file_handle, streaming, min_cluster_size):




        ################ CLUSTER BY BIN #######################
        #cluster fwd intervals
        fwd_read_pairs = [read_pair for read_pair in self.read_pair_list if read_pair.interval_direction == "fwd"]
        fwd_clusters_by_bin = cluster_read_pairs_by_chr_and_bin(fwd_read_pairs, bin_size)
        print("********************* total fwd non-overlapping clusters found by bin: %d" % sum([len(chr_list) for chr_list in fwd_clusters_by_bin.values()]))



        #cluster rev intervals
        rev_read_pairs = [read_pair for read_pair in self.read_pair_list if read_pair.interval_direction == "rev"]
        rev_clusters_by_bin = cluster_read_pairs_by_chr_and_bin(rev_read_pairs, bin_size)
        print("********************* total rev non-overlapping clusters found by bin: %d" % sum([len(chr_list) for chr_list in rev_clusters_by_bin.values()]))


        ############################ END CLUSTER BY BIN ###########################






        ##########TODO THIS IS DEBUG

#        #cluster fwd intervals
#        fwd_read_pairs = [read_pair for read_pair in self.read_pair_list if read_pair.interval_direction == "fwd"]
#        fwd_clusters = cluster_read_pairs_all(fwd_read_pairs)
#
#        print "******************total fwd clusters found: %d" %  len(fwd_clusters)
#        non_overlapping_fwd_clusters = remove_overlapping_clusters(fwd_clusters)
#        print "******************total fwd non-overlapping clusters found: %d" %  len(non_overlapping_fwd_clusters)
#
#
#        #cluster rev intervals
#        rev_read_pairs = [read_pair for read_pair in self.read_pair_list if read_pair.interval_direction == "rev"]
#        rev_clusters = cluster_read_pairs_all(rev_read_pairs)
#
#        print "******************total rev clusters found: %d" % len(rev_clusters)
#        non_overlapping_rev_clusters = remove_overlapping_clusters(rev_clusters)
#        print "******************total rev non-overlapping clusters found: %d" %  len(non_overlapping_rev_clusters)
#        print fwd_clusters_by_chr.keys()
#        print rev_clusters_by_chr.keys()

        ############################## END THIS IS DEBUG




        ############## START NEW PARALLEL VERSION #################################

        input_arg_list = []

        #print fwd_clusters_by_bin.keys()
        #print rev_clusters_by_bin.keys()
        for key in fwd_clusters_by_bin.keys():
            if key in rev_clusters_by_bin.keys():
                input_arg_list.append((key, fwd_clusters_by_bin[key], rev_clusters_by_bin[key], psorted_bamfile_name, verbose, streaming, min_cluster_size))

        #print input_arg_list[0][0][0].readpair_list[0].read1

        print("sending %d jobs to %d processes" % (len(input_arg_list), num_CPUs))

        mp.set_start_method('spawn')

        pool = Pool(num_CPUs)

        # dummy_arg_list = [("string1", "string2")] * len(input_arg_list)

        # all_clusters_by_bin = pool.imap(dummy_func, input_arg_list)

        all_clusters_by_bin = pool.map(pair_clusters_by_bin, input_arg_list)

        pool.close()
        pool.join()


        #for mem debug
        #return list(all_clusters_by_bin)


        ################ END NEW PARALLEL VERSION #################################
        # (paired, fwd, rev, bed_strings) = all_clusters_by_bin[0:3]

        # print all_clusters_by_bin

        if streaming:

            cluster_counts = [(len(p), len(f), len(r)) for (p,f,r,s) in all_clusters_by_bin]
            print("******************total fwd single clusters found: %d" %  sum([f for (p,f,r) in cluster_counts]))
            print("******************total rev single clusters found: %d" %  sum([r for (p,f,r) in cluster_counts]))
            print("******************total cluster pairs found: %d" %  sum([p for (p,f,r) in cluster_counts]))

            bed_string = "\n".join([s for (p,f,r,s) in all_clusters_by_bin if s != ""])


            # print bed_string
            bed_file_handle.write(bed_string)
            bed_file_handle.close()

            

        else:
            cluster_counts = [(len(p), len(f), len(r)) for (p,f,r) in all_clusters_by_bin]
            print("******************total fwd single clusters found: %d" %  sum([f for (p,f,r) in cluster_counts]))
            print("******************total rev single clusters found: %d" %  sum([r for (p,f,r) in cluster_counts]))
            print("******************total cluster pairs found: %d" %  sum([p for (p,f,r) in cluster_counts]))

            

        return all_clusters_by_bin


##################### END PARALLEL VERSION #############################################


    def generate_clusters(self, verbose, psorted_bamfile_name, bed_file_handle, streaming, min_cluster_size):
##################### BEGIN NON PARALLEL VERSION ######################################
        #cluster fwd intervals
        fwd_read_pairs = [read_pair for read_pair in self.read_pair_list if read_pair.interval_direction == "fwd"]
        fwd_clusters = cluster_read_pairs_all(fwd_read_pairs)

        print("******************total fwd clusters found: %d" %  len(fwd_clusters))
        non_overlapping_fwd_clusters = remove_overlapping_clusters(fwd_clusters)
        print("******************total fwd non-overlapping clusters found: %d" %  len(non_overlapping_fwd_clusters))


        #cluster rev intervals
        rev_read_pairs = [read_pair for read_pair in self.read_pair_list if read_pair.interval_direction == "rev"]
        rev_clusters = cluster_read_pairs_all(rev_read_pairs)

        print("******************total rev clusters found: %d" % len(rev_clusters))
        non_overlapping_rev_clusters = remove_overlapping_clusters(rev_clusters)
        print("******************total rev non-overlapping clusters found: %d" %  len(non_overlapping_rev_clusters))

        #bam_file_name = output_prefix + ".proper_pair.sorted.bam"
        psorted_bamfile = pysam.Samfile(psorted_bamfile_name, "rb")


        #pair clusters by genomic location, keeping track of which indices in the array have been paired, so that you can pick out the unpaired ones after
        cluster_pairs = []
        paired_fwd_clusters_indices = []
        paired_rev_clusters_indices = []
        bed_string = ""

        # iterate over combinations of fwd and rev clusters, skipping if clusters dont meet min size requirements
        for fwd_index, fwd_cluster in enumerate(non_overlapping_fwd_clusters):
            if fwd_cluster.num_reads < min_cluster_size:
                continue

            for rev_index, rev_cluster in enumerate(non_overlapping_rev_clusters):
                if rev_cluster.num_reads < min_cluster_size:
                    continue

                if fwd_cluster.is_overlapping_strict(rev_cluster):
                    new_cluster_pair = ClusterPair(fwd_cluster, rev_cluster)
                    if not streaming:
                        reads = psorted_bamfile.fetch(new_cluster_pair.get_chr(), new_cluster_pair.get_insertion_int_start(), new_cluster_pair.get_insertion_int_end())
                        new_cluster_pair.calc_zygosity(reads)
                    else:
                        bed_line = new_cluster_pair.to_bed()
                        if bed_string == "":
                            bed_string = bed_line
                        else:
                            bed_string = bed_string + "\n" + bed_line
                    if new_cluster_pair.insertion_int_end < new_cluster_pair.insertion_int_start:
                        if True:
                            print("cluster pair not paired!")
                    else:
                        cluster_pairs.append(new_cluster_pair)
                        paired_fwd_clusters_indices.append(fwd_index)
                        paired_rev_clusters_indices.append(rev_index)

        #make lists of unpaired clusters
        # unpaired_fwd_clusters = []
        # unpaired_rev_clusters = []
        # for fwd_index in range(len(non_overlapping_fwd_clusters)):
        #     if fwd_index not in paired_fwd_clusters_indices:
        #         unpaired_fwd_clusters.append(non_overlapping_fwd_clusters[fwd_index])

        # for rev_index in range(len(non_overlapping_rev_clusters)):
        #     if rev_index not in paired_rev_clusters_indices:
        #         unpaired_rev_clusters.append(non_overlapping_rev_clusters[rev_index])


        # print bed_string
        if streaming: 
            bed_file_handle.write(bed_string)
            bed_file_handle.close()


        print("******************total cluster pairs found: %d" %  len(cluster_pairs))
        if verbose:
            for (fwd_cluster, rev_cluster) in cluster_pairs:
                print("*************************cluster_pair:**************************************")
                print("fwd cluster:")
                print("cluster coordinates: %s %d %d" % (fwd_cluster[0].interval_chr, fwd_cluster[0].interval_start, fwd_cluster[-1].interval_end ))
                print(" ".join(read.str_int() for read in fwd_cluster))
                print(" ".join(read.str_TE_annot_list() for read in fwd_cluster))
                print("rev cluster:")
                print("cluster coordinates: %s %d %d" % (rev_cluster[0].interval_chr, rev_cluster[0].interval_start, rev_cluster[-1].interval_end ))
                print(" ".join(read.str_int() for read in rev_cluster))
                print(" ".join(read.str_TE_annot_list() for read in rev_cluster))

        return (cluster_pairs, [], [], bed_string)


############################### END NON PARALLEL VERSION ########################################################





def pair_clusters_by_bin(param_list):

    (key, fwd_clusters, rev_clusters, bam_file_name, verbose, streaming, min_cluster_size) = param_list


    print("processing cluster pairs on %s" % (key))
    #print "pairing clusters in parallel for chr %s" % fwd_clusters[0].chr
    non_overlapping_fwd_clusters = remove_overlapping_clusters(fwd_clusters)
    if verbose:
        print("non overlapping fwd clusters\t%d" % (len(non_overlapping_fwd_clusters)))
    non_overlapping_rev_clusters = remove_overlapping_clusters(rev_clusters)
    if verbose:
        print("non overlapping rev clusters\t%d" % (len(non_overlapping_rev_clusters)))
    if not streaming:
        proper_pair_bam = pysam.AlignmentFile(bam_file_name, "rb")
    #print "haha"


    #print "ok1"
    #pair clusters by genomic location, keeping track of which indices in the array have been paired, so that you can pick out the unpaired ones after
    cluster_pairs = []
    paired_fwd_clusters_indices = []
    paired_rev_clusters_indices = []
    bed_string = ""
    for fwd_index, fwd_cluster in enumerate(non_overlapping_fwd_clusters):
        if fwd_cluster.num_reads < min_cluster_size:
                continue

        for rev_index, rev_cluster in enumerate(non_overlapping_rev_clusters):
            if rev_cluster.num_reads < min_cluster_size:
                    continue

            if fwd_cluster.is_overlapping_strict(rev_cluster):
                new_cluster_pair = ClusterPair(fwd_cluster, rev_cluster)
                #print new_cluster_pair.get_chr()
                if not streaming:
                    reads = proper_pair_bam.fetch(new_cluster_pair.get_chr(), new_cluster_pair.get_insertion_int_start(), new_cluster_pair.get_insertion_int_end())
                    new_cluster_pair.calc_zygosity(reads)
                else:
                    bed_line = new_cluster_pair.to_bed()
                    bed_string = bed_string + "\n" + bed_line
                #print "poop"
                if new_cluster_pair.get_insertion_int_end() < new_cluster_pair.get_insertion_int_start():
                    if True:
                        print("cluster pair not paired!")
                else:
                    cluster_pairs.append(new_cluster_pair)
                    paired_fwd_clusters_indices.append(fwd_index)
                    paired_rev_clusters_indices.append(rev_index)

    #make lists of unpaired clusters
    # unpaired_fwd_clusters = []
    # unpaired_rev_clusters = []
    # for fwd_index in range(len(non_overlapping_fwd_clusters)):
    #     if fwd_index not in paired_fwd_clusters_indices:
    #         unpaired_fwd_clusters.append(non_overlapping_fwd_clusters[fwd_index])

    # for rev_index in range(len(non_overlapping_rev_clusters)):
    #     if rev_index not in paired_rev_clusters_indices:
    #         unpaired_rev_clusters.append(non_overlapping_rev_clusters[rev_index])


    # if streaming:
    #     return (cluster_pairs, unpaired_fwd_clusters, unpaired_rev_clusters, bed_string)
    # else:
    #     return (cluster_pairs, unpaired_fwd_clusters, unpaired_rev_clusters)

    if streaming:
        return (cluster_pairs, [], [], bed_string)
    else:
        return (cluster_pairs, [], [])





#helper functions
def cluster_read_pairs_by_chr(read_pair_list):
    """this generates a list of  maximal  clusters, ie sets of overlapping read pairs. note: these clusters can be themselves overlapping.
    returns a disctionary of lists of Cluster objects, one entry per chromosome"""

    #sort according to end position then chromosome. sort is stable so the second sort will not unsort the positions
    read_pair_list.sort(key=lambda read_pair: read_pair.interval_end)
    read_pair_list.sort(key=lambda read_pair: read_pair.interval_chr)



    #store each a list of current Cluster objects, which contains a list of AlignedReadPair objects
    cluster_list = []

    #store a seperate cluster_list for each chromosome
    chr_cluster_lists = {}

    #read_pair_Q stores a list of currently overlapping read pair intervals
    read_pair_Q = deque([read_pair_list[0]])

    for read_pair in read_pair_list:
        #print read_pair_Q

        #if you can add the next interval to the current list of overlapping intervals, do so
        if read_pair.interval_chr == read_pair_Q[0].interval_chr and read_pair.interval_start <= read_pair_Q[0].interval_end:
            read_pair_Q.append(read_pair)

        #if the current read is from another chromosome, save the list of currently overlapping intervals as a cluster and empty it
        elif read_pair.interval_chr != read_pair_Q[0].interval_chr:
            new_cluster = Cluster(list(read_pair_Q))
            cluster_list.append(new_cluster)
            chr_cluster_lists[read_pair_Q[0].interval_chr] = cluster_list

            #empty queue and current cluster list since we are starting with a new chromosome
            cluster_list = []
            read_pair_Q.clear()
            read_pair_Q.append(read_pair)

        #otherwise, save the list of currently overlapping intervals as a cluster
        else:
            new_cluster = Cluster(list(read_pair_Q))
            cluster_list.append(new_cluster)
            # and pop off intervals in the Q as long as they do not overlap with your current interval -- these cannot constitute another maximal cluster
            while len(read_pair_Q) != 0 and read_pair.interval_start > read_pair_Q[0].interval_end:
                read_pair_Q.popleft()
            #then add your current read to the Q
            read_pair_Q.append(read_pair)
    #for cluster in cluster_list:
    #    print " ".join(read.str_int() for read in cluster)
        #print " ".join(read.str_TE_annot_list() for read in cluster)
    last_cluster = Cluster(list(read_pair_Q))
    cluster_list.append(last_cluster)
    chr_cluster_lists[read_pair_Q[0].interval_chr] = cluster_list

    return chr_cluster_lists



def cluster_read_pairs_by_chr_and_bin(read_pair_list, bin_size):
    """this generates a list of  maximal  clusters, ie sets of overlapping read pairs. note: these clusters can be themselves overlapping.
    returns a disctionary of lists of Cluster objects, one entry per bin"""

    # if no bin is set, cluster by chromosome
    if bin_size == None:
        return cluster_read_pairs_by_chr(read_pair_list)

    #sort according to end position then chromosome. sort is stable so the second sort will not unsort the positions
    read_pair_list.sort(key=lambda read_pair: read_pair.interval_end)
    read_pair_list.sort(key=lambda read_pair: read_pair.interval_chr)



    current_bin_start = 0
    current_bin_end = bin_size - 1

    current_chr = read_pair_list[0].interval_chr

    current_bin_key = "%s %d-%d" % (current_chr, current_bin_start, current_bin_end)




    #store each a list of current Cluster objects, which contains a list of AlignedReadPair objects
    current_bin_cluster_list = []

    #store a seperate cluster_list for each bin
    bin_cluster_lists = {}

    #read_pair_Q stores a list of currently overlapping read pair intervals
    read_pair_Q = deque([read_pair_list[0]])

    for read_pair in read_pair_list:
        #print read_pair_Q

        #if you can add the next interval to the current list of overlapping intervals, do so
        if read_pair.interval_chr == read_pair_Q[0].interval_chr and read_pair.interval_start <= read_pair_Q[0].interval_end:
            read_pair_Q.append(read_pair)

        #if the current read is from another chromosome, or from another bin, save the list of currently overlapping intervals as a cluster and empty it
        elif read_pair.interval_chr != current_chr or read_pair.interval_start > current_bin_end:
            new_cluster = Cluster(list(read_pair_Q))
            current_bin_cluster_list.append(new_cluster)
            bin_cluster_lists[current_bin_key] = current_bin_cluster_list

            #empty queue and current cluster list since we are starting with a new bin
            current_bin_cluster_list = []
            read_pair_Q.clear()
            read_pair_Q.append(read_pair)

            #update the current chromosome and bins

            if read_pair.interval_chr != current_chr:
                current_bin_start = 0
                current_bin_end = bin_size - 1

                current_chr = read_pair.interval_chr
            else:
                current_bin_start = current_bin_end + 1
                current_bin_end = current_bin_start + bin_size - 1

            current_bin_key = "%s %d-%d" % (current_chr, current_bin_start, current_bin_end)

        #otherwise, save the list of currently overlapping intervals as a cluster
        else:
            new_cluster = Cluster(list(read_pair_Q))
            current_bin_cluster_list.append(new_cluster)
            # and pop off intervals in the Q as long as they do not overlap with your current interval -- these cannot constitute another maximal cluster
            while len(read_pair_Q) != 0 and read_pair.interval_start > read_pair_Q[0].interval_end:
                read_pair_Q.popleft()
            #then add your current read to the Q
            read_pair_Q.append(read_pair)
    #for cluster in cluster_list:
    #    print " ".join(read.str_int() for read in cluster)
        #print " ".join(read.str_TE_annot_list() for read in cluster)
    last_cluster = Cluster(list(read_pair_Q))
    current_bin_cluster_list.append(last_cluster)
    bin_cluster_lists[current_bin_key] = current_bin_cluster_list

    return bin_cluster_lists

def cluster_read_pairs_all(read_pair_list):
    """this generates a list of  maximal  clusters, ie sets of overlapping read pairs. note: these clusters can be themselves overlapping.
    returns a disctionary of lists of Cluster objects, one entry per chromosome"""

    #sort according to end position then chromosome. sort is stable so the second sort will not unsort the positions
    read_pair_list.sort(key=lambda read_pair: read_pair.interval_end)
    read_pair_list.sort(key=lambda read_pair: read_pair.interval_chr)



    #store each cluster as a list of AlignedReadPair objects
    cluster_list = []



    #read_pair_Q stores a list of currently overlapping read pair intervals
    read_pair_Q = deque([read_pair_list[0]])

    for read_pair in read_pair_list:
        #print read_pair_Q

        #if you can add the next interval to the current list of overlapping intervals, do so
        if read_pair.interval_chr == read_pair_Q[0].interval_chr and read_pair.interval_start <= read_pair_Q[0].interval_end:
            read_pair_Q.append(read_pair)

        #if the current read is from another chromosome, save the list of currently overlapping intervals as a cluster and empty it
        elif read_pair.interval_chr != read_pair_Q[0].interval_chr:
            new_cluster = Cluster(list(read_pair_Q))
            cluster_list.append(new_cluster)
            read_pair_Q.clear()
            read_pair_Q.append(read_pair)

        #otherwise, save the list of currently overlapping intervals as a cluster
        else:
            new_cluster = Cluster(list(read_pair_Q))
            cluster_list.append(new_cluster)
            # and pop off intervals in the Q as long as they do not overlap with your current interval -- these cannot constitute another maximal cluster
            while len(read_pair_Q) != 0 and read_pair.interval_start > read_pair_Q[0].interval_end:
                read_pair_Q.popleft()
            #then add your current read to the Q
            read_pair_Q.append(read_pair)

    last_cluster = Cluster(list(read_pair_Q))
    cluster_list.append(last_cluster)
    #for cluster in cluster_list:
    #    print " ".join(read.str_int() for read in cluster)
        #print " ".join(read.str_TE_annot_list() for read in cluster)

    return cluster_list



def remove_overlapping_clusters(cluster_list):
    """returns a list of clusters that do not overlap with any other. input is a list of lists of AlignedReadPair objects, sorted by end position.
    thus the start coordinate of the cluster will be the start coordinate of its first element: cluster[0]
    and the end coordinate of the cluster will be the end coordinate of its last element: cluster[-1]"""

    non_overlapping_clusters = []

    current_cluster = cluster_list[0]
    current_cluster_is_overlapped = False
    next_cluster_is_overlapped = False

    for next_cluster in cluster_list[1:]:
        #if the current cluster does not overlap teh next one,
        if current_cluster.cluster_end < next_cluster.cluster_start:
            next_cluster_is_overlapped = False
            #and is not overlapped itself
            if not current_cluster_is_overlapped:
                #add it to the list
                non_overlapping_clusters.append(current_cluster)
        #otherwise, flag the next cluster as overlapped
        else:
            next_cluster_is_overlapped = True
        #update current to next
        current_cluster = next_cluster
        current_cluster_is_overlapped = next_cluster_is_overlapped

    # add last cluster if it is not overlapped
    if not current_cluster_is_overlapped:
        #add it to the list
        non_overlapping_clusters.append(current_cluster)

    #for cluster in non_overlapping_clusters:
    #    print " ".join(read.str_int() for read in cluster)

    return non_overlapping_clusters










def table_header(library_name,  bam_file_name, te_annot):
    param_string = "#this table describes the read clusters identified in the bam file %s and corresponding to the transposon annotations in %s\n" % (bam_file_name, te_annot)
    title_string = "#this table contains three types of lines:\
#** insertion lines: one per predicted insertion site, corresponding to a pair of overlapping clusters, one fwd, one rev\n\
I\tcluster_pair_ID\tlib\tchrom\tstart\tend\tnum_fwd_reads\tnum_rev_reads\tfwd_span\trev_span\tbest_sc_pos_st\tbest_sc_pos_end\tsc_pos_support\n\
#here the start and end are defined as the intersection of the intervals predicted by the leftmost forward read and the rightmost reverse read.\n\n\
#** cluster lines (two per insertion, one fwd and one rev):\n\
C\tcluster_pair_ID\tlib\tdirection\tstart\tend\tchrom\tnum_reads\tspan\n\
#span is defined as the range of start positions in the cluster. A span of 0 means that all the reads originate at the same start site, and are probably an artifact. \
#a span the size of the fragment length indicates good coverage. \n\n\
#**read lines (fwd reads consitute the fwd clusters, rev reads the rev clusters)\n\
#the reads that are \"anchor\" are those that consitute the cluster, the reads that are \"mate\" are the anchors' mates, which map to a TE\n\
R\tcluster_pair_ID\tlib\tdirection\tinterval_start\tinterval_end\tchrom\tstatus\tbam_line\n\n\
#this file is meant to be easily manipulated with tools like grep and sed, for example\n\
#grep ^C table_file > cluster_pairs.out\n\
#will give you a list of all the clusters pairs\n\
#the R lines sharing the same ID all come from the same cluster pair, with itself the same ID, corresponding to the the insertion of that same ID, thus \n\
#grep -w cluster_pair_ID_X table_file > predicted_insertion_X.table.out\n\
#will give you the definition line of the the predicted insertion site, the fwd and rev clusters comprising insertion X and the description of the reads that constitute them. \n"

    return param_string + title_string
















