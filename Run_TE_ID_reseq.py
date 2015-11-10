
# elizabeth henaff
# Feb 2012
# elizabeth.henaff@cragenomica.es

import pysam
from BamReader import *
from ClusterList import *
from AlignedReadPair import *
import datetime
import subprocess
from pysam import csamtools
import os

import resource
import gc

from memory_profiler import profile


def reportResource(point=""):
    usage=resource.getrusage(resource.RUSAGE_SELF)
    return '''%s: usertime=%s systime=%s mem=%s mb
           '''%(point,usage[0],usage[1],
                (usage[2]*resource.getpagesize())/1000000.0 )


# @profile 
def run_jitterbug(psorted_bamfile_name, already_calc_discordant_reads, valid_discordant_reads_file_name, verbose, te_annot, \
    te_seqs, library_name, num_sdev, output_prefix, TE_name_tag, parallel, num_CPUs, bin_size, min_mapq, generate_test_bam, print_extra_output, conf_lib_stats, mem, min_cluster_size):


    mem_debug = False

    # print te_annot
    # print min_mapq

    

    #NOTE: comment this later !!!!!!!!!!!!!!!!
    #sorted_bam_reader = BamReader(output_prefix + ".proper_pair.sorted.bam", output_prefix)

    if mem_debug:
        reportResource("1")

    print("processing " + psorted_bamfile_name)
    if not output_prefix:
        output_prefix = psorted_bamfile_name



    # Make BamReader object with bam file of mapped and sorted reads
    # NOTE: uncomment this later !!!!!!!!!!!!!
    psorted_bam_reader = BamReader(psorted_bamfile_name, output_prefix)

    if generate_test_bam:
            print( "generating test bam")
            psorted_bam_reader.output_one_chr_reads()
            return None

    start_time = datetime.datetime.now()
    print( "starting at %s" % (str(start_time)))

    if conf_lib_stats:
        # Get the mean and sdev of insert size from supplied config file
        stats = {}
        for line in open(conf_lib_stats):
            line = line.strip()
            (tag, val) = line.split("\t")
            stats[tag] = (float(val))
        isize_mean = stats["fragment_length"]
        isize_sdev = stats["fragment_length_SD"]
        rlen_mean = stats["read_length"]
        rlen_sdev = stats["read_length_SD"]
        print( "mean fragment length taken from config file: %.2f" % (isize_mean))
        print( "standard deviation of fragment_length: %.2f" % (isize_sdev))
        print( "mean read length: %.2f" % (rlen_mean))
        print( "standard deviation of read length: %.2f" % (rlen_sdev))

    else:
        # Get the mean and sdev of insert size from the original bam file
        print("calculating mean insert size...")
        iterations = 1000000
        (isize_mean, isize_sdev, rlen_mean, rlen_sdev) = psorted_bam_reader.calculate_mean_sdev_isize(iterations)
        print( "mean fragment length over %d reads: %.2f" % (iterations, isize_mean))
        print( "standard deviation of fragment_length: %.2f" % (isize_sdev))
        print( "mean read length: %.2f" % (rlen_mean))
        print( "standard deviation of read length: %.2f" % (rlen_sdev))

        stats_file = open(output_prefix + ".read_stats.txt", "w")
        stats_file.write("fragment_length\t%.2f\n" % (isize_mean))
        stats_file.write("fragment_length_SD\t%.2f\n" % (isize_sdev))
        stats_file.write("read_length\t%.2f\n" % (rlen_mean))
        stats_file.write("read_length_SD\t%.2f" % (rlen_sdev))

        stats_file.close()

        #if fragment sdev is much larger than expected, there might be a problem with the reads of the mapping. Set as default 0.1*fragment length as a reasonable guess. 
        # This is necessary because aberrant values for sdev will mess up the interval overlap calculation and the filtering 
        if isize_sdev > 0.2*isize_mean:
            isize_sdev = 0.1*isize_mean
            print( "WARNING: fragment length standard deviation seems way too large to be realistic.\\n\
            There is maybe something weird with the flags in your bam mapping, or a very large number of large SV \\n\
            that are messing up the count.\\n\
            Setting the stdev to 0.1*fragment_length = %.2f for downstream calculations" % isize_sdev)



        time = datetime.datetime.now()
        print( "elapsed time: " + str(time - start_time))

    ################# Find valid discordant reads in given sorted bam file ################
    # This will print bam file(s) with the set(s) of valid discordant reads meaning
    # that the reads in  apair are mapped to distances greater than expected or to
    # two different chromosomes, and with at least one read that is mapped uniquely

    # if strict_repetitive is TRUE, will print out two bam files:
    # <bam_file_name>.valid_discordant_pairs_strict_rep.bam which is all valid discordant
    # read pairs with exactly one uniquely mapping and one repetitively mapping read pair
    #
    # <bam_file_name>.valid_discordant_pairs.bam which are all discordant valid discordant
    # read pairs with two uniquely mapping reads

    # if strict_repetitive is TRUE, will output both sets to a file named
    # <bam_file_name>.valid_discordant_pairs.bam

    #if have already run the prgm and calculated the discordant reads, dont do it again and look for a file called
    #<bam_file_name>.valid_discordant_pairs or <bam_file_name>.valid_discordant_pairs_strict_rep depending on the value of -s
    if not already_calc_discordant_reads:

        if mem_debug:
            reportResource("2")
        
        valid_discordant_reads_file_name = output_prefix + ".valid_discordant_pairs.bam"

        print( "selecting discordant reads...")
        #this writes the bam file of discordant reads to disc to be used later, and return the counts of different types of reads  
        (bam_stats, ref_lengths, ref_names) = psorted_bam_reader.select_discordant_reads_psorted( verbose, isize_mean, valid_discordant_reads_file_name)
        #print ref_names, ref_lengths
        coverage = (bam_stats["total_reads"] * rlen_mean )/ sum(ref_lengths)

        filter_conf_file = open(output_prefix + ".filter_config.txt", "w")
        filter_conf_file.write("cluster_size\t2\t%d\n" % (5*coverage))
        filter_conf_file.write("span\t2\t%d\n" % isize_mean)
        filter_conf_file.write("int_size\t%d\t%d\n" % (rlen_mean, 2*(isize_mean + 2*isize_sdev - (rlen_mean - rlen_sdev))) )
        filter_conf_file.write("softclipped\t2\t%d\n" % (5*coverage))
        filter_conf_file.write("pick_consistent\t0\t-1")

        filter_conf_file.close()

        bam_stats_file = open(output_prefix + ".bam_stats.txt", "w")
        for key, value in bam_stats.items():
            bam_stats_file.write("%s\t%d\n" % (key, value))

        if mem_debug:
            reportResource("3")
        bam_stats_file.close()

        time = datetime.datetime.now()
        if mem_debug:
            reportResource("4")
        print( "elapsed time: " + str(time - start_time))


#        print "sorting proper pair bam file in the background... "
#        args = ["samtools", "sort", output_prefix + ".proper_pair.bam", output_prefix + ".proper_pair.sorted"]
#        proper_pair_sort = subprocess.Popen(args)
    else:
        print( "using already selected discordant reads in %s" % (valid_discordant_reads_file_name))



    ################### Select valid discordant reads that match a TE #####################
    # Of the valid discordant read pairs, select those that have exactly one side that
    # overlaps a TE. This will
    # return a list of AlignedReadPair objects
    # the interval size in calculated as the INSIDE interval between the two reads ?? is this true? , plus numsdev * the sd of the insert size

    print( "selecting discordant read pairs where exactly one maps to a TE...")
    interval_size = isize_mean + num_sdev * isize_sdev

    if te_annot:
        discordant_bam_reader = BamReader(valid_discordant_reads_file_name, output_prefix)
        read_pair_one_overlap_TE_list = discordant_bam_reader.select_read_pair_one_overlap_TE_annot(te_annot, interval_size, min_mapq)

        if not (print_extra_output or already_calc_discordant_reads):
            os.remove(valid_discordant_reads_file_name)

    else:

        #here you would map mate reads to TE sequences and whatnot.
        pass

    if mem_debug:
            reportResource("5")
    time = datetime.datetime.now()
    print("elapsed time: " + str(time - start_time))

    ######################## wait till the proper pair bam file is sorted, and index it ###########################




    ##################### Cluster discordant read pairs that match TE #####################
    # Cluster the list of AlignedReadPair objects that have one read overlapping a TE
    # according to the uniquely mapping non-TE genomic location
    # then pair a fwd and rev cluster if they are overlapping
    # for the clusters that can be paired, calculate the softclipped support and the core reads, which will indicate heterozygosity of predicted TE insertion

    print( "generating clusters...")
    if len(read_pair_one_overlap_TE_list) == 0:
        print( "there might be an error, no discordant reads mapped to a TE location. please check the gff file... are the chromosome names the same as in the reference?")
        sys.exit(2)
    cluster_list = ClusterList(read_pair_one_overlap_TE_list)

    if not parallel:
        (cluster_pairs, unpaired_fwd_clusters, unpaired_rev_clusters, bed_string) = cluster_list.generate_clusters(verbose, psorted_bamfile_name, "", False, min_cluster_size)
        all_clusters = [(cluster_pairs, unpaired_fwd_clusters, unpaired_rev_clusters)]
    #parallel version:
    else:
        # last tow args are bed file handle and streaming, unnecessary and False, in this version 
        all_clusters = cluster_list.generate_clusters_parallel(verbose, num_CPUs, bin_size, psorted_bamfile_name, "", False, min_cluster_size)

    del read_pair_one_overlap_TE_list
    gc.collect()
    if mem_debug:
            reportResource("5")
    time = datetime.datetime.now()
    print( "elapsed time: " + str(time - start_time))


    ###################### print reads that were clustered to bam, output to gff and table ########################################


    print( "writing clustered reads to bam file, writing to gff and tables... ")


    pair_gff_output_file = open(output_prefix + ".TE_insertions_paired_clusters.gff3", "w")
    pair_table_output_file = open(output_prefix + ".TE_insertions_paired_clusters.supporting_clusters.table", "w")
    pair_table_output_file.write(table_header(library_name, library_name, te_annot))

    # if print_extra_output: 
    #     single_gff_output_file = open(output_prefix + ".TE_insertions_single_cluster.gff3", "w")
    #     single_table_output_file = open(output_prefix + ".TE_insertions_single_cluster.supporting_clusters.table", "w")
    #     single_table_output_file.write(table_header(library_name, library_name, te_annot))

    print(len(all_clusters))
    cluster_ID = 0
    for (cluster_pairs, unpaired_fwd_clusters, unpaired_rev_clusters) in all_clusters:


        # unpaired clusters are no longer reported
        # if print_extra_output:
        #     for fwd_cluster in unpaired_fwd_clusters:
        #         single_gff_output_file.write(fwd_cluster.to_gff(cluster_ID, library_name, TE_name_tag) + "\n")
        #         single_table_output_file.write(fwd_cluster.to_table(cluster_ID, library_name))

        #         cluster_ID += 1

        #     for rev_cluster in unpaired_rev_clusters:
        #         single_gff_output_file.write(rev_cluster.to_gff(cluster_ID, library_name, TE_name_tag) + "\n")
        #         single_table_output_file.write(rev_cluster.to_table(cluster_ID, library_name))
        #         cluster_ID += 1

        #print cluster_ID
        for cluster_pair in cluster_pairs:
            pair_gff_output_file.write(cluster_pair.to_gff(cluster_ID, library_name, TE_name_tag) + "\n")
            pair_table_output_file.write(cluster_pair.to_table(cluster_ID, library_name))
            cluster_ID += 1

    #clustered_reads_bam_file.close()
    time = datetime.datetime.now()
    print( "elapsed time: " + str(time - start_time))

    pair_gff_output_file.close()
    pair_table_output_file.close()

    # if print_extra_output:
    #     single_gff_output_file.close()
    #     single_table_output_file.close()

    
    end_time = str(datetime.datetime.now())
    print( "done! at " + end_time)

    run_stats = open(output_prefix + ".run_stats.txt", "w")
    run_stats.write("lib\t%s\n" % (library_name))
    if not already_calc_discordant_reads:
        run_stats.write("coverage\t%s\n" % (coverage))
    run_stats.write("runtime\t%s\n" % ( datetime.datetime.now() - start_time))
    run_stats.write("numCPUs\t%s\n" % (num_CPUs))
    run_stats.close()





def write_cluster_pairs_reads_to_bam(bam_file, cluster_pairs):

    for cluster_pair in cluster_pairs:
        read_pair_list = cluster_pair.fwd_cluster.readpair_list + cluster_pair.rev_cluster.readpair_list
        print(len(read_pair_list))
        for read_pair in read_pair_list:
            print( "read pair")
            print( str(read_pair.read1))
            print( str(read_pair.read2))
            bam_file.write(read_pair.read1)
            bam_file.write(read_pair.read2)



def write_single_clusters_to_bam(bam_file, fwd_clusters, rev_clusters):
    for cluster in fwd_clusters + rev_clusters:
        for read_pair in cluster.readpair_list:
            bam_file.write(read_pair.read1)
            bam_file.write(read_pair.read2)

