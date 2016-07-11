#!/usr/bin/env python

# for CRG cluster users:
##!/software/so/el6.3/PythonPackages-2.7.6/bin/python

# old CRG cluster paths:
######!/software/so/el6.3/PythonPackages-2.7.3-VirtualEnv/bin/python
#####!/software/so/el6.3/PythonPackagesVirtualEnv/bin/python

# path for loki: 
##!/soft/bin/python

# elizabeth henaff
# Feb 2012
# elizabeth.m.henaff@gmail.com



import getopt,sys
import argparse
from Run_TE_ID_reseq import *
from Run_TE_ID_reseq_streaming import *



def main(argv):



    parser = argparse.ArgumentParser()

    #required arguments
    parser.add_argument("mapped_reads", help="reads mapped to reference genome in bam format, sorted by position")
    parser.add_argument("TE_annot", help="annotation of transposable elements in reference genome, in gff3 format")

    #optional arguments - flags
    parser.add_argument("-v", "--verbose", help="print more output to the terminal", action="store_true")
    #flags for debug
    # parser.add_argument("--keep_temp", help="do not delete temp files: bam of discordant reads", action="store_true")
    parser.add_argument("--pre_filter", help="pre-filter reads with samtools, and save intermediate filtered read subset", action="store_true")
    # parser.add_argument("--mem", help="bam file has been mapped with bwa-mem and you want to use split (=chimeric=supplementary) reads. If unset, bam files \
    #    mapped with either bwa or bwa-mem can be used, but split reads will be ignored", action="store_true")

    #optional arguments
    parser.add_argument("-l", "--lib_name", help="sample or library name, to be included in final gff output")
    parser.add_argument("-d","--sdev_mult", type=int, help="use SDEV_MULT*fragment_sdev + fragment_length when calculating insertion intervals. Best you don't touch this." , default=2)
    parser.add_argument("-o", "--output_prefix", help="prefix of output files. Can be path/to/directory/file_prefix", default="jitterbug")
    parser.add_argument("-n", "--numCPUs", type=int, help="number of CPUs to use", default=1)
    parser.add_argument("-b", "--bin_size", type=int, help="If parallelized, size of bins to use, in bp. If numCPUs > 1 and bin_size == 0, will parallelize by entire chromosomes", default=50000000)
    parser.add_argument("-q", "--minMAPQ", type=int, help="minimum read mapping quality to be considered", default=15)
    parser.add_argument("-t", "--TE_name_tag", help="name of tag in TE annotation gff file to use to record inserted TEs", default="Name")
    parser.add_argument("-s", "--conf_lib_stats", help="tabulated config file that sets the values to use for fragment length and sdev, read length and sdev. \
        4 tab-deliminated lines: key\tvalue, keys are:fragment_length, fragment_length_SD, read_length,  read_length_SD")
    parser.add_argument("-c", "--min_cluster_size", type=int, help="min number of both fwd and rev reads to predict an insertion", default=2)
    parser.add_argument("--step_one_only", dest='step_one_only', action='store_true', help="if specified it stops execution after the first BAM filtering", default=False)
    parser.add_argument("--step_two_only", dest='step_two_only', action='store_true', help="if specified it starts exection after the first BAM filtering", default=False)


    #optional arguments for debug
    parser.add_argument("--disc_reads_bam", help="for debug. Use as input bam file of discordant reads only (generated by running with --pre_filter), and skip the step of perusing the input bam for discordant reads")

    args = parser.parse_args()

    fail_string = "use jitterbug.py -h for complete option list"

    #check args

    if not os.path.exists(args.mapped_reads):
        parser.error("error in required argument mapped_reads: file %s cannot be found. " % (args.mapped_reads))


    bamFileName, bamFileExtension = os.path.splitext(args.mapped_reads)
    if bamFileExtension != ".bam":
        parser.error("the file specified for argument mapped_reads does not end in .bam, are you sure it's a BAM file?")

    #print bamFileName, bamFileExtension
    if not os.path.exists(args.mapped_reads + ".bai"):
        parser.error("it does not seem that the bam file specified for argument mapped_reads is indexed. Please index using samtools index.")

    if not os.path.exists(args.mapped_reads + ".bai"):
        parser.error("error in required argument mapped_reads: file %s cannot be found. " % (args.mapped_reads))
 

    if not os.path.exists(args.TE_annot):
        parser.error("error in required argument TE_annot: file %s cannot be found. " % (args.TE_annot))


    if os.path.dirname(args.output_prefix) and not os.path.exists(os.path.dirname(args.output_prefix)):
        parser.error("error in optional argument --output_prefix: directory %s does not exist" % (args.output_prefix))

    if args.disc_reads_bam and not os.path.exists(args.disc_reads_bam):
        parser.error("error in optional argument --disc_reads_bam: file %s does not exist" % (args.disc_reads_bam))

    if args.conf_lib_stats and not os.path.exists(args.conf_lib_stats):
        parser.error("error in optional argument --conf_lib_stats: file %s does not exist" % (args.conf_lib_stats))

    # if args.pre_filter and not (args.output_prefix and args.conf_lib_stats):
    #     parser.error("error: if --pre_filter is set, --output_prefix and --conf_lib_stats must be too")

    
    ##set appropriate booleans. This is slightly redundant and could be collapsed, but for legacy reasons it's this way
    #TODO: streamline this
    parallel = args.numCPUs >= 1
    already_calc_discordant_reads = args.disc_reads_bam != None

    ##legacy vars from previous versions, that are not worth putting in public opts but I don't want to get rid of quite yet

    #option of re-mapping to a set of TE (or other) sequences. This could be useful in the future for looking at viral insertions, for example
    te_seqs = None

    #to generate a test bam file of input bam for testing: extracts read pairs where either read1 or read2 map to chr1 or chr2
    generate_test_bam = False

    mem = False

    if args.pre_filter:
        run_jitterbug_streaming(args.mapped_reads, args.verbose, args.TE_annot, te_seqs, \
            args.lib_name, args.sdev_mult, args.output_prefix, args.TE_name_tag, parallel, \
            args.numCPUs, args.bin_size, args.minMAPQ, generate_test_bam, args.pre_filter, args.conf_lib_stats, mem, args.min_cluster_size,args.step_one_only,args.step_two_only)

    else:
        run_jitterbug(args.mapped_reads, already_calc_discordant_reads, \
            args.disc_reads_bam, args.verbose, args.TE_annot, te_seqs, \
            args.lib_name, args.sdev_mult, args.output_prefix, args.TE_name_tag, parallel, \
            args.numCPUs, args.bin_size, args.minMAPQ, generate_test_bam, args.pre_filter, args.conf_lib_stats, mem, args.min_cluster_size)


if __name__ == "__main__":
    main(sys.argv[1:])
