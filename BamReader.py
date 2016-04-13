import pysam
import pybedtools
import sys
import numpy
from AlignedReadPair import *
import gc
#from memory_profiler import profile
import y_serial_v060 as y_serial
import cPickle as pickle


class BamReader:

    def __init__(self, bam_file_name, prefix):
        self.bam_file_name = bam_file_name
        if prefix:
            self.prefix = prefix
        else:
            self.prefix = bam_file_name


    def calculate_mean_sdev_isize(self, num_itr):
        """calculate the mean insert size and std dev for num_itr first reads in the bam file"""
        bam_file = pysam.Samfile(self.bam_file_name, "rb")
        # itr = bam_file.fetch(bam_file.references[0], start=50000)
        print "blah"
        # print itr.tell()
        counter  = 0
        isize_array = []
        read_length_array = []

        while counter < num_itr - 1:
            read = bam_file.next()
            #print str(read)
            #print read.fancy_str()
            if read.is_proper_pair and read.mapq > 30 and read.isize >0:
                # print read.pos
                isize_array.append(abs(read.isize))
                read_length_array.append(read.rlen)
                counter += 1
        mean = numpy.median(isize_array)
        sdev = numpy.std(isize_array)
        rlen_mean = numpy.mean(read_length_array)
        rlen_sdev = numpy.std(read_length_array)
        bam_file.close()
        return(mean, sdev, rlen_mean, rlen_sdev)


    def output_repetitive_reads(self):
        bam_file = pysam.Samfile(self.bam_file_name, "rb")
        output_bam_file = pysam.Samfile(self.prefix + ".repetitive.bam", mode="wb", referencenames=bam_file.references, referencelengths=bam_file.lengths)

        read = bam_file.next()
        while 1:

            if is_mapped_mult_times(read):
                output_bam_file.write(read)

            try:
                read = bam_file.next()
            except StopIteration:
                break
        bam_file.close()
        output_bam_file.close()

    def output_one_chr_reads(self):
        bam_file = pysam.Samfile(self.bam_file_name, "rb")
        output_bam_file = pysam.Samfile(self.prefix + ".chr1_chr2.bam", mode="wb", referencenames=bam_file.references, referencelengths=bam_file.lengths)
        read1 = bam_file.next()

        read_pairs_dict = {}



        while 1 :
            #if verbose:
            #    print read1
            #    print read2


            if read_pairs_dict.has_key(read1.qname):
                if read1.is_read1:
                    read_pairs_dict[read1.qname][0] = read1
                else:
                    read_pairs_dict[read1.qname][1] = read1

            elif bam_file.getrname(read1.rname) == "Chr1" or bam_file.getrname(read1.rname) == "Chr2" :
                
                if read1.is_read1:
                    read_pairs_dict[read1.qname] = [read1,None]
                else:
                    read_pairs_dict[read1.qname] = [None, read1]


                #shift to next pair
            try:
                read1 = bam_file.next()
            except StopIteration:
                break

        for [read1, read2] in read_pairs_dict.itervalues():
            bam_file.write(read1)
            bam_file.write(read2)
        bam_file.close()
    #@profile 
    def select_read_pair_one_overlap_TE_annot(self, TE_annot, int_size, min_mapq,db,bin_size=50000000):
	
        """ output bam file of read pairs where exactly one read overlaps with an annotation in supplied gff file\
            also returns a BedTools object """

       # print "selecting discordant reads that overlap with a TE in annotation " + TE_annot + " ..."

        #use pysam to open the bam file because it has better object definition for the reads
        valid_discordant_bam = pysam.Samfile(self.bam_file_name, "rb")

        #file to save the discordant read pairs where exactly one read overlaps a TE annotation
        #overlap_TE_bam_file = pysam.Samfile(self.prefix + ".one_read_overlap_TE.bam", mode="wb", referencenames=valid_discordant_bam.references, referencelengths=valid_discordant_bam.lengths)

        #use pybedtools to look up the overlap of the Interval defined by the read with the Intervals defined by the gff file
        TE_annot_intervals = pybedtools.IntervalFile(TE_annot)


        #make a list of AlignedReadPair objects for each read pair in the list that has exactly one read overlapping a TE
        read_pairs_xor_overlap_TE = []
	print db
	read_pair_database = y_serial.Main( db )
	bin_list=list()

        try:
        	read1 = valid_discordant_bam.next()
        	read2 = valid_discordant_bam.next()
        except StopIteration:
            print "ERROR: no reads are found in %s, exiting" % (bam_file_name)
            sys.exit(2)



        while 1 :
            #if verbose:
            #    print read1
            #    print read2

            #check that the reads are truly a pair:
            #if not, scoot down one in the iteration

            if read1.qname != read2.qname:
                print "unmatched pair in valid discordant reads. Problem!!"
                #sys.exit(2)
                read1 = read2
                try:
                    read2 = valid_discordant_bam.next()
                except StopIteration:
                    break
                continue

            read_pair = AlignedReadPair(read1, read2)


            #see if read1 is a TE
            read1_all_mappings = get_all_mapping_pos(read1, valid_discordant_bam)
            for (chr, start, end) in read1_all_mappings:
                map_interval = pybedtools.Interval(chr, start, end, strand='+')
                overlapping_TE_annots = TE_annot_intervals.all_hits(map_interval)
                if len(overlapping_TE_annots) > 0:
                    read_pair.TE_annot_attr_list.extend([ gff_interval.attrs for gff_interval in overlapping_TE_annots ])
                    read_pair.TE_map_gff_list.extend([ str(gff_interval) for gff_interval in overlapping_TE_annots ])
                    read_pair.read1_is_TE = True
                    #print "read1 TE"
                    #if read1 is TE, then read2 is the anchor, so set the interval chr to the chr of read2
                    read_pair.interval_chr = valid_discordant_bam.getrname(read2.rname)


            #see if read2 is a TE
            read2_all_mappings = get_all_mapping_pos(read2, valid_discordant_bam)
            for (chr, start, end) in read2_all_mappings:
                map_interval = pybedtools.Interval(chr, start, end, strand='+')
                overlapping_TE_annots = TE_annot_intervals.all_hits(map_interval)
                if len(overlapping_TE_annots) > 0:
                    read_pair.TE_annot_attr_list.extend([ gff_interval.attrs for gff_interval in overlapping_TE_annots ])
                    read_pair.TE_map_gff_list.extend([ str(gff_interval) for gff_interval in overlapping_TE_annots ])
                    read_pair.read2_is_TE = True
                    #print "read2 TE"
                    #if read2 is TE, then read1 is the anchor, so set the interval chr to the chr of read1
                    read_pair.interval_chr = valid_discordant_bam.getrname(read1.rname)


            #only add the AlignedRead to the list if exactly one read maps to a TE location, and the anchor is not repetitive

            if read_pair.read1_is_TE and not read_pair.read2_is_TE and not is_mapped_mult_times(read2):
                if min_mapq:
                    if read2.mapq >= min_mapq:

                        read_pair.calculate_outside_interval(int_size, read1, read2)
                        read_pair.calc_anchor_is_softclipped(read1, read2)

                        ###### TODO: this is where the Aligned ReadPair objects are deleted

                        #read_pair.read1 = None
                        #read_pair.read2 = None

                        #print read_pair.read1

                        #############
                        #read_pairs_xor_overlap_TE.append(read_pair)
			tabname="[c%s_%d_%d_%s]"%(read_pair.interval_chr,bin_size*(int(read_pair.interval_start)/bin_size),bin_size*(1+int(read_pair.interval_start)/bin_size) ,read_pair.interval_direction  )
			#read_pair_database.insert(read_pair,'',tabname)
			read_pairs_xor_overlap_TE.append([read_pair,tabname])
                        ##overlap_TE_bam_file.write(read_pair.read1)
                        ##overlap_TE_bam_file.write(read_pair.read2)
                else:
                    read_pair.calculate_outside_interval(int_size, read1, read2)
                    read_pair.calc_anchor_is_softclipped(read1, read2)
                    #read_pairs_xor_overlap_TE.append(read_pair)
		    tabname="[c%s_%d_%d_%s]"%(read_pair.interval_chr,bin_size*(int(read_pair.interval_start)/bin_size),bin_size*(1+int(read_pair.interval_start)/bin_size) ,read_pair.interval_direction  )
		    #read_pair_database.insert(read_pair,'',tabname)
		    read_pairs_xor_overlap_TE.append([read_pair,tabname])

            elif read_pair.read2_is_TE and not read_pair.read1_is_TE and not is_mapped_mult_times(read1):
                if min_mapq:
                    if read1.mapq >= min_mapq:

                        read_pair.calculate_outside_interval(int_size, read1, read2)
                        read_pair.calc_anchor_is_softclipped(read1, read2)

                        #read_pairs_xor_overlap_TE.append(read_pair)
			tabname="[c%s_%d_%d_%s]"%(read_pair.interval_chr,bin_size*(int(read_pair.interval_start)/bin_size),bin_size*(1+int(read_pair.interval_start)/bin_size) ,read_pair.interval_direction  )
			
			#read_pair_database.insert(read_pair,'',tabname)
			read_pairs_xor_overlap_TE.append([read_pair,tabname])
                        ##overlap_TE_bam_file.write(read_pair.read1)
                        ##overlap_TE_bam_file.write(read_pair.read2)
                else:
                    read_pair.calculate_outside_interval(int_size, read1, read2)
                    read_pair.calc_anchor_is_softclipped(read1, read2)

                    #read_pairs_xor_overlap_TE.append(read_pair)
		    tabname="[c%s_%d_%d_%s]"%(read_pair.interval_chr,bin_size*(int(read_pair.interval_start)/bin_size),bin_size*(1+int(read_pair.interval_start)/bin_size) ,read_pair.interval_direction  )
		    #read_pair_database.insert(read_pair,'',tabname)
		    read_pairs_xor_overlap_TE.append([read_pair,tabname])
	    if len(read_pairs_xor_overlap_TE) >= 100000:
		read_pair_database.ingenerator(read_pairs_xor_overlap_TE,'read_pairs')
		#try to get the unique
		tmp = list()
		tmp = [bin_ for  pair,bin_ in read_pairs_xor_overlap_TE]
		bin_list.extend(list(set(tmp)))		
		del read_pairs_xor_overlap_TE
		read_pairs_xor_overlap_TE=list()
	    
            #shift to next pair
            try:
                read1 = valid_discordant_bam.next()
                read2 = valid_discordant_bam.next()
            except StopIteration:
                break
	if len(read_pairs_xor_overlap_TE) >0:
		read_pair_database.ingenerator(read_pairs_xor_overlap_TE,'read_pairs')
		#try to get the unique
		tmp = list()
		tmp = [bin_ for  pair,bin_ in read_pairs_xor_overlap_TE]
		bin_list.extend(list(set(tmp)))
		del read_pairs_xor_overlap_TE
		read_pairs_xor_overlap_TE=list()
	bin_list = list(set(bin_list))
	read_pair_database.insert(bin_list,'','bin_list')
	print bin_list
        print "number discordant read pairs with exactly one read overlapping a TE: %d" % len(read_pairs_xor_overlap_TE)
        #print "\n".join(pair.str() for pair in read_pairs_xor_overlap_TE)
        #overlap_TE_bam_file.close()

        valid_discordant_bam.close()
        return read_pairs_xor_overlap_TE





    def select_discordant_reads(self, strict_repetitive, verbose, isize, outfile_name):
        """ This function selects discordant read pairs from a bam file that are putatively predictive of a transposable element insertion/deletion in the resequenced sample\
        valid discordant read pairs are those that:\
        have a mapping distance greater that the norm expected (as calculated by bwa)or are mapped to two different chromosomes, \
        AND \
        have at least one uniquely mapping read per pair\
        these valid discordant reads can be further separated between those that have two uniquely mapping reads \
        and those that have exactly one uniquely and one repetitively mapping read using the -s option\
        parameters: \
        -strict_repetitive if set will force to separate valid discordant read pairs between those that have two uniquely mapping reads (output to <sorted_bam_file>.valid_discordant_pairs.bam) \
        and those that have exactly one uniquely mapping read and one repetitively mapping read (output to <sorted_bam_file>.valid_discordant_pairs_strict_rep.bam) \
        -verbose, print to std out each read pair and how it is categorized (proper mapped, discordant too small insert size, valid discordant both unique, valid discordant unique/rep, unmapped)"""


       # print "selecting discordant reads..."

        #open file in r (read) b (bam) mode
        bam_file = pysam.Samfile(self.bam_file_name, "rb")


        #file to save the softclipped reads that are NOT multiple maps for both sides
        #soft_clipped_bam_file = pysam.Samfile(self.prefix + ".softclipped.bam", mode="wb", referencenames=bam_file.references, referencelengths=bam_file.lengths)
        #proper_pair_bam_file = pysam.Samfile(self.prefix + ".proper_pair.bam", mode="wb", referencenames=bam_file.references, referencelengths=bam_file.lengths)

        #file to save the valid discordant pairs: with at least one uniquely mapping read and all insert sizes are greater than expected.
        valid_discordant_pairs = pysam.Samfile(outfile_name, mode="wb", referencenames=bam_file.references, referencelengths=bam_file.lengths)


        #keep track of what kind of reads were found
        single_read_count = 0
        multiple_map_pair_count = 0
        discordant_pair_too_small_isize = 0
        valid_discordant_pair_count = 0
        proper_pair_count = 0
        unmapped_pair_count = 0
        valid_discordant_pair_count_strict = 0
        #total_reads = bam_file.mapped + bam_file.unmapped

        read1 = bam_file.next()
        read2 = bam_file.next()

        #for testing restrict while loop to counter
        #counter = 0
        #print "TESTING!! only first 100000 reads analyzed! "
        #while counter < 100000:
        #    counter += 1
        while 1 :
            if verbose:
                print "next pair"

            if verbose:
                print read1
                print read2
            #check that the reads are truly a pair:
            #if not, scoot down one in the iteration

            if read1.qname != read2.qname:
                if verbose:
                    print "unmatched pair! "
                single_read_count += 1
               # if is_softclipped(read1):
                   # soft_clipped_bam_file.write(read1)
                read1 = read2
                try:
                    read2 = bam_file.next()
                except StopIteration:
                    break
                continue

            #do not keep pair if it is properly mapped but save it if either read is softclipped
            if read1.is_proper_pair:
                if verbose:
                    print "proper pair!"
                proper_pair_count +=1
               # if is_softclipped(read1) or is_softclipped(read2):
               #     soft_clipped_bam_file.write(read1)
               #     soft_clipped_bam_file.write(read2)
               # else:
               #     proper_pair_bam_file.write(read1)
               #     proper_pair_bam_file.write(read2)
               #  proper_pair_bam_file.write(read1)
               #  proper_pair_bam_file.write(read2)

                try:
                    read1 = bam_file.next()
                    read2 = bam_file.next()
                except StopIteration:
                    break
                continue


            if read1.is_unmapped or read2.is_unmapped:
                if verbose:
                    print "unmapped!"
                unmapped_pair_count += 1
                try:
                    read1 = bam_file.next()
                    read2 = bam_file.next()
                except StopIteration:
                    break
                continue

           # if min_mapq:
           #     if read1.mapq < min_mapq or read2.mapq < min_mapq:
           #         if verbose:
           #             print "qual less than %d! " % (min_mapq)
           #         unmapped_pair_count += 1
           #         try:
           #             read1 = bam_file.next()
           #             read2 = bam_file.next()
           #         except StopIteration:
           #             break
           #         continue

            #note: isize is deprecated, change to tlen when new version installed
            # isize is 0 when mapped to a different chromosome, so want to keep those
            if 0 < read1.isize < abs(isize):
                if verbose:
                    print "invalid discordant, insert size too small! isize = %d" % read1.isize
                discordant_pair_too_small_isize += 1
               # if is_softclipped(read1) or is_softclipped(read2):
               #     soft_clipped_bam_file.write(read1)
               #     soft_clipped_bam_file.write(read2)
                try:
                    read1 = bam_file.next()
                    read2 = bam_file.next()
                except StopIteration:
                    break
                continue



            #throw away pairs where both map multiple times


            if is_mapped_mult_times(read1) and is_mapped_mult_times(read2):
                if verbose:
                    print "invalid discordant, both reads repetitively map!"
                multiple_map_pair_count += 1
               # if is_softclipped(read1) or is_softclipped(read2):
               #     soft_clipped_bam_file.write(read1)
               #     soft_clipped_bam_file.write(read2)
                try:
                    read1 = bam_file.next()
                    read2 = bam_file.next()
                except StopIteration:
                    break
                continue

            #if strict: separate valid discordant reads between those that have two uniquely mapping reads and and those with one repetitive one
            if strict_repetitive:
                if is_mapped_mult_times(read1) or is_mapped_mult_times(read2):
                    if verbose:
                        print "discordant and exactly one read repetitive! map distance = %d" % read1.isize
                    valid_discordant_pair_count_strict += 1
                    valid_discordant_pairs_strict.write(read1)
                    valid_discordant_pairs_strict.write(read2)
                   # if is_softclipped(read1) or is_softclipped(read2):
                   #     soft_clipped_bam_file.write(read1)
                   #     soft_clipped_bam_file.write(read2)

                    try:
                        read1 = bam_file.next()
                        read2 = bam_file.next()
                    except StopIteration:
                        break
                else:
                    if verbose:
                        print "discordant and both reads map uniquely! map distance = %d" % read1.isize
                    valid_discordant_pair_count += 1
                    valid_discordant_pairs.write(read1)
                    valid_discordant_pairs.write(read2)
                   # if is_softclipped(read1) or is_softclipped(read2):
                   #     soft_clipped_bam_file.write(read1)
                   #     soft_clipped_bam_file.write(read2)

                    try:
                        read1 = bam_file.next()
                        read2 = bam_file.next()
                    except StopIteration:
                        break


            #otherwise, output to same file whether both unique or one unique
            else:
                if verbose:
                    print "discordant! map distance = %d" % read1.isize
                valid_discordant_pair_count += 1
                valid_discordant_pairs.write(read1)
                valid_discordant_pairs.write(read2)
               # if is_softclipped(read1) or is_softclipped(read2):
               #     soft_clipped_bam_file.write(read1)
               #     soft_clipped_bam_file.write(read2)

                try:
                    read1 = bam_file.next()
                    read2 = bam_file.next()
                except StopIteration:
                    break

        bam_stats = {}
        bam_stats["single_read_count"] = single_read_count
        bam_stats["multiple_map_pair_count"] = multiple_map_pair_count
        bam_stats["valid_discordant_pair_count"] = valid_discordant_pair_count
        bam_stats["proper_pair_count"] = proper_pair_count
        bam_stats["unmapped_pair_count"] = unmapped_pair_count
        bam_stats["discordant_pair_too_small_insert_size"] = discordant_pair_too_small_isize
        bam_stats["valid_discordant_pair_count_strict"] = valid_discordant_pair_count_strict
        bam_stats["total_reads"] = single_read_count + 2*(multiple_map_pair_count + valid_discordant_pair_count + proper_pair_count + unmapped_pair_count + discordant_pair_too_small_isize + valid_discordant_pair_count_strict)

        #proper_pair_bam_file.close()
        valid_discordant_pairs.close()
        if strict_repetitive:
            valid_discordant_pairs_strict.close()

        # print "done selecting discordant reads."
        lengths = bam_file.lengths
        refs = bam_file.references
        bam_file.close()

        return (bam_stats, lengths, refs)


    def select_discordant_reads_psorted(self, verbose, isize, outfile_name):

        """ This function selects discordant read pairs from a bam file that are putatively predictive of a transposable element insertion/deletion in the resequenced sample\
        valid discordant read pairs are those that:\
        have a mapping distance greater that the norm expected (as calculated by bwa)or are mapped to two different chromosomes, \
        AND \
        have at least one uniquely mapping read per pair\
        these valid discordant reads can be further separated between those that have two uniquely mapping reads \
        and those that have exactly one uniquely and one repetitively mapping read using the -s option\
        parameters: \
        -strict_repetitive if set will force to separate valid discordant read pairs between those that have two uniquely mapping reads (output to <sorted_bam_file>.valid_discordant_pairs.bam) \
        and those that have exactly one uniquely mapping read and one repetitively mapping read (output to <sorted_bam_file>.valid_discordant_pairs_strict_rep.bam) \
        -verbose, print to std out each read pair and how it is categorized (proper mapped, discordant too small insert size, valid discordant both unique, valid discordant unique/rep, unmapped)"""


       # print "selecting discordant reads..."

        #open file in r (read) b (bam) mode
        bam_file = pysam.Samfile(self.bam_file_name, "rb")


        #file to save the softclipped reads that are NOT multiple maps for both sides
        #soft_clipped_bam_file = pysam.Samfile(self.prefix + ".softclipped.bam", mode="wb", referencenames=bam_file.references, referencelengths=bam_file.lengths)
        #proper_pair_bam_file = pysam.Samfile(self.prefix + ".proper_pair.bam", mode="wb", referencenames=bam_file.references, referencelengths=bam_file.lengths)

        #file to save the valid discordant pairs: with at least one uniquely mapping read and all insert sizes are greater than expected.
        valid_discordant_pairs = pysam.Samfile(outfile_name, mode="wb", referencenames=bam_file.references, referencelengths=bam_file.lengths)

        read_names_set = set([])
        read_pairs_dict = {}

        #keep track of what kind of reads were found
        single_read_count = 0
        multiple_map_pair_count = 0
        discordant_pair_too_small_isize = 0
        valid_discordant_pair_count = 0
        proper_pair_count = 0
        unmapped_pair_count = 0
        valid_discordant_pair_count_strict = 0
        total_reads_count = 0

        #total_reads = bam_file.mapped + bam_file.unmapped

        read1 = bam_file.next()
        

        #for testing restrict while loop to counter
        # counter = 0
        # print "TESTING!! only first 100000 reads analyzed! "
        # while counter < 1000:
        #     counter += 1
        while 1 :

            total_reads_count += 1
            # print read1.qname

            if verbose:
                print "next pair"

            if verbose:
                print read1

            

            #do not keep if it is properly mapped 
            if read1.is_proper_pair:
                if verbose:
                    print "proper pair!"
                proper_pair_count +=1
                try:
                    read1 = bam_file.next()
                except StopIteration:
                    break
                continue


            if read1.is_unmapped:
                if verbose:
                    print "unmapped!"
                unmapped_pair_count += 1
                try:
                    read1 = bam_file.next()
                except StopIteration:
                    break
                continue


            #note: isize is deprecated, change to tlen when new version installed
            # isize is 0 when mapped to a different chromosome, so want to keep those
            if 0 < read1.isize < abs(isize):
                if verbose:
                    print "invalid discordant, insert size too small! isize = %d" % read1.isize
                discordant_pair_too_small_isize += 1
                try:
                    read1 = bam_file.next()
                except StopIteration:
                    break
                continue



            
            if verbose:
                print "discordant! map distance = %d" % read1.isize
            valid_discordant_pair_count += 1

            # if this read is mate of one that was already selected, write it
            if read1.qname in read_pairs_dict:
                if read1.is_read1:
                    valid_discordant_pairs.write(read1)
                    valid_discordant_pairs.write(read_pairs_dict[read1.qname])
		    del read_pairs_dict[read1.qname]
                else:
                    valid_discordant_pairs.write(read_pairs_dict[read1.qname])
                    valid_discordant_pairs.write(read1)
		    del read_pairs_dict[read1.qname]
                    

            # otherwise, make a new entry in the dictionary
            else:
                
                read_pairs_dict[read1.qname] = read1
                
                #read_names_set.add(read1.qname)
            


            try:
                read1 = bam_file.next()
            except StopIteration:
                break

	read_pairs_dict.clear()
	gc.collect()
        bam_stats = {}
        bam_stats["single_read_count"] = single_read_count
        bam_stats["valid_discordant_pair_count"] = valid_discordant_pair_count
        bam_stats["proper_pair_count"] = proper_pair_count
        bam_stats["unmapped_pair_count"] = unmapped_pair_count
        bam_stats["discordant_pair_too_small_insert_size"] = discordant_pair_too_small_isize
        bam_stats["total_reads"] = total_reads_count


        # print "writing discordant reads ..."
        # for reads in read_pairs_dict.itervalues():
        #     print reads
        #     valid_discordant_pairs.write(reads[0])
        #     valid_discordant_pairs.write(reads[1])

        valid_discordant_pairs.close()


        print "done selecting discordant reads."
        lengths = bam_file.lengths
        refs = bam_file.references
        bam_file.close()

	del read_names_set
	del read_pairs_dict
	gc.collect()
        return (bam_stats, lengths, refs)

    


    def select_discordant_reads_psorted2(self, verbose, isize, outfile_name):
    
	    """ This function selects discordant read pairs from a bam file that are putatively predictive of a transposable element insertion/deletion in the resequenced sample\
	    valid discordant read pairs are those that:\
	    have a mapping distance greater that the norm expected (as calculated by bwa)or are mapped to two different chromosomes, \
	    AND \
	    have at least one uniquely mapping read per pair\
	    these valid discordant reads can be further separated between those that have two uniquely mapping reads \
	    and those that have exactly one uniquely and one repetitively mapping read using the -s option\
	    parameters: \
	    -strict_repetitive if set will force to separate valid discordant read pairs between those that have two uniquely mapping reads (output to <sorted_bam_file>.valid_discordant_pairs.bam) \
	    and those that have exactly one uniquely mapping read and one repetitively mapping read (output to <sorted_bam_file>.valid_discordant_pairs_strict_rep.bam) \
	    -verbose, print to std out each read pair and how it is categorized (proper mapped, discordant too small insert size, valid discordant both unique, valid discordant unique/rep, unmapped)"""
    
    
	   # print "selecting discordant reads..."
    
	    #open file in r (read) b (bam) mode
	    bam_file = pysam.Samfile(self.bam_file_name, "rb")
    
    
	    #file to save the softclipped reads that are NOT multiple maps for both sides
	    #soft_clipped_bam_file = pysam.Samfile(self.prefix + ".softclipped.bam", mode="wb", referencenames=bam_file.references, referencelengths=bam_file.lengths)
	    #proper_pair_bam_file = pysam.Samfile(self.prefix + ".proper_pair.bam", mode="wb", referencenames=bam_file.references, referencelengths=bam_file.lengths)
    
	    #file to save the valid discordant pairs: with at least one uniquely mapping read and all insert sizes are greater than expected.
	    valid_discordant_pairs = pysam.Samfile(outfile_name, mode="wb", referencenames=bam_file.references, referencelengths=bam_file.lengths)
    
	    read_names_set = set([])
	    read_pairs_dict = dict()
    
	    #keep track of what kind of reads were found
	    single_read_count = 0
	    multiple_map_pair_count = 0
	    discordant_pair_too_small_isize = 0
	    valid_discordant_pair_count = 0
	    proper_pair_count = 0
	    unmapped_pair_count = 0
	    valid_discordant_pair_count_strict = 0
	    total_reads_count = 0
    
	    #total_reads = bam_file.mapped + bam_file.unmapped
    
	    read1 = bam_file.next()
	    while 1 :
    
		total_reads_count += 1
		#do not keep if it is properly mapped 
		if read1.is_proper_pair:
		    proper_pair_count +=1
		    try:
			read1 = bam_file.next()
		    except StopIteration:
			break
		    continue
    
    
		if read1.is_unmapped:
		    unmapped_pair_count += 1
		    try:
			read1 = bam_file.next()
		    except StopIteration:
			break
		    continue
		#note: isize is deprecated, change to tlen when new version installed
		# isize is 0 when mapped to a different chromosome, so want to keep those
		if 0 < read1.isize < abs(isize):
		    if verbose:
			print "invalid discordant, insert size too small! isize = %d" % read1.isize
		    discordant_pair_too_small_isize += 1
		    try:
			read1 = bam_file.next()
		    except StopIteration:
			break
		    continue
		valid_discordant_pair_count += 1
		# write the read in its corresponding file
		readID = read1.qname
		fileid = readID.split(':')[2]
		if fileid in read_pairs_dict:
		    pickle.dump(read1,read_pairs_dict[fileid],pickle.HIGHEST_PROTOCOL)
		else:
		    read_pairs_dict[fileid]=open(outfile_name+'_'+fileid+'pkl','wb')
		    pickle.dump(read1,read_pairs_dict[fileid],pickle.HIGHEST_PROTOCOL)
		try:
		    read1 = bam_file.next()
		except StopIteration:
		    break
    
    #### Retrieve information for the name-split files and process them to extract pairs
	    for id_ in read_pairs_dict.keys():
		read_pairs_dict[id_].close()
	    print read_pairs_dict.keys()
	    # Now launch 
	    for id_ in read_pairs_dict.keys():
		print id_
		rp_dict = dict()
		with open(outfile_name+'_'+id_+'pkl', "rb") as f:
		    go = True
		    while go:
			try:
			    print 'test'
			    read1 = pickle.load(f)
			    print 'tost'
			    print read1
			    if read1.qname in rp_dict.keys():
				print 'aha'
				valid_discordant_pairs.write(read1)
				valid_discordant_pairs.write(rp_dict[read1.qname])
				del valid_discordant_pairs[read1.qname]
			    else:
			       rp_dict[read1.qname] = read1
			except EOFError:
			    print 'eoferrorno'
			    go = False
		rp_dict.clear()
		#os.unlink(outfile_name+'_'+id_+'pkl')
	    
    
    
    
    #### Final round where stats get collected
	    read_pairs_dict.clear()
	    gc.collect()
	    bam_stats = {}
	    bam_stats["single_read_count"] = single_read_count
	    bam_stats["valid_discordant_pair_count"] = valid_discordant_pair_count
	    bam_stats["proper_pair_count"] = proper_pair_count
	    bam_stats["unmapped_pair_count"] = unmapped_pair_count
	    bam_stats["discordant_pair_too_small_insert_size"] = discordant_pair_too_small_isize
	    bam_stats["total_reads"] = total_reads_count
    
    
	    # print "writing discordant reads ..."
	    # for reads in read_pairs_dict.itervalues():
	    #     print reads
	    #     valid_discordant_pairs.write(reads[0])
	    #     valid_discordant_pairs.write(reads[1])
    
	    valid_discordant_pairs.close()
    
    
	    print "done selecting discordant reads."
	    lengths = bam_file.lengths
	    refs = bam_file.references
	    bam_file.close()
    
	    del read_names_set
	    del read_pairs_dict
	    gc.collect()
	    return (bam_stats, lengths, refs)





## helper functions
def is_softclipped(read):
    if read.cigar == None:
        return 0
    for op, num in read.cigar:
        if op == 4:
            return 1

    return 0




def is_mapped_mult_times(read):
    if read.tags != None:
        try:
            if read.opt('XT') == 'R':
                return 1
        except KeyError:
            pass
        try:
            if read.opt('X0') != None:
                if read.opt('X0') > 1:
                    return 1
        except KeyError:
            pass
        try:
            if read.opt('X1') != None:
                if read.opt('X1') > 0:
                    return 1
        except KeyError:
            pass
        return 0
    else:
        return 0

def get_all_mapping_pos(read, bam_file_obj):
    """ read is an AlignedRead object as defined in pysam. return a list of position intervals corresponding to alternate mappings found in the XA tag\
    note: position is considered to be leftmost one-based, regardless of strand."""
    refname = bam_file_obj.getrname(read.rname)
    #add primary mapping position to list
    positions_list = [(refname, read.pos, read.pos + read.alen)]
    #then add alternate mapping positions
    try:
        alt_locations = read.opt('XA').strip(";").split(";")
        for location in alt_locations:
            (chr, pos, cigar, edit_dist) = location.split(",")
            pos = abs(int(pos))
            positions_list.append( (chr, pos, pos + read.qlen))
    except KeyError:
        return positions_list
    #print len(positions_list)
    return positions_list


