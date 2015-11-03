#!/software/so/el6.3/PythonPackages-2.7.3-VirtualEnv/bin/python

import sys
import getopt
#from GffAnnot import *
import math
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import *
from pylab import *


def main(argv):

    # get options passed at command line
    try:
        opts, args = getopt.getopt(argv, "g:G:c:o:t:")
    except getopt.GetoptError:
        #print helpString
        sys.exit(2)


    tag_list = None
    gff_input_2 = None
    gff_file_2=None
    tag_conf_file = False
    gff_file_2 = None
    label_2 = None
    #print opts

    for opt, arg in opts:
        if opt == '-g':
            gff_input_1 = arg
        elif opt == '-G':
            gff_input_2 = arg
        elif opt == '-c':
            tag_conf_file = arg
        elif opt == '-o':
            outfile_name = arg
        elif opt == '-t':
            tag_conf_file = arg


    

    gff_file_1, label_1 = gff_input_1.split(",")
    

    if gff_input_2:
        gff_file_2, label_2 = gff_input_2.split(",")
        

    plot(gff_file_1, label_1, gff_file_2, label_2, tag_conf_file, outfile_name)


def plot(gff_file_1, label_1, gff_file_2, label_2, tag_conf_file, outfile_name):

    

    IN_GFF_FILE_1 = open(gff_file_1, "r")

    feature_list_1 = []


    for line in IN_GFF_FILE_1:
        line = line.strip()
        #ignore empty line or comment lines
        if line[0] == '#':
            continue
        if line == "":
            continue
        feature = GffAnnot(line)
        feature_list_1.append(feature)




    if tag_conf_file:
        tag_keys_conf = open(tag_conf_file)
        tag_conf_list = []
        for line in tag_keys_conf:
            #print line
            if line == "":
                continue
            if line[0] == '#':
                continue
            (tag, xmax, bin_size, ymax, tick_step ) = line.split("\t")
            if float(bin_size) < 1:
                bin_size = float(bin_size)
            else:
                bin_size = int(bin_size)
            tag_conf_list.append((tag, int(xmax), bin_size, int(ymax), int(tick_step) ))
            tot_plots = len(tag_conf_list)

    else:
        tag_keys = []
        all_tag_keys = feature_list_1[0].tags.keys()
        #select only keys that have values convertible to float
        for key in all_tag_keys:
            #print key
            try:
                key_val = float(feature_list_1[0].tags[key])
                tag_keys.append(key)
            except ValueError:
                if key == "softclipped_pos":
                    tag_keys.append(key)
                else:
                    continue
        tot_plots = len(tag_keys) + 1

    num_rows = math.ceil(tot_plots / 2.0)


    figure = plt.figure(figsize=(7, 6), dpi=300)

    # set font aspects
    font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 6}

    rc('font', **font)  # pass in the font dict as kwargs


    sizes1 = [feature.stop - feature.start for feature in feature_list_1]

    if gff_file_2:
        IN_GFF_FILE_2 = open(gff_file_2, "r")
        feature_list_2 = []

        for line in IN_GFF_FILE_2:
            line = line.strip()
            #ignore empty line or comment lines
            if line[0] == '#':
                continue
            if line == "":
                continue
            feature = GffAnnot(line)
            feature_list_2.append(feature)

        sizes2 = [feature.stop - feature.start for feature in feature_list_2]
        print len(sizes2)

       ############################################ CONF FILE PROVIDED ##############################################################
      #
    if tag_conf_file:


        plot_count = 1


        for (tag, set_xmax, set_bin_size, set_ymax, set_tick_step ) in tag_conf_list:
            if tag == "softclipped_pos":
                key_vals_1 = [float(get_softclipped_support_val(feature.tags[tag])) for feature in feature_list_1]

                if gff_file_2:
                    key_vals_2 = [float(get_softclipped_support_val(feature.tags[tag])) for feature in feature_list_2]

            elif tag == "length":
                key_vals_1 = sizes1
                if gff_file_2:
                        key_vals_2 = sizes2

            else:
                try:
                    key_vals_1 = [float(feature.tags[tag]) if float(feature.tags[tag]) < set_xmax else set_xmax for feature in feature_list_1]
                    if gff_file_2:
                        key_vals_2 = [float(feature.tags[tag]) if float(feature.tags[tag]) < set_xmax else set_xmax for feature in feature_list_2]
                except (ValueError, KeyError):
                    continue
            key_plot = figure.add_subplot(num_rows, 3, plot_count)

            # if set_xmax > max(key_vals_1):
            #     set_xmax = max(key_vals_1)
            hist_bins = arange(0, set_xmax + set_bin_size, set_bin_size)
            hist_bins = list(hist_bins)
            hist_bins.append(hist_bins[-1] + set_bin_size)
            #print hist_bins


            #PRINT GFF FILE 1
            hist, edges = np.histogram(key_vals_1, hist_bins)
            # print hist

            hist = list(hist)

            range_all = np.arange(len(hist))

            key_plot.bar(range_all[:-1], hist[:-1], width=0.5, lw=0, color = "grey", alpha = 0.5, label = label_1)
            key_plot.bar(range_all[-1], hist[-1], width=0.5, lw=0, color = "grey", alpha = 0.5)

            # PRITN GFF FILE 2
            if gff_file_2 and len(sizes2) > 0:
                hist2, edges2 = np.histogram(key_vals_2, hist_bins)
                #print hist2
                hist2 = list(hist2)

                range_all = np.arange(len(hist2))
                shift_range = [pos + 0.5 for pos in range_all]
                key_plot.bar(shift_range[:-1], hist2[:-1], width=0.5, lw=0, color = "red",  alpha = 0.5, label = label_2)
                key_plot.bar(shift_range[-1], hist2[-1], width=0.5, lw=0, color = "red",  alpha = 0.5)

            key_plot.set_xlabel(tag)
            key_plot.set_ylabel("count")



            ticks = range(0, set_ymax + set_tick_step, set_tick_step)
            yticks(ticks)

            hist_bins[-2] = ">%d" % hist_bins[-2]
            xticks(range_all, hist_bins, rotation=45)

            plot_count += 1
        ############################################ CONF FILE NOT PROVIDED ##############################################################
    else:
        size_plot = figure.add_subplot(num_rows, 3, 1 )

        n, bins, patches = size_plot.hist(sizes1,  bins=10, color="blue", alpha=0.5, rwidth=0.3, align="left", label=label_1 )
        if gff_file_2 and len(sizes2) > 0:
            n, bins, patches = size_plot.hist(sizes2,  bins=bins, color="red", alpha=0.5, rwidth=0.3, align="mid", label=label_2)

        size_plot.set_xlabel("length")
        size_plot.set_ylabel("count")
        legend()
        plot_count = 2


        for key in tag_keys:
            if key == "softclipped_pos":
                continue
                key_vals_1 = [float(get_softclipped_support_val(feature.tags[key])) for feature in feature_list_1]

                if gff_file_2 and len(sizes2) > 0:
                    key_vals_2 = [float(get_softclipped_support_val(feature.tags[key])) for feature in feature_list_2]

            else:
                try:
                    key_vals_1 = [float(feature.tags[key]) for feature in feature_list_1]
                    if gff_file_2 and len(sizes2) > 0:
                        key_vals_2 = [float(feature.tags[key]) for feature in feature_list_2]
                except (ValueError, KeyError):
                    continue
            key_plot = figure.add_subplot(num_rows, 3, plot_count)

            n, bins, patches = key_plot.hist(key_vals_1, bins=20, facecolor="blue", alpha=0.5, rwidth=0.3, align="left", label = label_1)

            if gff_file_2 and len(sizes2) > 0:
                n, bins, patches = key_plot.hist(key_vals_2, bins=bins, facecolor="red", alpha=0.5, rwidth=0.3, align="mid", label = label_2)

            key_plot.set_xlabel(key)
            key_plot.set_ylabel("count")

            plot_count += 1



    plt.tight_layout()
    legend(loc="upper left", bbox_to_anchor=(0.45,0.85))


    #plt.show()
    if not outfile_name:
        outfile_name = "%s.metrics.pdf" % (gff_file_1)
        if gff_file_2:
            outfile_name = "%s.metrics.pdf" % (gff_file_1)
    plt.savefig(outfile_name, format='pdf')

def get_softclipped_support_val(tag_value):
    values = tag_value[1:-1].split(",")
    return values[2]

class GffAnnot:


    def __init__(self, line):
                    (chrom, so_type, method, start, stop, score, strand, frame, tags) = line.split("\t")
                    self.chrom = chrom
                    self.type = so_type
                    self.method = method
                    self.start = int(start) - 1
                    self.stop = int(stop)
                    self.score = score
                    self.strand = strand
                    self.frame = frame
                    self.tags = {}
                    self.gff_line = line
                    self.length = (self.stop - self.start) + 1

                    #print line
                    tag_pairs = tags.split(";")
                    for pair in tag_pairs:
                    #                print pair
                                    pair.strip()
                                    if pair == "":
                                        continue
                                    if pair == ".":
                                        continue
                                    #print pair
                                    (key, value) = pair.split("=")
                                    key = key.strip()
                                    value = value.strip()
                                    self.tags[key] = value

    def __eq__(self, other):
        if self.chrom != other.chrom:
            return False
        if self.type != other.type:
            return False
        if self.method != other.method:
            return False
        if self.start != other.start:
            return False
        if self.stop != other.stop:
            return False
        if self.score != other.score:
            return False
        if self.strand != other.strand:
            return False
        if self.frame != other.frame:
            return False
        for (key, value) in self.tags.iteritems():
            if key not in other.tags:
                return False
            elif other.tags[key] != value:
                return False
        return True

    def to_string(self):
        tag_string = ""
        for tag, value in enumerate(self.tags):
            tag_string += "=".join(tag, value)
        return "\t".join(self.chrom, self.type, self.method, self.start, self.end, self.score, self.strand, self.frame, tag_string)


if __name__ == "__main__":
    main(sys.argv[1:])



















