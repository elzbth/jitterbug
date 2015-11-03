#!/usr/bin/python

import sys
import getopt
import pybedtools
import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles



def main(argv):

    # get options passed at command line
    try:
        opts, args = getopt.getopt(argv, "1:2:3:o:")
    except getopt.GetoptError:
        #print helpString
        sys.exit(2)

    gff3 = None

    #print opts
    outfile_name = None

    for opt, arg in opts:
        if opt == '-1':
            gff1, label1 = arg.split(",")
        elif opt == '-2':
            gff2, label2 = arg.split(",")
        elif opt == '-3':
            gff3, label3  = arg.split(",")
        elif opt == '-o':
            outfile_name = arg




def venn_diag(gff1, label1, gff2, label2, gff3, label3, library_name):

    outfile_name = library_name + "_" + label1 + "." + label2
    if gff3:
        outfile_name = outfile_name + "." + label3


    figure = plt.figure()

    a = pybedtools.BedTool(gff1)
    b = pybedtools.BedTool(gff2)
    if gff3:
        c = pybedtools.BedTool(gff3)

        Abc = (a - b - c)
        aBc = (b - a - c)
        ABc = (a + b - c)
        abC = (c - a - b)
        AbC = (a + c - b)
        aBC = (b + c - a)
        ABC = (a + b + c)
        # print "Abc: %d\naBc: %d\nABc:%d\nabC: %d\nAbC: %d\naBC: %d\nABC: %d" % ( Abc.count(), aBc.count(), ABc.count(), abC.count(), AbC.count(), aBC.count(), ABC.count() )

        results = "%s\t%s\t%s\n" % (label1, label2, label3)
        results = results + "x\t \t \t%d\n" % (Abc.count())
        results = results + " \tx\t \t%d\n" % (aBc.count())
        results = results + "x\tx\t \t%d\n" % (ABc.count())
        results = results + " \t \tx\t%d\n" % (abC.count())
        results = results + "x\t \tx\t%d\n" % (AbC.count())
        results = results + " \tx\tx\t%d\n" % (aBC.count())
        results = results + "x\tx\tx\t%d" % (ABC.count())


        v = venn3(subsets=(Abc.count(), aBc.count(), ABc.count(), abC.count(), AbC.count(), aBC.count(), ABC.count()), set_labels = (label1, label2, label3))
        text(0,0,results, transform=ax.transAxes)

        Abc.saveas("%s.not%s.not%s.gff" % (label1, label2, label3))
        aBc.saveas("not%s.%s.not%s.gff" % (label1, label2, label3))
        abC.saveas("not%s.not%s.%s.gff" % (label1, label2, label3))

    else:
        Ab = (a - b)
        aB = (b - a)
        AB = (a + b)

        results = "%s\t%s\n" % (label1, label2)
        results = results + "x\t \t%d\n" % (Ab.count())
        results = results + " \tx\t%d\n" % (aB.count())
        results = results + "x\tx\t%d\n" % (AB.count())

        
        v = venn2(subsets=(Ab.count(), aB.count(), AB.count()), set_labels = (label1, label2))
        #plt.text(0,0,results)

        Ab_name = "%s_%s.not%s.gff" % (library, label1, label2)
        Ab.saveas(Ab_name)

        aB_name = "not%s.%s.gff" % (label1, label2)
        # aB.saveas(aB_name)

        AB_name = "%s.%s.gff" % (label1, label2)
        # AB.saveas(AB_name)


    #plt.show()

    # out_report = open(outfile_name + "_venn_diagram_report.txt" , "w")
    # out_report.write(results)
    # out_report.close()

    plt.savefig(outfile_name + ".pdf", format='pdf')
    return (Ab_name, AB_name, aB_name)



if __name__ == "__main__":
    main(sys.argv[1:])



















