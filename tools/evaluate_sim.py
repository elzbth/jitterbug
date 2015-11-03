#!/software/so/el6.3/PythonPackages-2.7.6/bin/python

import sys
import getopt
import pybedtools

from plot_gff_annots_2files import *



def main(argv):

    
    # get options passed at command line
    try:
        opts, args = getopt.getopt(argv, "g:i:o:s:h:")
    except getopt.GetoptError:
        print "evaluate_sim.py opts wrong"
        sys.exit(2)

    helpString = "evaluate_sim.py -g predicted_insertions.gff -i real_insertions.gff -o output_prefix -s [True|False] (save intersection files)"
    
    out_prefix = ""
    save_files = True
    gffFile = ""

    for opt, arg in opts:
        if opt == '-g':
            gffFile = arg
        elif opt == '-i':
            real_insertions = arg
        elif opt == '-o':
            out_prefix = arg
        elif opt == '-s':
            save_files = eval(arg)
        elif opt == '-h':
            print helpString

    if not gffFile:
        print helpString
        sys.exit(2)


    predicted = pybedtools.BedTool(gffFile)
    real = pybedtools.BedTool(real_insertions)


    # TP = predicted.intersect(real, wa=True, u=True)


    TP = predicted.intersect(real, u=True, wa=True)
    FP = predicted.intersect(real, v=True)
    FN = real.intersect(predicted, v=True)
    detected = real.intersect(predicted, u=True, wa=True)
    


    # TP = predicted + real
    # FP = predicted - real
    # FN = real - TP
    # detected = real + predicted

    if save_files:
        TP.saveas( out_prefix + "TP_" + gffFile )
        FP.saveas( out_prefix + "FP_" + gffFile )
        FN.saveas( out_prefix + "FN_" + gffFile)
        detected.saveas(out_prefix + "detected_" + gffFile)

    TP_count = float(len(TP))
    FP_count = float(len(FP))
    FN_count = float(len(FN))

    real_count = TP_count + FN_count

    

    PPV =  TP_count / (TP_count + FP_count)
    sensitivity = TP_count / (real_count)

    # print len(detected)
    print "real: %d\t%s" % (len(real), real_insertions)
    print "predicted: %d\t%s\n" % (len(predicted), gffFile)

    # print gffFile
    print "num_sims\tTP\tFP\tFN\tPPV\tsens" 
    print "%d\t%d\t%d\t%d\t%.2f\t%.2f\t"% (real_count, TP_count, FP_count, FN_count, PPV*100, sensitivity*100) 

    if save_files:
        plot(out_prefix + "TP_" + gffFile, "TP", out_prefix + "FP_" + gffFile, "FP", False, None)



if __name__ == "__main__":
    main(sys.argv[1:])
