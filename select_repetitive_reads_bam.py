#!/soft/bin/python
from BamReader import *
import getopt,sys

def main(argv):


    try:
        opts, args = getopt.getopt(argv, "i:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)



    # Get parameters from command line
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt == '-i':
            bam_file_name = arg



    bam_reader = BamReader(bam_file_name, None)
    bam_reader.output_repetitive_reads()


if __name__ == "__main__":
    main(sys.argv[1:])
