import sys
import traceback
import math
import numpy
import scipy.stats


import pytransit.transit_tools as transit_tools
import pytransit.norm_tools as norm_tools
import pytransit.tnseq_tools as tnseq_tools


def error(x):
    print "Error: %s" % x

def usage():
    print """python %s -pt <Annotation .prot_table or .gff3>  -f <Comma seperated list of.wig files> [options]

        Optional Arguments:
            -n <string>     :=  Normalization method. Default: -n TTR
            --help          :=  Prints this help message and quits.
""" % sys.argv[0]




def main(args, kwargs):

    # Get arguments
    norm_method = kwargs.get("n", "TTR")

    # Data processing
    wiglist = kwargs["f"].split(",")
    (data, position) = tnseq_tools.get_data(wiglist)
    (norm_data, factors) = norm_tools.normalize_data(data, method=norm_method)
    G = tnseq_tools.Genes(wiglist, annotation, minread=1, ignoreCodon=False, nterm=0.0, cterm=0.0, data=wiglist, position=position)

    


if __name__ == "__main__":

     # Helper function to parse inputs
    (args, kwargs) = transit_tools.cleanargs(sys.argv[1:])

    # If asking for help, print help and exit
    if "-help" in kwargs:
        usage()
        sys.exit()

    # Check required arguments
    missingArgs = False
    if "f" not in kwargs:
        missingArgs = True
        error("Missing -f argument")
    if "pt" not in kwargs:
        missingArgs = True
        error("Missing -pt argument")

    if missingArgs:
        usage()
        sys.exit()

    main(args, kwargs)

    
