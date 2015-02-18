#!/usr/local/bin/python3

import sys
import sequenceIO
import kazusaIO
import optimizing
import seqstats
import math

import re

# colored output
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

def main(argv):

    if (len(argv) != 4):
        print("Usage: optimize.py [sequence file] [source taxid] [target taxid]")
        sys.exit()

    seqfile = argv[1]
    src_taxid = argv[2]
    trg_taxid = argv[3]

    seq = sequenceIO.readFile(seqfile)
    seqr = sequenceIO.D2R(seq)
    print("== INPUT SEQUENCE, GC: ", round(seqstats.getGC(seq)*100), "%")
    # print(seq)

    cuSrc = kazusaIO.getCU(src_taxid)
    cuTrg = kazusaIO.getCU(trg_taxid)
    opt = optimizing.RelativeUsageCodonOptimizer(cuSrc, cuTrg)

    inCodons = re.findall('...', seqr) # list of codons
    outCodons = list()

    for c in inCodons:
        outCodons.append(opt.optimize(c))

    for i in range(0, len(outCodons)):
        if outCodons[i] == inCodons[i]:
            print(sequenceIO.R2D(inCodons[i]), end="")
        else:
            print(bcolors.OKGREEN, sequenceIO.R2D(inCodons[i]), bcolors.ENDC, end="", sep="")

    print("") #newline

    res = "".join(outCodons)
    resd = sequenceIO.R2D(res)
    print("== OPTIMIZED SEQUENCE, GC: ", round(seqstats.getGC(resd)*100), "%")

    for i in range(0, len(outCodons)):
        if outCodons[i] == inCodons[i]:
            print(sequenceIO.R2D(outCodons[i]), end="")
        else:
            print(bcolors.WARNING, sequenceIO.R2D(outCodons[i]), bcolors.ENDC, end="", sep="")

    print("") #newline

    print("== STATS: CHANGED CODONS: ", seqstats.getChanged(seq, resd)[0], "/", seqstats.getChanged(seq, resd)[1])

if __name__ == '__main__':
    main(sys.argv)
