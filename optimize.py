#!/usr/local/bin/python3

import sys
import sequenceIO
import kazusaIO
import optimizing

import re



def main(argv):
    
    if (len(argv) != 4):
        print("Usage: optimize.py [sequence file] [source taxid] [target taxid]")
    
    seqfile = argv[1]
    src_taxid = argv[2]
    trg_taxid = argv[3]    
    
    seq = sequenceIO.readFile(seqfile)
    seqr = sequenceIO.D2R(seq)
    print("== INPUT SEQUENCE:")
    print(seqr)
    
    cuSrc = kazusaIO.getCU(src_taxid)
    cuTrg = kazusaIO.getCU(trg_taxid)
    opt = optimizing.RelativeUsageCodonOptimizer(cuSrc, cuTrg)

    inCodons = re.findall('...', seqr) # list of codons
    outCodons = list()
    
    for c in inCodons:
        outCodons.append(opt.optimize(c))
        
    print("== OPTIMIZED SEQUENCE:")
    print("".join(outCodons))

if __name__ == '__main__':
    main(sys.argv) 