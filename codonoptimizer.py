#!/usr/local/bin/python3

import sys

def D2Rcomp(codon):
    res = ""
    for c in codon:
        res += nucD2Rcomp(c)
    return res
    
def nucD2Rcomp(nuc):
    if (nuc == 'A'):
         return 'U'
    elif (nuc == 'T'):
        return 'A'
    elif (nuc == 'G'):
         return 'C'
    elif (nuc == 'C'):
        return 'G' 
    else:
        sys.exit("Unknown nucleotide! ERROR!") 
    
def R2Dcomp(codon):
    res = ""
    for c in codon:
        res += nucR2Dcomp(c)
    return res
    
def nucR2Dcomp(nuc):
    if (nuc == 'U'):
         return 'A'
    elif (nuc == 'A'):
        return 'T'
    elif (nuc == 'C'):
         return 'G'
    elif (nuc == 'G'):
        return 'C' 
    else:
        sys.exit("Unknown nucleotide! ERROR!")
        
def D2R(seq):
    seq2 = list(seq)
    for i in range(0, len(seq2)):
        if (seq2[i] == 'T'):
            seq2[i] = 'U'
    return "".join(seq2)

def R2D(seq):
    seq2 = list(seq)
    for i in range(0, len(seq2)):
        if (seq2[i] == 'U'):
            seq2[i] = 'T'
    return "".join(seq2)

def getAA(codon, table):
    return table[codon]

def getCodonRank(codon, table, aa_table):
    aa_name = getAA(codon, aa_table)
    
    for i in range(0, len(table[aa_name])):
        if (table[aa_name][i] == codon):
            return(i)
            
    return -1
    
def getCodonWithRank(aa, rank, table):
    return table[aa][rank]

def main(argv):

    sequence_path = argv[1]
    source_codon_path = argv[2]
    target_codon_path = argv[3]
    
    aa_table = dict()
    src_table = dict()
    trg_table = dict()
    
    src_file = open(source_codon_path)
    
    for l in src_file.readlines():
        ls = l.split(" ")
        aa_table[ls[0]] = ls[1]
        if (ls[1] in src_table.keys()):
            src_table[ls[1]].append((ls[0],ls[2]))
        else:
            src_table[ls[1]] = list()
            src_table[ls[1]].append((ls[0],ls[2]))
            
    src_file.close()
    
    src_table_2 = dict()
    
    for (k, v) in src_table.items():
        src_table_2[k] = list()
        for (k2, v2) in sorted(v, key=lambda x: (x[1])):
            src_table_2[k].append(k2)
            
        
    trg_file = open(target_codon_path)
    
    for l in trg_file.readlines():
        ls = l.split(" ")
        if (ls[1] in trg_table.keys()):
            trg_table[ls[1]].append((ls[0], ls[2]))
        else:
            trg_table[ls[1]] = list()
            trg_table[ls[1]].append((ls[0], ls[2]))
            
    trg_file.close()
    
    trg_table_2 = dict()
    
    for (k, v) in trg_table.items():
        trg_table_2[k] = list()
        for (k2, v2) in sorted(v, key=lambda x: (x[1])):
            trg_table_2[k].append(k2)
    
    
    seqf = open(sequence_path)
    
    sourceseq = ""
    resseq = ""
    
    print(src_table_2)
    print("\n\n")
    print(trg_table_2)
    
    while True:
        codon = seqf.read(3)
        
        if not codon:
          # print "End of file"
          break       
        
        sourceseq += codon.upper()
        codonR = D2R(codon.upper())
        
        AA = getAA(codonR, aa_table)
        rnk = getCodonRank(codonR, src_table_2, aa_table)
        resCodon = getCodonWithRank(AA, rnk, trg_table_2)
        
        resseq += R2D(resCodon)     
        
          
    print("Source sequence:    ", sourceseq)
    print("Optimized sequence: ", resseq)    
    

if __name__ == '__main__':
    main(sys.argv) 

