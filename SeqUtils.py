'''
Created on 16.02.2015

@author: David
'''

from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import Seq
from Bio.Data.CodonTable import generic_by_id

def stripCharsNotInList(seq, toKeep):
    charsToReplace = list()
    for c in seq:
        if not c in toKeep:
            charsToReplace.append(c)
    for c in charsToReplace:
        seq = seq.replace(c, '')        
    return seq

def getRemainderSuffix(seq):
    '''
    if sequence length is not a multiple of 3,
    return the last characters that are not part of any codon
    '''
    remainderLen = len(seq)%3
    if remainderLen == 0:
        return ""
    return seq[-remainderLen:]

def expandAmbiguous(seq):
    '''
    expand an ambigouous DNA sequence
    to a set of all possible unambiguous sequences
    '''
    res = list()
    res.append(seq)
    
    for i in range(len(seq)):
        resTmp = list()
        for s in res:
            for n in ambiguous_dna_values[s[i]]:
                sTmp = s[:i] + n + s[i+1:]
                resTmp.append(sTmp)
        res = resTmp
        
    return set(res)
        
#print(expandAmbiguous("NAG"))

def expandAmbiguousMult(seqs):
    '''
    expand mutiple (ambiguous) input sequences
    to all possible unambiguois sequences
    '''
    
    res = set()
    
    for s in seqs:
        rTemp = expandAmbiguous(s)
        res = res | rTemp
    
    return res

#print(expandAmbiguousMult({"NAA", "ANT"}))

def searchSubseq(seq, query):
    '''
    find all occurences of a query sequence,
    both on the coding and noncoding strand of a DNA sequence
    returns a set of (from, to) pairs identifying the found positions
    '''
    
    res = set()
    
    comp = str(Seq(seq).complement())
    revQuery = query[::-1]
    #print(revQuery)
    #print(comp)
    
    for i in range(len(seq)):
        if (seq[i:i+len(query)] == query) or (comp[i:i+len(query)] == revQuery):
            res.add((i, i+len(query)-1))
    
    return res

#print(searchSubseq("AAAGT", "AC"))

def searchSubseqs(seq, queries):
    '''
    search for multiple subsequences in DNA seq,
    both on coding and noncoding strand
    return set of positions (from, to) of all hits (regardless of query)
    '''    
    
    res = set()
    
    for q in queries:
        rTemp = searchSubseq(seq, q)
        res = res | rTemp
        
    return res

def getCodonsForRange(pos):
    f, t = pos
    res = range(int(f/3),int(t/3)+1)
    return set(res)

def getCodonsForRanges(ranges):
    res = set()
    for r in ranges:
        resT = getCodonsForRange(r)
        res = res | resT
        
    return res