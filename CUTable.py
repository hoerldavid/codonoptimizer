'''
Created on 08.02.2015

@author: david
'''

from Bio.Data.CodonTable import standard_dna_table
from Bio.Seq import Seq

class CUTable(object):
    '''
    classdocs
    '''


    def __init__(self, SPSUM_LABEL, spsum, code=standard_dna_table):
        '''
        Constructor
        SPSUM_LABEL is a string of whitespace separated codons
        spsusm is a string of whitespace separated numbers(frequencies)
        as present in the .spsum files of Kazusa
        Code is a Biopython codon-table object (TODO?: only standard at the moment)
        '''
        
        # codon >> aa map
        self.codon2AA = dict()
        # aa >> (codon, freq) list
        self.aa2CodonAndFreq = dict()
        
        labels = SPSUM_LABEL.split(" ")
        freqs = spsum.split(" ")
        
        for i in range(len(labels)):
            ## RNA labels to DNA labels
            lt = Seq(labels[i]).back_transcribe()
            labels[i] = str(lt)
            
            ## convert frequencies to ints
            freqs[i] = int(freqs[i])
            
        #print(code.forward_table)
        
        for i in range(len(labels)):
            if labels[i] in code.stop_codons:
                self.codon2AA[labels[i]] = 'STOP'
            else:
                self.codon2AA[labels[i]] = code.forward_table[labels[i]]
              
            tAA = self.codon2AA[labels[i]]
            if not self.aa2CodonAndFreq.get(tAA):
                self.aa2CodonAndFreq[tAA] = list()
            
            self.aa2CodonAndFreq[tAA].append((labels[i], freqs[i]))
        
        #print(self.codon2AA)
        #print(self.aa2CodonAndFreq)
        
    def print(self):
        '''
        print a (relatively) readable vesion of the CUTable to stdout
        '''
        for k, v in self.aa2CodonAndFreq.items():
            print(k + ":")
            for cd, fr in self.getCodonsForAA(k):
                print("\t" + cd + ": " + str(fr))
    
    def getAAForCodon(self, codon):
        '''
        get the AA corresponding to codon
        '''
        return self.codon2AA[codon]
    
    def getCodonRelativeUsage(self, codon):
        aa = self.codon2AA[codon]
        for cd, fr in self.getCodonsForAARelative(aa):
            if cd == codon:
                return fr
        
    def getCodonUsage(self, codon):
        aa = self.codon2AA[codon]
        for cd, fr in self.getCodonsForAA(aa):
            if cd == codon:
                return fr
            
    def getCodonsForAA(self, aa):
        '''
        get a list of (codon, rel. usage) for all codons coding for aa
        '''
        res = list()
        sumFr = 0
        for cd, fr in self.aa2CodonAndFreq[aa]:
            sumFr += fr
            
        for cd, fr in self.aa2CodonAndFreq[aa]:
            res.append((cd, fr/sumFr))
        
        return res
    
    def getCodonsForAARelative(self, aa):
        '''
        get a list of (codon, rel. usage %of max) for all codons coding for aa
        '''
        res = list()
        maxFr = 0
        for cd, fr in self.aa2CodonAndFreq[aa]:
            if fr > maxFr:
                maxFr = fr
            
        for cd, fr in self.aa2CodonAndFreq[aa]:
            res.append((cd, fr/maxFr))
        
        return res