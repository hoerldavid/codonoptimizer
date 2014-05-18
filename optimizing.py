

class CodonOptimizer:
    '''
    base class of codon optimizer.
    methods for loading tables are implemented here
    '''
    
    def __init__(self, cuSource, cuTarget):
        '''
        initialize with 2 codon usage tables
        '''
        self.aa_table = dict()
        self.src_table = dict()
        self.trg_table = dict()
    
        for l in cuSource.split("\n"):
            ls = l.split(" ") # split line into: codon, aa, relative usage, ...
            
            self.aa_table[ls[0]] = ls[1] # map codon->aa constructed from source table
            
            # create map aa->(codon, usage) for source organism
            if (ls[1] in self.src_table.keys()):
                self.src_table[ls[1]].append((ls[0],ls[2]))
            else:
                self.src_table[ls[1]] = list()
                self.src_table[ls[1]].append((ls[0],ls[2]))
                
                
        for l in cuTarget.split("\n"):
            ls = l.split(" ") # split line into: codon, aa, relative usage, ...
        
            # create map aa->(codon, usage) for source organism
            if (ls[1] in self.trg_table.keys()):
                self.trg_table[ls[1]].append((ls[0],ls[2]))
            else:
                self.trg_table[ls[1]] = list()
                self.trg_table[ls[1]].append((ls[0],ls[2]))
                
class MostUsedCodonOptimizer(CodonOptimizer):
    '''
    optimize sequence to use only the most used codons in target organism
    '''    
    def optimize(self, codon):
        aa = self.aa_table[codon]
        target_codons = self.trg_table[aa] #codons for same aa in target
        codons_sorted = sorted(target_codons, reverse=True, key=lambda x: (x[1])) 
        return(codons_sorted[0][0])
        
class RelativeUsageCodonOptimizer(CodonOptimizer):
    '''
    optimize sequence to minimize the difference (percent-wise) between relative usage in source and target organism.
    '''
    def optimize(self, codon):
        aa = self.aa_table[codon]
        
        for c in self.src_table[aa]:
            if c[0] == codon:
                srcUsage = float(c[1])
                
        # create list of possible codons (codon, rel.usage difference (max/min))
        targetCodons = list()
        for c in self.trg_table[aa]:
            if (srcUsage < float(c[1])):
                targetCodons.append((c[0], float(c[1])/srcUsage))
            else:
                targetCodons.append((c[0], srcUsage/float(c[1])))
                
        return(sorted(targetCodons, key=lambda x: (x[1]))[0][0])