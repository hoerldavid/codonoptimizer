import re

def getGC(seq):
    gccount = 0;
    for i in range(0, len(seq)):
        if (seq[i] == 'G' or seq[i] == 'C'):
            gccount += 1
            
    return float(gccount)/len(seq)
    
def getChanged(seq1, seq2):
    seq1c = re.findall('...', seq1)
    seq2c = re.findall('...', seq2)    
    changecount = 0;
    for i in range(0, len(seq1c)):
        if (seq1c[i] != seq2c[i]):
            changecount+=1
            
    return (changecount, len(seq1c))
        
    
    