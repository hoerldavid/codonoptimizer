##### FILE IO #####

def readFile(filename):
    '''
    Parse sequence from file to string.
    File may be in FASTA format or raw sequence
    '''
    fd = open(filename)
    c = fd.read()
    fd.close()
    
    # if file starts with '>', try to parse FASTA format
    if (c.startswith(">")):
        ret = parseFASTA(c)
    # else parse raw sequence
    else:
        ret = parseRaw(c)
    return ret


def parseFASTA(content):
    '''
    return first sequence in a FASTA-formated file
    '''    
    seqs = content.split(">")
    lines = seqs[1].split("\n") # consider only first sequence
    del lines[0] # ignore identifier line
    return("".join(lines).upper())
    
    
def parseRaw(content):
    '''
    return as uppercase with whitespace and newlines removed
    '''
    return content.strip("\n ").upper()
    
    
def D2R(seq):
    '''
    convert (sense) DNA to (m)RNA sequence
    '''
    seq2 = list(seq)
    for i in range(0, len(seq2)):
        # replace Ts with Us
        if (seq2[i] == 'T'):
            seq2[i] = 'U'
    return "".join(seq2)


def R2D(seq):
    '''
    convert (m)RNA to (sense) DNA sequence
    '''
    seq2 = list(seq)
    for i in range(0, len(seq2)):
        # replace Us with Ts
        if (seq2[i] == 'U'):
            seq2[i] = 'T'
    return "".join(seq2)