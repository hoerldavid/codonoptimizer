##### FILE READING #####

def readFile(filename):
    file = open(filename)
    c = file.read()
    file.close()
    
    if (c.startswith(">")):
        ret = parseFASTA(c)
    else:
        ret = parseRAW(c)


def parseFASTA(content):    
    lines = content.split("\n")
    del lines[0]
    print("".join(lines))
    
def parseRaw(content):
    return 0