#!/usr/local/bin/python3

import sys
import sequenceIO
import kazusaIO

def main(argv):
    seq = sequenceIO.readFile(argv[1])
    #print(seq)
    print(kazusaIO.getCU("1423")) # b.subtilis

if __name__ == '__main__':
    main(sys.argv) 