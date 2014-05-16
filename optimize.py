#!/usr/local/bin/python3

import sys
import sequenceIO

def main(argv):
    sequenceIO.readFile(argv[1])


if __name__ == '__main__':
    main(sys.argv) 