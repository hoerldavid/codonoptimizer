'''
Created on 07.02.2015

@author: david
'''

from os import listdir
from os.path import join

from CUTable import CUTable

class KazusaSPSUMHandler(object):
    '''
    classdocs
    '''


    def __init__(self, path=""):
        '''
        Constructor
        '''
        
        self.taxidToCodonUsage = dict()
        self.descToTaxidAndNCDS = dict()

        if not path:
            pass
        else:
            self.loadSPSUMFromPath(path)
        
        
    def loadSPSUMFromPath(self, path):
        '''
        
        '''
        
        files = listdir(path)
        
        for f in files:
            if f.endswith(".spsum"):
                self.readSPSUMFile(join(path, f))
            elif f == "SPSUM_LABEL":
                self.readSPSUM_LABEL(join(path, f))
        
        
    def readSPSUMFile(self, file):
        fd = open(file)
        
        while True:
            
            descLine = fd.readline()
            if not descLine: break
            
            cu = fd.readline()
            
            # split into taxid:desc:nrCDS
            taxid = descLine.split(":", 1)[0].strip()
            desc = descLine.split(":", 1)[1].strip()
            ncds = desc.rsplit(":", 1)[1].strip()
            desc = desc.rsplit(":", 1)[0].strip()
            
            self.taxidToCodonUsage[taxid] = cu.strip()
            self.descToTaxidAndNCDS[desc] = (taxid, ncds)
            
        fd.close()
        
    def readSPSUM_LABEL(self, file):
        
        fd = open(file)
        fd.readline() # omit first line
        
        self.SPSUM_LABEL = fd.readline().strip()
        
        fd.close()
    
        
    def search(self, query):
        '''
        '''
        
        res = dict()
        
        for desc, taxidNCS in self.descToTaxidAndNCDS.items():
            if desc.find(query) != -1:
                res[desc] = taxidNCS
                
        return res
        
    def getCUTable(self, taxid):
        '''
        '''
        res = CUTable(self.SPSUM_LABEL, self.taxidToCodonUsage[taxid])
        return res
        
        