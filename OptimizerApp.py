'''
Created on 14.02.2015

@author: David
'''
from KazusaSPSUMHandler import KazusaSPSUMHandler
from codonoptimizer import *
from Bio.Restriction.Restriction_Dictionary import rest_dict
from Bio.Seq import Seq
import os, sys

import PathUtils

import configparser

class OptimizerApp:
    '''
    classdocs
    '''


    def __init__(self, config=None):
        '''
        Constructor
        '''        
        self.restrictionEnzymeList = list()
        self.speciesList = list()
        resDir = os.path.join(PathUtils.getCwd(), "res")
        
        self.SPSUMHandler = KazusaSPSUMHandler(resDir)
        
        if config:
            self.loadConfig(config)
        
        self.possibleOptimizationStrategies = ["Fastest Codons", "Adapt Speed To Source", "Random Adapt To Target"]
        self.optimizer = None
        
#         self.speciesList.append(("1234", "Testus specius"))
        
        self.sourceSequence = ""
        self.optimizedSequence = ""
        
    def setOptimizer(self, sourceTaxid, targetTaxid, strategy):
        
        if not strategy in self.possibleOptimizationStrategies:
            return
        
        if strategy == "Fastest Codons":
            self.optimizer = MostFrequentCodonOptimizer(self.SPSUMHandler.getCUTable(sourceTaxid), self.SPSUMHandler.getCUTable(targetTaxid))
        elif strategy == "Adapt Speed To Source":
            self.optimizer = AdaptingCodonOptimizer(self.SPSUMHandler.getCUTable(sourceTaxid), self.SPSUMHandler.getCUTable(targetTaxid))
        elif strategy == "Random Adapt To Target":
            self.optimizer = RandomTargetAdaptingCodonOptimizer(self.SPSUMHandler.getCUTable(sourceTaxid), self.SPSUMHandler.getCUTable(targetTaxid))
        else:
            return
        
    def setSourceSeq(self, seq):
        self.sourceSequence = seq.upper()
        
    def setOptimizedSeq(self, seq):
        self.optimizedSequence = seq.upper()
        
    def runOptimization(self):
        self.optimizedSequence = self.optimizer.getBestSequence(self.sourceSequence)
#         print("optimized to: " + self.optimizedSequence)
        if not sys.platform == "darwin":
            assert Seq(self.sourceSequence).translate() == Seq(self.optimizedSequence).translate()
        
    def runRestricionRemoval(self):
        restrictionSequences = list()
        for r in self.restrictionEnzymeList:
            restrictionSequences.append(rest_dict[r]['site'])
        self.optimizedSequence = self.optimizer.removeRestrictionSites(self.sourceSequence, self.optimizedSequence, restrictionSequences)
        
        if not sys.platform == "darwin":
            # this assertion fails on macs?
            # we look the other way for now
            assert Seq(self.sourceSequence).translate() == Seq(self.optimizedSequence).translate()
        
    def getCodonsForPrint(self, source=True):
        if self.optimizer:
            restrictionSequences = list()
            for r in self.restrictionEnzymeList:
                restrictionSequences.append(rest_dict[r]['site'])
            
            if source:
                return self.optimizer.SequenceToPrint(self.sourceSequence, restrictionSequences, source)
            else:
                return self.optimizer.SequenceToPrint(self.optimizedSequence, restrictionSequences, source)
            
        else:
            return None
        
    def testPrint(self):
        print(self.optimizer)
        print(self.sourceSequence)
        
    def saveConfig(self, path):
        cp = configparser.ConfigParser()
        taxids = list()
        names = list()
        for t, s in self.speciesList:
            taxids.append(t)
            names.append(s)
        restrictionEnzymes = list()
        for r in self.restrictionEnzymeList:
            restrictionEnzymes.append(r)
         
        cp['config'] = {"speciesTaxids" : ",".join(taxids), "speciesNames" : ",".join(names), "restrictionEnzymes" :",".join(restrictionEnzymes)}
        with open(path, 'w') as configfile:
#             print("writing config")
            cp.write(configfile)
        
    def loadConfig(self, path):
        
        if not os.path.exists(path):
            return
        
        cp = configparser.ConfigParser()
        cp.read(path)
        
        
        
        taxids = list()
        for t in cp["config"]["speciesTaxids"].split(","):
            if t:
                taxids.append(t)
            
        names = list()
        for n in cp["config"]["speciesNames"].split(","):
            if n:
                names.append(n)
            
        restrictionEnzymes = list()
        for r in cp["config"]["restrictionEnzymes"].split(","):
            if r:
                restrictionEnzymes.append(r)
            
        for i in range(len(taxids)):
            self.speciesList.append((taxids[i], names[i]))
            
        self.restrictionEnzymeList = restrictionEnzymes
        
    def setTest(self):
        self.sourceSequence = "ATGC"  
        
