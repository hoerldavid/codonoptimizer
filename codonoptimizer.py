'''
Created on 16.02.2015

@author: David
'''
from CUTable import CUTable
import re
import SeqUtils
import heapq
import random

class CodonOptimizer(object):
    '''
    classdocs
    '''


    def __init__(self, sourceCU, targetCU):
        '''
        Constructor
        '''
        self.sourceCU = sourceCU
        self.targetCU = targetCU
        
    def score(self, sourceCodon, targetCodon):
        '''
        empty method to be overridden by implementations
        '''
        pass
    
    def getBestCodon(self, sourceCodon):
        
        aa = self.sourceCU.getAAForCodon(sourceCodon)
        targetCodonsAndScore = list()
        
        for codon, usage in self.targetCU.getCodonsForAA(aa):
            targetCodonsAndScore.append((codon, self.score(sourceCodon, codon)))
            
        return sorted(targetCodonsAndScore, key=lambda x:(x[1]))[0][0]
    
    def getBestSequence(self, sourceSequence):
        sourceCodons = re.findall('...', sourceSequence)
        remainder = SeqUtils.getRemainderSuffix(sourceSequence)
        
        result = ""
        for co in sourceCodons:
            result += self.getBestCodon(co)
            
        result += remainder
        return result
    
    def scoreSequence(self, sourceSequence, resultSequence):
        sourceCodons = re.findall('...', sourceSequence) # list of codons
        targetCodons = re.findall("...", resultSequence)
        
        scoreSum=0
        for i in range(len(sourceCodons)):
            scoreSum += self.score(sourceCodons[i], targetCodons[i])
            
        return scoreSum
    
    def getNextBestCodon(self, sourceCodon, optimizedCodon):
        aa = self.sourceCU.getAAForCodon(sourceCodon)
        targetCodonsAndScore = list()
        
        for codon, usage in self.targetCU.getCodonsForAA(aa):
            targetCodonsAndScore.append((codon, self.score(sourceCodon, codon)))
            
        targetCodonsAndScore = sorted(targetCodonsAndScore, key=lambda x:(x[1]))
        
        for i in range(len(targetCodonsAndScore)-1):
            if targetCodonsAndScore[i][0] == optimizedCodon:
                return targetCodonsAndScore[i+1][0]
        # no "worse" codon could be found
        return None
    
    def getPossibleOneStepChanges(self, sourceSeq, optimizedSeq, codonsToConsider):
        '''
        get all possible sequences resulting from substituting the codons to consider
        in the optimized Seq with the next best codons
        '''
        sourceCodons = re.findall('...', sourceSeq) # list of codons
        optCodons = re.findall("...", optimizedSeq)
        remainder = SeqUtils.getRemainderSuffix(optimizedSeq)
        
        res = list()
        
        for i in range(len(optCodons)):
            if i in codonsToConsider:
                nextBestCodon = self.getNextBestCodon(sourceCodons[i], optCodons[i])
                if nextBestCodon:
                    tCodons = optCodons
                    tCodons[i] = nextBestCodon
                    tSeq = "".join(tCodons)
                    tSeq += remainder
                    res.append(tSeq)
        
        return res
    
    def removeRestrictionSites(self, sourceSeq, optimizedSeq, restrictionSites):
        '''
        get the best sequence that does not contain any restriction sites
        by substituting codons in the optimizedSeq in a shortest-paths manner
        '''
        
        checkedSeqs = set()
        worklist = list()
        
        restrictionSites = SeqUtils.expandAmbiguousMult(restrictionSites)
        
        heapq.heappush(worklist, (self.scoreSequence(sourceSeq, optimizedSeq), optimizedSeq))
        
        while len(worklist) > 0:
            
            tScore, tSeq = heapq.heappop(worklist)
            restrictionLocations = SeqUtils.searchSubseqs(tSeq, restrictionSites)
            restrictionCodons = SeqUtils.getCodonsForRanges(restrictionLocations)
            
            if not restrictionCodons:
                
                return tSeq
            else:
                checkedSeqs.add(tSeq)
                possibleChanges = self.getPossibleOneStepChanges(sourceSeq, tSeq, restrictionCodons)
                if possibleChanges:
                    for tNewSeq in possibleChanges:
                        if not tNewSeq in checkedSeqs:
                            heapq.heappush(worklist, (self.scoreSequence(sourceSeq, tNewSeq), tNewSeq))
        
        return None
    
    def SequenceToPrint(self, seq, restrictionSites, source=True):
        '''
        return tuples of the form (codon, usage, isRestrictionSite) for the sequence
        '''
        codons = re.findall('...', seq)
        remainder = SeqUtils.getRemainderSuffix(seq)
        
        result = list()
        restrictionSites = SeqUtils.expandAmbiguousMult(restrictionSites)
        restrictionLocations = SeqUtils.searchSubseqs(seq, restrictionSites)
        restrictionCodons = SeqUtils.getCodonsForRanges(restrictionLocations)
        
        for i in range(len(codons)):
            usage = self.sourceCU.getCodonRelativeUsage(codons[i]) if source else self.targetCU.getCodonRelativeUsage(codons[i])
            result.append((codons[i], usage, i in restrictionCodons))
        
        if remainder:
            result.append((remainder, None, None))
        
        return result
    
    
    
class MostFrequentCodonOptimizer(CodonOptimizer):
    
    def score(self, sourceCodon, targetCodon):
        '''
        returns the negative relative usage (% of most frequent codon for aa) of a target codon
        negative usage means that the most frequent codon will get the lowest = best score
        '''
        aa = self.sourceCU.getAAForCodon(sourceCodon)
        codonsAndUsage = self.targetCU.getCodonsForAARelative(aa)
        
        for codon, usage in codonsAndUsage:
            if codon == targetCodon:
                return -usage
            
class AdaptingCodonOptimizer(CodonOptimizer):
    
    def score(self, sourceCodon, targetCodon):
        '''
        returns %-wise difference in relative usage
        returns rel.usage of source codon / relative usage of target if source is more frequent
        and vice versa
        --> two codons with the same relative frequency get score of 1, others a (worse) score > 1
        '''
        aa = self.sourceCU.getAAForCodon(sourceCodon)
        codonsAndUsageSource = self.sourceCU.getCodonsForAA(aa)
        codonsAndUsageTarget = self.targetCU.getCodonsForAA(aa)
        
        usageSource = 0
        usageTarget = 0
        
        for codon, usage in codonsAndUsageSource:
            if codon == sourceCodon:
                usageSource = usage
                
        for codon, usage in codonsAndUsageTarget:
            if codon == targetCodon:
                usageTarget = usage
                
        if usageSource == 0:
            return float("inf")
        elif usageSource < usageTarget:
            return usageTarget/usageSource
        else:
            return usageSource/usageTarget
            
class RandomTargetAdaptingCodonOptimizer(CodonOptimizer):
    
    def __init__(self, sourceCU, targetCU, threshold=0.0):
        self.threshold = threshold
        super(RandomTargetAdaptingCodonOptimizer, self).__init__(sourceCU, targetCU)
        
    def getCodonsForAAThresholded(self, aa):
        '''
        return an ordered list of (codon, usage) pairs, excluding those where usage is below threshold
        usages are normalized by the sum, to ensure that they add up to 1
        '''
        codonsAndUsageTarget = self.targetCU.getCodonsForAA(aa)
        codonsAndUsageTarget = sorted(codonsAndUsageTarget, key=lambda x:x[1])
        
        res = list()
        sumUs = 0.0
        for cd, us in codonsAndUsageTarget:
            if us > self.threshold:
                res.append((cd, us))
                sumUs += us
                
        for cd, us in res:
            us /= sumUs
            
        return res
        
    def getRandomOptimizedCodon(self, codon):
        aa = self.sourceCU.getAAForCodon(codon)
        possibleCodons = self.getCodonsForAAThresholded(aa)
        
        targetSum = random.random()
        sumSoFar = 0.0
        for co, us in possibleCodons:
            sumSoFar += us
            if targetSum <= sumSoFar:
                return co
                
        
    def getBestSequence(self, sourceSequence):
        sourceCodons = re.findall('...', sourceSequence)
        remainder = SeqUtils.getRemainderSuffix(sourceSequence)
        
        result = ""
        for co in sourceCodons:
            result += self.getRandomOptimizedCodon(co)
            
        result += remainder
        return result
        
    def removeRestrictionSites(self, sourceSeq, optimizedSeq, restrictionSites):
        '''
        get the best sequence that does not contain any restriction sites
        by re-randomizing until no restriction site remains
        '''
        
        codons = re.findall('...', optimizedSeq)
        remainder = SeqUtils.getRemainderSuffix(optimizedSeq)
        
        restrictionSites = SeqUtils.expandAmbiguousMult(restrictionSites)
        restrictionLocations = SeqUtils.searchSubseqs(optimizedSeq, restrictionSites)
        restrictionCodons = SeqUtils.getCodonsForRanges(restrictionLocations)
        
        # try to re-randomize a finite amount of times
        ITERMAX = 10000
        iteration = 0
        
        while iteration < ITERMAX:
            
            for i in restrictionCodons:
                codons[i] = self.getRandomOptimizedCodon(codons[i])
                
            tSeq = "".join(codons) + remainder
            
            tRL = SeqUtils.searchSubseqs(tSeq, restrictionSites)
            tRC = SeqUtils.getCodonsForRanges(tRL)
            
            if not tRC:
                return tSeq
                
            iteration += 1
        
        return None
    
        