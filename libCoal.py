#! /usr/bin/python2.7

import numpy as np
from scipy.special import binom
from copy import deepcopy

class partition(object):
    """Encodes a partition of the set [n]={0,...,n-1}, i.e. a
    decomposition into non-intersectiong sets whose union is [n]
    the argument initialBlocks should be a list of lists encoding the
    initial block-state e.g. [[1],[2],[3],[4,5]]"""
    def __init__(self, initialBlocks):
        self.blocks = list(initialBlocks)
        self.sortBlocks()
        self.isOrdered = True #set to false, if blocks are altered
        self.n_blocks = self.countBlocks()
        self.n_elements = self.countElements()
    
    def __str__(self):
        return "P:"+str(self.blocks)
    
    def __repr__(self):
        return "%s\n(%r)" % (self.__class__, self.__dict__)
    
    def numberPartition(self):
        lam = []
        for B in self.blocks:
            lam.append(len(B))
        lam.sort()
        return lam[::-1] #return lam in reversed order
    
    def sortBlocks(self):
        '''Sorts the blocks by order of least elements, and sorts the
        elements of each block.'''
        
        'Sort the list of blocks by order of least elements'
        self.blocks.sort(compareMinima) # Sort the blocks
        #self.blocks.sort(lambda B1,B2:min(B1)-min(B2)) # Old Version
        
        'sort each individual block'
        for B in self.blocks:
            B.sort()

        self.isOrdered = True

    def countBlocks(self):
        return len(self.blocks)

    def countElements(self):
        return sum(map(len,self.blocks))
    
    def hasDistinctElements(self):

        listOfElements = []

        for B in self.blocks:
            listOfElements.extend(B)

        for i in listOfElements:
            if listOfElements.count(i) > 1:
                return False

        return True
    
    def mergeBlocksSingle(self,*argList):
        '''Merge all blocks with index in argList (a tuple of integers).
        Blocks are indexed (from 0) by order of least elements.
        For example, mergeBlocksSingle(0,2) would merge the 1st and 3rd
        blocks of the partition in question.
        '''

        #ensure that blocks are sorted
        if not(self.isOrdered):
            self.sortBlocks()

        #ensure that we do not get indexing errors when looping
        argList = list(argList) #argList would otherwise be immutable
        argList.sort()
        argList.reverse()

#        blocks_new = self.blocks
        B_new = []

        for i in argList:
            #blocks_new.remove(self.blocks[i])
            #B_new.extend(self.blocks[i])
            B_new.extend(self.blocks.pop(i))

        self.blocks.append(B_new)
        self.n_blocks = self.countBlocks()
        self.sortBlocks()
    
    def mergeBlocksMultiple(self,argList):
        '''arglist is a tuple (I_0,...,I_k) of integer-lists. This method
        merges all blocks with index j in the same I_i into one large
        block. It is written under the assumption that no index is out
        of range, and no index occurs more than once
        
        example input:
        >>> pi = partition([[1],[2],[3],[4],[5]])
        >>> print "blocks = ",pi.blocks
        blocks =  [[1], [2], [3], [4], [5]]
        >>> pi.mergerMultiple([0,2,4],[1,3])
        >>> print "blocks = ",pi.blocks
        blocks =  [[1, 3, 5], [2, 4]]
        '''

        unaffected = range(len(self.blocks))
        blocks_new = []

        #ensure that blocks are sorted
        if not(self.isOrdered):
            self.sortBlocks()

        #remove empty blocks from ArgList
        argList = filter(lambda x: len(x) > 0 , argList)
        
        #We first deal with the affected blocks
        for i,I in enumerate(argList):
            blocks_new.append([]) #Add an empty Block. It has index i
            for j in I:
                blocks_new[i].extend(self.blocks[j])
                unaffected.remove(j)

        blocks_new.extend([self.blocks[i] for i in unaffected])
        
        self.blocks = blocks_new
        self.n_blocks = self.countBlocks()
        self.sortBlocks()

class coalescent(object):
    '''Used to encode the jump-chain of a coalescent process. Supports updating
    when merger-events occur. Also supports adding mutations amd computing
    coalescent-statistics.'''
    
    def __init__(self,initialPartition,t_0 = 0.):
        self.jumps = [(t_0,initialPartition)]
        self.mutations = []
        self.t_lastJump = t_0
        self.k_current = initialPartition.countBlocks()
        self.times = self.getJumptimes()
        self.partitionChain = self.getPartitionChain()
    
    def __str__(self):
        l = []
        for event in self.jumps:
            l.append((event[0],str(event[1])))
        return str(l)
    
    def getJumptimes(self):
        times = []
        for x in self.jumps:
            times.append(x[0])
        return times
    
    def getPartitionChain(self):
        chain = []
        for x in self.jumps:
            chain.append(deepcopy(x[1]))
        return chain
    
    def getNumberPartitionChain(self):
        chain = []
        for x in self.jumps:
            chain.append(x[1].numberPartition())
        return chain
        
    def getMutations(self):
        return deepcopy(self.mutations)
        
    def addJump(self,t_new,newBlocks):
#        if t_new < self.t_lastJump:
        if t_new < self.jumps[-1][0] :
            pass
            print 'New event at time t=',t_new,' will occur in the past! t_max = ',self.t_lastJump,'. This is not supported!'
        else:
            self.jumps.append((t_new,newBlocks))
            self.k_current = newBlocks.countBlocks()
            self.t_lastJump = t_new
    
    def mergerEvent(self,t,mergedBlocks):
        '''mergedBlocks is a list of lists; all elememts of the same sublist
        get merged, and a jump to the resulting partition at time "t" is added
        to the list of jumps'''
        #self.jumps[-1][1] is the partition after the last jump
        newPartition = deepcopy(self.jumps[-1][1])
        newPartition.mergeBlocksMultiple(mergedBlocks)
        self.addJump(t,newPartition)

    def coagulationEvent(self,t,part2):
        '''(similar to mergerEvent) adds (t,coag(part1,part2)) to the list of
        jumps, where "part1" signifies the current state of the coalescent'''
        newPartition = coag(self.jumps[-1][1],part2)
        self.addJump(t,newPartition)
    
    def addMutation(self,mutation):
        '''adds a Mutation. Mutations are encoded as tuples of the form
        (time,lineage); signifying that at time t, a mutation occurs to the
        n-th lineage (w.r.t. order of least elements starting at 0) of the
        coalescent'''
        self.mutations.append(mutation)
    
    def getState(self,t):
        'return the value of the coalescent at time t'
        i = 0
        while i < len(self.jumps) - 1 and t >= self.jumps[i+1][0]:
            i = i+1
        return deepcopy(self.jumps[i][1])

    def __getStateNoCoppy(self,t):
        '''Returns the last object in self.jumps with a time-index <= t
        It does NOT return a copy of said object, and is only meant facilitate
        that other methods that do not alter the internal state can be sped up
        by skipping the copy-step'''
        i = 0
        while i < len(self.jumps) - 1 and t >= self.jumps[i+1][0]:
            i += 1
        return self.jumps[i][1]
    
    def computeSFS(self):
        '''Computes the Site Frquency Spectrum of the coalescent, and encodes
        it as an array. Indexing starts at 0 i.e. SFS[i-1] = xi^(n)_i.
        Mutations happening after T_MRCA are ignored in this implementation'''
        SFS = np.zeros(self.jumps[0][1].n_elements)
        for m in self.mutations:
#            partition = self.__getStateNoCoppy(m[0])
            partition = self.__getStateNoCoppy(m[0])
            
#            if 1 < partition.n_blocks and partition.n_blocks - 1 <= m[1]:
            if 1 < partition.n_blocks and m[1] < partition.n_blocks:
#            if True:
                SFS[ len(partition.blocks[m[1]]) - 1 ] += 1
        return SFS
        
    def computeNormalizedSFS(self,factor=False):
        '''
        Returns the SFS/factor. If no factor is passed, #mutations+1 is used
        '''
        SFS = self.computeSFS()
#        SFS = self.SFS
        
        if not(factor):        
#            '''This normalization-factor was suggested by Bjarki Eldon, and
#            seems reasonable. It is not the same factor as the one used in the
#            paper by Birkner, Blath and Eldon (2013)'''
#            normalizationFactor = (len(self.mutations)+1)**-1
            normalizationFactor = (max(len(self.mutations), 1) )**-1
        else:
            normalizationFactor = factor**-1
#        for i,xi in enumerate(SFS):
#            normSFS[i] = xi*normalizationFactor
        
        return normalizationFactor * SFS
    
    def computeTreeLength(self):
        '''returns the total length of lineages'''
        l = 0.
        for i in range(1,len(self.jumps)):
            # RHS: number of lineages multiplied by length of time between jumps
            l += self.jumps[i-1][1].n_blocks * (self.jumps[i][0] - self.jumps[i-1][0])
        return l
    
    def getJumpsInTimeRange(self,t_min,t_max):
        return filter(lambda x: t_min <= x[0] and x[0] <= t_max, self.jumps)
#        return deepcopy(filter(lambda x: t_min <= x[0] and x[0] <= t_max, self.jumps))
    def getMutationsInTimeRange(self,t_min,t_max):
        return filter(lambda x: t_min <= x[0] and x[0] <= t_max, self.mutations)

class simExchCoalWithMut(object):
    '''
    A parrent class for all coalescent simulations, containing all comon methods and fields.
    '''
    
    IgnoreInvisibleJumps = True
    
    def __init__(self,n,mutationRate = 0,T_max = float('inf'),*args):
        '''
        n = number of individuals
        mutationRate = Mutation rate of the coalescent.
        T_Max = Time-horizon of the coalescent.
        '''
        self.n = n
        self.mutationRate = mutationRate
#        self.mergerRate = mergerRate
        self.T_max = T_max
        self.T_MRCA = float('inf')
        self.args = args
        self.SFS = np.zeros(n)
        self.coal = coalescent(partition([[i] for i in range(n)]))
#        self.simulateCoalescent(self.args)
        self.simulateCoalescent()
        
    def simulateCoalescent(self):
        t_current = 0

        keepGoing = True
        while self.coal.k_current > 1 and keepGoing:
            
            t_old = t_current
            
            #Generate Next jump-event. SampleJumps returns a tuple of the form (t,Indexes). It simplementation differs in each subclass.
            jumpEvent = self.sampleJumps()
            t_current = jumpEvent[0]
            self.coal.t_lastJump = t_current
            
            #If the jumpEvent is invisble in the sense that no blocks merge, we ignore it (but add the time elapsed).
            while self.IgnoreInvisibleJumps and not self.validateMergers(jumpEvent[1]):
                jumpEvent = self.sampleJumps()
                t_current = jumpEvent[0]
                self.coal.t_lastJump = t_current
            
            if t_current > self.T_max:
                keepGoing = False
            
            if keepGoing:
                #Simulate all mutations that have happend in the time between the two jumps and add them.
                delta_t = t_current - t_old
                no_new_mutations = np.random.poisson(self.mutationRate * self.coal.k_current * delta_t)
                affectedLineages = np.random.randint(0,self.coal.k_current,no_new_mutations)
                mutationTimes = np.random.uniform(t_old,t_current,no_new_mutations)
                for i in range(no_new_mutations):
                    self.coal.addMutation((mutationTimes[i],affectedLineages[i]))
                
                #Add the jump itself
                self.coal.mergerEvent(jumpEvent[0],jumpEvent[1],)

            if self.coal.k_current == 1:
                self.T_MRCA = self.coal.t_lastJump
        
        self.SFS = self.coal.computeSFS()
        
    def sampleJumps(self):
        "Rewrite depending on what kind of coalescent we want to simulate"
        pass
    
    def sampleJumps_Kingman(self,rate=1.):
        'Samples a jump of kingmans coalescent, whereby "rate" is the pairwise merger rate.'
        currentMergerRate= rate*binom(self.coal.k_current,2)
        waitingTime=np.random.exponential(currentMergerRate**-1)
        t = self.coal.t_lastJump + waitingTime
        affectedBlocks = list(np.random.choice(self.coal.k_current,2,replace=False))
        mergers = self.split(affectedBlocks)
        return (t,mergers)
    
    def sampleJumpsLambdaCTMC(self,P,q):
        '''Samples a jump of an arbitrary exchangeable Lambda-coalescent given
        The distribution of the block-counting-process.
            P[i,j] = Probability of jump from i to j
            q[i] = total rate of jumps out of i
        '''
        t = self.coal.t_lastJump
        b = self.coal.k_current
        waitingTime = np.random.exponential(q[b]**-1)
        t += waitingTime
        b_New = np.random.choice(b,p=P[b,:b])
        affectedBlocks = list(np.random.choice(b,size=b-b_New+1,replace=False))
        mergers = self.split(affectedBlocks)
        return (t,mergers)
    
    def sampleJumps_LambdaPointMeasure(self,phi,rate=1.):
        "Samples Jumps of a Lambda-coalescent, where Lambda = rate * dirac_phi"
        t = self.coal.t_lastJump
        t += np.random.exponential(rate**-1)
        #TODO: Figure out if this is the correct rate

## OLD AND WRONG
#        keepGoing = True
#        i=0
#        while keepGoing:
#            i += 1        
#        affectedBlocks = []
#        for i in range(self.coal.k_current):
#            if np.random.uniform() <= phi:
#                affectedBlocks.append(i)
#            keepGoing = False
#            if self.validateMergers(mergers):
#                keepGoing = False
#            if i>5000:
#                print "more than 5000 loops. WTF? phi=",phi
#                break
##
        r = np.random.uniform(size=self.coal.k_current)
        affectedBlocks = filter(lambda i: r[i]<=phi,range(self.coal.k_current))
        mergers = self.split(affectedBlocks)
        return (t, mergers)
    
    def validateMergers(self,mergers):
        #Tests if the list of lists "mergers" contains at least one list with more than one element.
        return len(filter(lambda x: len(x) > 1, mergers)) > 0
    
    def split(self,affectedBlocks):
        #Auxilary function for sampleJumps_LambdaPointMeasure. Rewrite if dealing with coalescents admitting simultanious mergers.
        return [affectedBlocks]
    
class simulateKingman(simExchCoalWithMut):
    '''args[0] = merger-rate of two fixed lineages'''
    def sampleJumps(self):
        if len(self.args) > 0 :
            rate=self.args[0]
        else:
            rate = 1
        return self.sampleJumps_Kingman(rate)

class simulateLambdaPoint(simExchCoalWithMut):
    '''Simulates a lambda coalescent, where lambda = \phi^2 dirac_phi,
    Where phi is an element of (0,1]. The argiuments are as follows:
    args[0] = phi
    The total rate of mergers in this case is '''
    def sampleJumps(self):
        phi = self.args[0]
        return self.sampleJumps_LambdaPointMeasure(phi)

class simulateLambdaPoint_FourWay(simulateLambdaPoint):

    def split(self,affectedBlocks):
        return fourwaySplit(affectedBlocks)
    
class simulateLambdaBeta(simulateLambdaPoint):
    '''Simulate a beta-coalescent, i.e. lambda = beta(2-alpha,alpha), where
    0 < alpha < 2.
        arg[0] = alpha
        arg[1] = P-matrix of jumps
        arg[2] = jump-rates'''
    def sampleJumps(self):
#        #OLD AND WRONG
#        phi = np.random.beta(2-self.args[0],self.args[0])
#        return self.sampleJumps_LambdaPointMeasure(phi)
        #NEW AND LESS WRONG
        return self.sampleJumpsLambdaCTMC(self.args[1],self.args[2])
#        t = self.coal.t_lastJump
#        b = self.coal.k_current
#        t += np.random.exponential(self.args[2][b]**-1)
#        b_New = np.random.choice(b,p=self.args[1][b,:b])
#        affectedBlocks = list(np.random.choice(b,size=b-b_New+ 1,replace=False))
#        return (t,[affectedBlocks])

class simulateLambdaBeta_FourWay(simulateLambdaBeta):
    
    def split(self,affectedBlocks):
        return fourwaySplit(affectedBlocks)

class simulateLambdaEldonWakely(simExchCoalWithMut):
    'Here Lambda = 2/(2+phi^2) dirac_0 + phi^2/(2+pgi^2) dirac_phi. args[0] = phi'
    
    def sampleJumps(self):
        phi = self.args[0]
        kingmanJump = self.sampleJumps_Kingman(rate=2./(2+(phi**2)))
        nonKingmanJump = self.sampleJumps_LambdaPointMeasure(phi,rate=1./(2+(phi**2)))
        if kingmanJump[0] < nonKingmanJump[0]:
            return kingmanJump
        else:
            return nonKingmanJump

class simulateLambdaEldonWakely_FourWay(simulateLambdaEldonWakely):
    
    def split(self,affectedBlocks):
        return fourwaySplit(affectedBlocks)

class simulateXiDirichlet(simExchCoalWithMut):
    '''
    args[0] =   beta (non-negative constant)
    args[1] =   gamma (a function that takes no arguments and outputs non-
                negative numbers (typically random) 
    '''
    def simulateCoalescent(self):
        kingmanSim = simulateKingman(self.n,self.mutationRate,float('inf'),1.)
        subordinator = compensatedPossionProcess(T_max=min(kingmanSim.T_MRCA,self.T_max),beta=self.args[0],gamma=self.args[1])
        # since the subordinator has drift 1 and positive jumps, S(T_mrca(kingman)) >= T_mrca(Subordinated Kingman)        

        subordinatorJumps = [(0.,0.)] + subordinator.jumps #we add a jump of magnitude 0 in order to simplify the algorithm
        for i,J in enumerate(subordinatorJumps):
            sumOfJumpsSoFar = sum([x[1] for x in subordinatorJumps[:i]])
            subordinatorPriorToJump = subordinator.drift*J[0] + sumOfJumpsSoFar
            subordinatorAfterJump = subordinatorPriorToJump + J[1]

            stateAfterJump = kingmanSim.coal.getState(subordinatorAfterJump)
            if stateAfterJump.n_blocks < self.coal.jumps[-1][1].n_blocks:
                self.coal.addJump(J[0],stateAfterJump)    

            if i==len(subordinatorJumps) - 1:
                nextJump = float('inf')
            else:
                nextJump = subordinatorAfterJump + subordinator.drift * (subordinatorJumps[i+1][0] - J[0])
            
            for event in kingmanSim.coal.getJumpsInTimeRange(subordinatorAfterJump,nextJump):
                eventTime = (event[0]-sumOfJumpsSoFar)/subordinator.drift
                # regarding the above:
                # if    T = drift * t + Sum_{k=1}^{N_t} J_k
                # then, t = (T - Sum_{k=1}^{N_t} J_k) / drift
                if event[1].n_blocks < self.coal.jumps[-1][1].n_blocks:
                    self.coal.addJump(eventTime,deepcopy(event[1]))

            for mutation in kingmanSim.coal.getMutationsInTimeRange(subordinatorAfterJump,nextJump):
                mutationTime = (mutation[0]-sumOfJumpsSoFar)/subordinator.drift
                self.coal.addMutation((mutationTime,mutation[1]))
        
        if self.coal.k_current == 1:
            self.T_MRCA = self.coal.t_lastJump
        
        self.SFS = self.coal.computeSFS()

class compensatedPossionProcess(object):
    '''Simulate a compensated Possion process with levy-measure beta * gamma
    and drift 1. Beta is a constant scaling the frequency of jump-events, and
    gamma is a measure'''
    
    def __init__(self,T_max,beta=1.,gamma=np.random.exponential,drift=1.):
        '''beta should be a non-negative number
        gamma should be a fungtion that takes no arguments and returns a
        positive number
        The process is siumlated on the time-interval [0,T_max]
        drift is assumed to be a constant term in this case'''
        self.beta = beta
        self.gamma = gamma
        self.T_max = T_max
        self.drift = drift
        self.jumps = self.simulateJumps()
        # jumps is a list of tuples of the form (t_jump,size_jump)
    
    def simulateJumps(self):
        jumps = []
        t = 0
        while True:
            t += np.random.exponential(self.beta**-1)
            if t < self.T_max:
                jumps.append((t,self.gamma()))
            else:
                break
        return jumps

    def getState(self,t):
        if t > self.T_max:
            t = self.T_max # the process is killed at time T_max
        S_t = self.drift*t
        S_t += sum([J[1] for J in filter(lambda x: x[0] <= t,self.jumps)])
        return S_t
    
    def countJumpsUntil(self,t):
        return len(filter(lambda J: J[0] <= t,self.jumps))
    
    def countJumpsInInterval(self,t_min,t_max):
        return len(filter(lambda J: t_min <= J[0] and J[0] <= t_max,self.jumps))

    def sumOfJumpsInInterval(self,t_min,t_max):
        return sum([J[1] for J in filter(lambda x: t_min <= x[0] and x[0] <= t_max,self.jumps)])
    
    def timeOfNextJump(self,t):
        if len(self.jumps) == 0: #case: No jumps have occured
            return float('inf')
        if t > self.jumps[-1][0]: # case: t > largest jumptime.
            return float('inf')
        i = 0
        while self.jumps[i][0] <= t:
            i += 1
        return self.jumps[i][0]
    
    def resample(self):
        self.jumps = self.simulateJumps()
        
def compareMinima(B1,B2):
    '''An auxilary function used for sorting blocks by order of least
    elements. Returns min(B1)-min(B2), i.e. a positive result when
    B1 > B2 ; where ">" denotes order of least elements.'''
    return min(B1)-min(B2)

def coag(part1,part2):
    '''returns the new partition resulting from coagulating part1 w.r.t.
    part2 as described in the lecture-notes by Bertoin
    URL: http://www.fim.math.ethz.ch/lectures/Lectures_Bertoin.pdf'''
    part_new = deepcopy(part1)
    part_new.mergeBlocksMultiple(part2.blocks)
    return part_new
    
def fourwaySplit(affectedBlocks):
    '''
    IN: a list of numbers.
        e.g. [1,2,3,4,5,6,7,8,9,10]
    OUT: a uniformly chosen random partition of the input into four blocks.
        e.g. [[1], [4, 10], [3, 5, 7, 9], [2, 6, 8]]
    '''
    l = [[],[],[],[]]
    for i in affectedBlocks:
        l[np.random.randint(4)].append(i)
    return l