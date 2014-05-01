#! /usr/bin/python2.7

import numpy as np
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

        #We first deal with the affected blocks
        for i,I in enumerate(argList):
            blocks_new.append([]) #Add an empty Block. It has index i
            for j in I:
                blocks_new[i].extend(self.blocks[j])
                unaffected.remove(j)

        blocks_new.extend([self.blocks[i] for i in unaffected])
        
        self.blocks = blocks_new
        self.sortBlocks()

class coalescent(object):
    '''Used to encode the jump-chain of a coalescent process. Supports updating when merger-events occur. Also supports adding mutations amd computing coalescent-statistics.'''
    
    def __init__(self,initialPartition,t_0 = 0.):
        self.jumps = [(t_0,initialPartition)]
        self.mutations = []
        self.t_max = t_0
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
        return x
    
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
        
    def addJump(self,t_new,newBlocks):
        if t_new < self.t_max:
            pass
            print 'New event at time t=',t_new,' will occur in the past! t_max = ',self.t_max,'. This is not supported!'
        else:
            self.jumps.append((t_new,newBlocks))
            self.k_current = newBlocks.countBlocks()
            self.t_max = t_new
    
    def mergerEvent(self,t,mergedBlocks):
        '''mergedBlocks is a list of lists; all elememts of the same sublist get merged, and a jump to the resulting partition at time "t" is added to the list of jumps'''
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
        'adds a Mutation. Mutations are encoded as tuples of the form (time,lineage); signifying that at time t, a mutation occurs to the n-th lineage (w.r.t. order of least elements starting at 0) of the coalescent'
        self.mutations.append(mutation)
    
    def getState(self,t):
        'return the value of the coalescent at time t'
        #TODO: Reread this method. It was written while drunk.
        i = 0
        while i < len(self.jumps) - 1 and t >= self.jumps[i+1][0]:
            i = i+1
        return deepcopy(self.jumps[i][1])

    def __getStateNoCoppy(self,t): #should not be called ny user
        i = 0
        while i < len(self.jumps) - 1 and t >= self.jumps[i+1][0]:
            i += 1
        return self.jumps[i][1]
    
    def computeSFS(self):
        '''Computes the Site Frquency Spectrum of the coalescent, and encodes it as an array. Indexing starts at 0 i.e. SFS[i-1] = \xi_i. Mutations happening after T_MRCA are ignored in this implementation'''
        SFS = np.zeros(self.jumps[0][1].n_blocks)
        for m in self.mutations:
            partition = self.__getStateNoCoppy(m[0])
            if 1 < partition.n_blocks and partition.n_blocks - 1 <= m[1]:
                SFS[ len(partition.blocks[m[1]]) - 1 ] += 1
        return SFS

class exchCoalWithMut(object):
    
    def __init__(self,n,tetha,T_max = infty):
        self.mutationRate = theta
        self.T_MRCA = inf
        self.SFS = np.zeros(n)
        self.coal = coalescent(partition([[i] for i in range(n)]))
        self.simulateCoalescent()
        
    def simulateCoalescent(self):
        t_current = 0
        while self.coal.k_current > 1 and t_current < T_max:
            #finish THis
            pass
            
    
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
