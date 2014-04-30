#! /usr/bin/python2.7

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
        Blocks are indexed (from 0) by order of least elements'''

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
    
    def mergeBlocksMultiple(self,*argList):
        '''arglist is a tuple (I_0,...,I_k) of integers. This method
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

def compareMinima(B1,B2):
    '''An auxilary function used for sorting blocks by order of least
    elements. Returns min(B1)-min(B2), i.e. a positive result when
    B1 > B2 ; where ">" denotes order of least elements.'''
    return min(B1)-min(B2)
