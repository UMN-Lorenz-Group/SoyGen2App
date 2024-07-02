
import random
import numpy as np
import numba

try:
    from numba.experimental import jitclass
except ModuleNotFoundError:
    from numba import jitclass

from numba import njit, jit, int8, int32,int64, boolean, deferred_type, optional, float32
from collections import OrderedDict

try:
    profile
except:
    def profile(x): 
        return x

#####################################
###                               ###
###       Burrows Wheeler         ###
###                               ###
#####################################

class BurrowsWheelerLibrary():
    def __init__(self, haplotypeList):
        self.library = createBWLibrary(np.array(haplotypeList))
        self.hapList = haplotypeList
        self.nHaps = len(haplotypeList)
    def getHaplotypeMatches(self, haplotype, start, stop):
        nHaps, hapIndexes = getConsistentHaplotypes(self.library, haplotype, start, stop)
        haps = [(self.hapList[hapIndexes[index, 0]], hapIndexes[index, 1]) for index in range(nHaps)]
        return haps
    @profile
    def getBestHaplotype(self, weights, start, stop):
        index = getHaplotypesPlusWeights(self.library, weights, start, stop)

        return self.hapList[index][start:stop]

jit_BurrowsWheelerLibrary_spec = OrderedDict()
jit_BurrowsWheelerLibrary_spec['a'] = int64[:,:]
jit_BurrowsWheelerLibrary_spec['d'] = int64[:,:]
jit_BurrowsWheelerLibrary_spec['zeroOccPrev'] = int64[:,:]
jit_BurrowsWheelerLibrary_spec['nZeros'] = int64[:]
jit_BurrowsWheelerLibrary_spec['haps'] = int8[:,:]

@jitclass(jit_BurrowsWheelerLibrary_spec)
class jit_BurrowsWheelerLibrary():
    def __init__(self, a, d, nZeros, zeroOccPrev, haps):
        self.a = a
        self.d = d
        self.nZeros = nZeros
        self.zeroOccPrev = zeroOccPrev
        self.haps = haps

    def getValues(self):
        return (self.a, self.d, self.nZeros, self.zeroOccPrev, self.haps)


    def update_state(self, state, index):
        pass

    def get_null_state(self, value, index):
        if value == 0:
            lowerR = 0
            upperR = nZeros[stop-1]

        if value == 1:
            lowerR = nZeros[stop-1]
            upperR = nHaps

        pass

@njit
def createBWLibrary(haps):
    
    #Definitions.
    # haps : a list of haplotypes
    # a : an ordering of haps in lexographic order.
    # d : Number of loci of a[i,j+k] == a[i,-1, j+k]


    nHaps = haps.shape[0]
    nLoci = haps.shape[1]
    a = np.full(haps.shape, 0, dtype = np.int64)
    d = np.full(haps.shape, 0, dtype = np.int64)

    nZerosArray = np.full(nLoci, 0, dtype = np.int64)

    zeros = np.full(nHaps, 0, dtype = np.int64)
    ones = np.full(nHaps, 0, dtype = np.int64)
    dZeros = np.full(nHaps, 0, dtype = np.int64)
    dOnes = np.full(nHaps, 0, dtype = np.int64)
    
    nZeros = 0
    nOnes = 0
    for j in range(nHaps):
        if haps[j, nLoci-1] == 0:
            zeros[nZeros] = j
            if nZeros == 0:
                dZeros[nZeros] = 0
            else:
                dZeros[nZeros] = 1
            nZeros += 1
        else:
            ones[nOnes] = j    
            if nOnes == 0:
                dOnes[nOnes] = 0
            else:
                dOnes[nOnes] = 1
            nOnes += 1
    if nZeros > 0:
        a[0:nZeros, nLoci-1] = zeros[0:nZeros]
        d[0:nZeros, nLoci-1] = dZeros[0:nZeros]

    if nOnes > 0:
        a[nZeros:nHaps, nLoci-1] = ones[0:nOnes]
        d[nZeros:nHaps, nLoci-1] = dOnes[0:nOnes]

    nZerosArray[nLoci-1] = nZeros

    for i in range(nLoci-2, -1, -1) :
        zeros = np.full(nHaps, 0, dtype = np.int64)
        ones = np.full(nHaps, 0, dtype = np.int64)
        dZeros = np.full(nHaps, 0, dtype = np.int64)
        dOnes = np.full(nHaps, 0, dtype = np.int64)
    
        nZeros = 0
        nOnes = 0

        dZerosTmp = -1 #This is a hack.
        dOnesTmp = -1

        for j in range(nHaps) :

            dZerosTmp = min(dZerosTmp, d[j,i+1])
            dOnesTmp = min(dOnesTmp, d[j,i+1])
            if haps[a[j, i+1], i] == 0:
                zeros[nZeros] = a[j, i+1]
                dZeros[nZeros] = dZerosTmp + 1
                nZeros += 1
                dZerosTmp = nLoci
            else:
                ones[nOnes] = a[j, i+1]
                dOnes[nOnes] = dOnesTmp + 1
                nOnes += 1
                dOnesTmp = nLoci


        if nZeros > 0:
            a[0:nZeros, i] = zeros[0:nZeros]
            d[0:nZeros, i] = dZeros[0:nZeros]

        if nOnes > 0:
            a[nZeros:nHaps, i] = ones[0:nOnes]
            d[nZeros:nHaps, i] = dOnes[0:nOnes]
        nZerosArray[i] = nZeros


    #I'm going to be a wee bit sloppy in creating zeroOccPrev
    #Not defined at 0 so start at 1.
    zeroOccPrev = np.full(haps.shape, 0, dtype = np.int64)

    for i in range(1, nLoci):
        count = 0
        for j in range(0, nHaps):
            if haps[a[j, i], i-1] == 0:
                count += 1
            zeroOccPrev[j, i] = count


    library = jit_BurrowsWheelerLibrary(a, d, nZerosArray, zeroOccPrev, haps)
    return library

@jit(nopython=True, nogil=True)
def getConsistentHaplotypes(bwLibrary, hap, start, stop):
    a, d, nZeros, zeroOccPrev, haps = bwLibrary.getValues()
    nHaps = a.shape[0]
    nLoci = a.shape[1]


    intervals = np.full((nHaps, 2), 0, dtype = np.int64)
    intervals_new = np.full((nHaps, 2), 0, dtype = np.int64)
    
    nIntervals = 0
    nIntervals_new = 0
    
    #Haps go from 0 to nHaps-1. Loci go from start to stop-1 (inclusive).
    #The first hap with one is nZeros. The last hap with zero is nZeros -1.
    #Last loci is stop -1
    #These are split out because they represent *distinct* haplotypes.
    #Maybe could do this with tuple and list append but *shrug*.

    if hap[stop-1] == 0 or hap[stop-1] == 9:
        lowerR = 0
        upperR = nZeros[stop-1]
        if upperR >= lowerR:
            intervals[nIntervals, 0] = lowerR
            intervals[nIntervals, 1] = upperR
            nIntervals += 1
    
    if hap[stop-1] == 1 or hap[stop-1] == 9:
        lowerR = nZeros[stop-1]
        upperR = nHaps
        if upperR >= lowerR:
            intervals[nIntervals, 0] = lowerR
            intervals[nIntervals, 1] = upperR
            nIntervals += 1

    #Python indexing is annoying.
    #Exclude stop and stop-1, include start.
    #Intervals are over haplotypes.
    for i in range(stop-2, start-1, -1):
        # print(intervals[0:nIntervals,:])

        nIntervals_new = 0

        #Doing it on interval seems to make marginally more sense.
        for interval in range(nIntervals):

            int_start = intervals[interval, 0]
            int_end = intervals[interval, 1]

            if hap[i] == 0 or hap[i] == 9:
                if int_start == 0:
                    lowerR = 0
                else:
                    lowerR = zeroOccPrev[int_start-1, i+1] 
                upperR = zeroOccPrev[int_end-1, i+1] #Number of zeros in the region. 
                if upperR > lowerR: #Needs to be greater than. Regions no longer inclusive.
                    # print("Added 0:", int_start, int_end, "->>", lowerR, upperR)
                    intervals_new[nIntervals_new, 0] = lowerR
                    intervals_new[nIntervals_new, 1] = upperR
                    nIntervals_new += 1
            if hap[i] == 1 or hap[i] == 9:

                # of ones between 0 and k (k inclusive) is k+1 - number of zeros.

                if int_start == 0:
                    lowerR = nZeros[i]
                else:
                    lowerR = nZeros[i] + (int_start - zeroOccPrev[int_start-1, i+1]) 
                upperR = nZeros[i] + (int_end - zeroOccPrev[int_end-1, i+1])
                if upperR > lowerR:
                    # print("Added 1:", int_start, int_end, "->>", lowerR, upperR)
                    intervals_new[nIntervals_new, 0] = lowerR
                    intervals_new[nIntervals_new, 1] = upperR
                    nIntervals_new += 1
                # else:
                    # print(i, "interval rejected:", int_start, int_end, "->", upperR, lowerR)
        #This is basically intervals = intervals_new
        for j in range(nIntervals_new):
            intervals[j, 0] = intervals_new[j, 0]
            intervals[j, 1] = intervals_new[j, 1]
        nIntervals = nIntervals_new
        # print("Finished", i, "->", nIntervals)
    # print(intervals[0:nIntervals,:])

    hapIndexes = np.full((nHaps, 2), 0, dtype = np.int64)
    nHapsAssigned = 0
    for i in range(nIntervals):

        int_start = intervals[i, 0]
        int_end = intervals[i, 1]


        hapIndexes[nHapsAssigned, 0] = a[int_start,start]
        hapIndexes[nHapsAssigned, 1] = int_end - int_start
        nHapsAssigned +=1

    return (nHapsAssigned, hapIndexes)

# def printSortAt(loci, library):
#     a, d, nZeros, zeroOccPrev, haps = library.getValues()
#     vals = haps[a[:,loci],:]
#     for i in range(vals.shape[0]):
#         print(i, ' '.join(map(str, vals[i,:])) )
    
#     # print(zeroOccPrev[:,:])

# hapLib = [np.array([1, 0, 0, 0, 0, 0, 1], dtype = np.int8),
#           np.array([0, 1, 0, 0, 0, 1, 0], dtype = np.int8),
#           np.array([0, 1, 0, 0, 0, 1, 0], dtype = np.int8),
#           np.array([0, 1, 0, 0, 0, 1, 0], dtype = np.int8),
#           np.array([0, 0, 1, 0, 1, 0, 0], dtype = np.int8),
#           np.array([1, 1, 1, 0, 0, 0, 0], dtype = np.int8),
#           np.array([0, 0, 1, 0, 1, 0, 0], dtype = np.int8),
#           np.array([1, 1, 1, 0, 1, 0, 0], dtype = np.int8),
#           np.array([0, 0, 0, 1, 0, 0, 0], dtype = np.int8),
#           np.array([0, 1, 1, 1, 0, 0, 0], dtype = np.int8),
#           np.array([0, 1, 1, 1, 0, 0, 0], dtype = np.int8),
#           np.array([1, 1, 0, 1, 0, 0, 0], dtype = np.int8),
#           np.array([0, 0, 1, 0, 1, 0, 0], dtype = np.int8),
#           np.array([0, 0, 1, 0, 1, 0, 0], dtype = np.int8),
#           np.array([0, 0, 1, 0, 1, 0, 0], dtype = np.int8),
#           np.array([0, 0, 1, 0, 1, 0, 0], dtype = np.int8),
#           np.array([0, 0, 1, 0, 1, 0, 0], dtype = np.int8)]

# bwlib = BurrowsWheelerLibrary(hapLib)

# # printSortAt(0, bwlib.library)
# printSortAt(6, bwlib.library); print("")
# printSortAt(5, bwlib.library); print("")
# printSortAt(4, bwlib.library); print("")
# printSortAt(3, bwlib.library); print("")
# printSortAt(2, bwlib.library); print("")
# printSortAt(1, bwlib.library); print("")
# printSortAt(0, bwlib.library); print("")

# # print(bwlib.getHaplotypeMatches(haplotype = np.array([0, 0, 0], dtype = np.int8), start = 0, stop = 3))
# tmp = (bwlib.getHaplotypeMatches(haplotype = np.array([9, 9, 9, 9, 9, 9, 9], dtype = np.int8), start = 0, stop = 7))
# for key, value in tmp:
#     print(key, value)

@njit
def getConsistentHaplotypesBruteForce(bwLibrary, hap, start, stop):
    hap = hap[start:stop]

    a, d, nZeros, zeroOccPrev, haps = bwLibrary.getValues()
    recodedHaps = haps[:, start:stop]
    
    nHaps = recodedHaps.shape[0]
    nLoci = recodedHaps.shape[1]

    consistent = np.full(nHaps, 0, dtype = np.int64)

    #Last correct index
    lastIndex = -1
    for j in range(nHaps):

        #Otherwise, we have not enough information and need to search.
        add = True
        for i in range(nLoci):
            if hap[i] != 9 :
                if recodedHaps[j, i] != hap[i]:
                    add = False
        if add:
            consistent[j] = 1
    
    hapIndexes = np.full((nHaps, 2), 0, dtype = np.int64)
    nHapsAssigned = 0
    for i in range(nHaps):
        if consistent[i] > 0:
            # hapIndexes[nHapsAssigned, 0] = a[i,start]
            hapIndexes[nHapsAssigned, 0] = i
            hapIndexes[nHapsAssigned, 1] = consistent[i]
            nHapsAssigned +=1

    return (nHapsAssigned, hapIndexes)
@njit
def getHaplotypesPlusWeights(bwLibrary, weights, start, stop):
    #Weights are weights to the original haplotypes (haps)

    a, d, nZeros, zeroOccPrev, haps = bwLibrary.getValues()
    recodedWeights = weights[a[:, start]]

    nHaps = d.shape[0]
    nLoci = stop - start
    bestLoci = -1
    bestWeight = -1
    currentLoci = 0
    currentWeight = 0

    for j in range(nHaps):
        #Will need to double check this. This code will be slow!

        if d[j, start] < nLoci:
            #Process the last haplotype before moving on.
            if currentWeight > bestWeight :
                bestLoci = currentLoci
                bestWeight = currentWeight
            
            currentLoci = j
            currentWeight = 0

        currentWeight += recodedWeights[j]

    #Make sure we check the last haplotype.
    if currentWeight > bestWeight :
        bestLoci = currentLoci
        bestWeight = currentWeight


    #REMEMBER TO RECODE
    return a[bestLoci, start]

#Older version that doesn't use all of the meta data we have.
# @jit(nopython=True, nogil=True)
# def getConsistentHaplotypes(bwLibrary, hap, start, stop):
#     a, d, nZeros, zeroOccPrev, haps = bwLibrary.getValues()
#     hap = hap[start:stop]
#     recodedHaps = haps[a[:, start], start:stop]
    
#     nHaps = recodedHaps.shape[0]
#     nLoci = recodedHaps.shape[1]
#     consistent = np.full(nHaps, 0, dtype = np.int64)

#     lastCorrect = -1
#     firstError = nLoci + 1

#     #Last correct index
#     lastIndex = -1
#     for j in range(nHaps):
#         #Basically, we know how much overlap there was with the previous haplotype.
#         #We can use that to get a better estimate of where this one will be correct.

#         #By definition, all of 0 -> d[j, start]-1 inclusive is the same.
#         #All of 0 -> lastCorrect (inclusive) is correct.
#         #First error is the position of the first error. If firstError < nLoci, this is a problem. (nLoci is out of our bounds)
        
#         lastCorrect = min(lastCorrect, d[j, start]-1)
#         if firstError > d[j,start]-1: firstError = nLoci
        
#         #Two elif statements.
#         #First can we spot an error?
#         #Second if no error's exist, are we sure that the entire haplotype is right.
#         if firstError < nLoci: #or equal?
#             consistent[j] = 0
#             lastIndex = -1
#         #nLoci is the last position we care about (nLoci is out of our bounds). If nLoci is correct, then we're good).
#         elif lastCorrect >= nLoci-1:
#             #Adding in some short circuit code to prevent duplication
#             # lastIndex = -1 ###THIS LINE TOGGLES DUPLICATION PREVENTION
#             if lastIndex != -1:
#                 consistent[lastIndex] += 1
#             else:
#                 consistent[j] = 1
#                 lastIndex = j
#         else:
#             #Otherwise, we have not enough information and need to search.
#             #Last correct is the last known correct loci.
#             stopIters = False
#             i = lastCorrect
#             while not stopIters:
#                 i = i+1 #We know that the last correct loci is correct. So we're safe to start at last correct +1 
#                 if hap[i] != 9 and recodedHaps[j, i] != hap[i]:
#                     lastCorrect = i-1
#                     firstError = i
#                     stopIters = True
#                 elif i == nLoci-1:
#                     stopIters = True
#                     lastCorrect = nLoci-1
#                     firstError = nLoci 

#             if firstError < nLoci: 
#                 consistent[j] = 0
#                 lastIndex = -1

#             elif lastCorrect >= nLoci-1: #This will probably be nLoci-1 since that is where our search stops.
#                 consistent[j] = 1
#                 lastIndex = j
    
#     hapIndexes = np.full((nHaps, 2), 0, dtype = np.int64)
#     nHapsAssigned = 0
#     for i in range(nHaps):
#         if consistent[i] > 0:
#             hapIndexes[nHapsAssigned, 0] = a[i,start]
#             hapIndexes[nHapsAssigned, 1] = consistent[i]
#             nHapsAssigned +=1

#     return (nHapsAssigned, hapIndexes)
