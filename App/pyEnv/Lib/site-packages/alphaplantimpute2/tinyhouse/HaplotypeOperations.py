import numba
from numba import jit, int8, int32, boolean, jitclass, float32, int64
import numpy as np

def setup_individual(ind):
        fillInPhaseFromGenotypes(ind.haplotypes[0], ind.genotypes)
        fillInPhaseFromGenotypes(ind.haplotypes[1], ind.genotypes)
        fillInGenotypesFromPhase(ind.genotypes, ind.haplotypes[0], ind.haplotypes[1])

def align_individual(ind):
    # Note: The note below is no longer true.
    # Note: We never directly set genotypes so no need to go from genotypes -> phase

    # Fill in phase from genotypes.
    fillInPhaseFromGenotypes(ind.haplotypes[0], ind.genotypes)
    fillInPhaseFromGenotypes(ind.haplotypes[1], ind.genotypes)

    # Fill in genotypes from phase, and fill in haplotypes from each other.
    fillInGenotypesFromPhase(ind.genotypes, ind.haplotypes[0], ind.haplotypes[1])
    fillInCompPhase(ind.haplotypes[0], ind.genotypes, ind.haplotypes[1])
    fillInCompPhase(ind.haplotypes[1], ind.genotypes, ind.haplotypes[0])

def ind_fillInGenotypesFromPhase(ind):
    #Note: We never directly set genotypes so no need to go from genotypes -> phase
    fillInGenotypesFromPhase(ind.genotypes, ind.haplotypes[0], ind.haplotypes[1])

def fillFromParents(ind):
    if ind.sire is not None:
        if ind.sire.genotypes is not None:
            sirePhase = getPhaseFromGenotypes(ind.sire.genotypes)
            fillIfMissing(ind.haplotypes[0], sirePhase)
    
    if ind.dam is not None:
        if ind.dam.genotypes is not None:
            damPhase = getPhaseFromGenotypes(ind.dam.genotypes)
            fillIfMissing(ind.haplotypes[0], damPhase)

@jit(nopython=True)
def fillIfMissing(orig, new):
    for i in range(len(orig)):
        if orig[i] == 9:
            orig[i] = new[i]

@jit(nopython=True)
def fillInGenotypesFromPhase(geno, phase1, phase2):
    for i in range(len(geno)):
        if geno[i] == 9:
            if phase1[i] != 9 and phase2[i] != 9:
                geno[i] = phase1[i] + phase2[i]
@jit(nopython=True)
def fillInCompPhase(target, geno, compPhase):
    for i in range(len(geno)):
        if target[i] == 9:
            if geno[i] != 9:
                if compPhase[i] != 9:
                    target[i] = geno[i] - compPhase[i]
    
@jit(nopython=True)
def fillInPhaseFromGenotypes(phase, geno):
    for i in range(len(geno)):
        if phase[i] == 9 :
            if geno[i] == 0: phase[i] = 0
            if geno[i] == 2: phase[i] = 1

@jit(nopython=True)
def getPhaseFromGenotypes(geno):
    phase = np.full(len(geno), 9, dtype = np.int8)
    for i in range(len(geno)):
        if phase[i] == 9 :
            if geno[i] == 0: phase[i] = 0
            if geno[i] == 2: phase[i] = 1
    return phase


# def ind_randomlyPhaseRandomPoint(ind):
#     maxLength = len(ind.genotypes)
#     midpoint = np.random.normal(maxLength/2, maxLength/10)
#     while midpoint < 0 or midpoint > maxLength:
#         midpoint = np.random.normal(maxLength/2, maxLength/10)
#     midpoint = int(midpoint)
#     randomlyPhaseMidpoint(ind.genotypes, ind.haplotypes, midpoint)

# def ind_randomlyPhaseMidpoint(ind, midpoint= None):
#     if midpoint is None: midpoint = int(len(ind.genotypes)/2)
#     randomlyPhaseMidpoint(ind.genotypes, ind.haplotypes, midpoint)

# @jit(nopython=True)
# def randomlyPhaseMidpoint(geno, phase, midpoint):

#     index = 0
#     e = 1
#     changed = False
#     while not changed :
#         if geno[midpoint + index * e] == 1:
#             phase[0][midpoint + index * e] = 0
#             phase[1][midpoint + index * e] = 1
#             changed = True
#         e = -e
#         if e == -1: index += 1
#         if index >= midpoint: changed = True

