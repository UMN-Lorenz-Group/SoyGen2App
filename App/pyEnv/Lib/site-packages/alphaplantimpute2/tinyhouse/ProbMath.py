from numba import jit
import numpy as np
import collections

def getGenotypesFromMaf(maf) :
    nLoci = len(maf)
    mafGenotypes = np.full((4, nLoci), .25, dtype = np.float32)

    mafGenotypes[0,:] = (1-maf)**2
    mafGenotypes[1,:] = maf*(1-maf)
    mafGenotypes[2,:] = (1-maf)*maf
    mafGenotypes[3,:] = maf**2

    return mafGenotypes

def getGenotypeProbabilities_ind(ind, args = None, log = False):
    if args is None:
        error = 0.01
        seqError = 0.001
        sexChromFlag = False
    else:
        error = args.error
        seqError = args.seqerror
        sexChromFlag = getattr(args, "sexchrom", False) and ind.sex == 0 #This is the sex chromosome and the individual is male.

    if ind.reads is not None:
        nLoci = len(ind.reads[0])
    if ind.genotypes is not None:
        nLoci = len(ind.genotypes)
    if not log:
        return getGenotypeProbabilities(nLoci, ind.genotypes, ind.reads, error, seqError, sexChromFlag)
    else:
        return getGenotypeProbabilities_log(nLoci, ind.genotypes, ind.reads, error, seqError, sexChromFlag)


def getGenotypeProbabilities(nLoci, genotypes, reads, error = 0.01, seqError = 0.001, useSexChrom=False):
    vals = np.full((4, nLoci), .25, dtype = np.float32)
    if type(error) is float:
        error = np.full(nLoci, error)
    if type(seqError) is float:
        seqError = np.full(nLoci, seqError)

    errorMat = generateErrorMat(error)

    if genotypes is not None:
        setGenoProbsFromGenotypes(genotypes, errorMat, vals)
        
    if reads is not None:
        seqError = seqError
        log1 = np.log(1-seqError)
        log2 = np.log(.5)
        loge = np.log(seqError)
        valSeq = np.array([log1*reads[0] + loge*reads[1],
                                log2*reads[0] + log2*reads[1],
                                log2*reads[0] + log2*reads[1],
                                log1*reads[1] + loge*reads[0]])
        maxVals = np.amax(valSeq, 0)
        valSeq = valSeq - maxVals
        valSeq = np.exp(valSeq)
        vals *= valSeq

    if useSexChrom:
        #Recode so we only care about the two homozygous states, but they are coded as 0, 1.
        vals[1,:] = vals[3,:]
        vals[2,:] = 0
        vals[3,:] = 0

    return vals/np.sum(vals,0)


def getGenotypeProbabilities_log(nLoci, genotypes, reads, error = 0.01, seqError = 0.001, useSexChrom=False):
    vals = np.full((4, nLoci), .25, dtype = np.float32)
    if type(error) is float:
        error = np.full(nLoci, error)
    if type(seqError) is float:
        seqError = np.full(nLoci, seqError)

    errorMat = generateErrorMat(error)

    if genotypes is not None:
        setGenoProbsFromGenotypes(genotypes, errorMat, vals)
    
    vals = np.log(vals)

    if reads is not None:
        log1 = np.log(1-seqError)
        log2 = np.log(.5)
        loge = np.log(seqError)

        ref_reads = reads[0]
        alt_reads = reads[1]

        val_seq = np.full((4, nLoci), 0, dtype = np.float32)
        val_seq[0,:] = log1*ref_reads + loge*alt_reads
        val_seq[1,:] = log2*ref_reads + log2*alt_reads
        val_seq[2,:] = log2*ref_reads + log2*alt_reads
        val_seq[3,:] = loge*ref_reads + log1*alt_reads

        vals += val_seq

    output = np.full((4, nLoci), 0, dtype = np.float32)
    apply_log_norm_1d(vals, output)
    return output

@jit(nopython=True, nogil = True)
def apply_log_norm_1d(vals, output):
    nLoci = vals.shape[-1]
    for i in range(nLoci):
        output[:,i] = log_norm_1D(vals[:, i])


@jit(nopython=True, nogil = True)
def log_norm_1D(mat):
    log_exp_sum = 0
    first = True
    maxVal = 100
    for a in range(4):
        if mat[a] > maxVal or first:
            maxVal = mat[a]
        if first:
            first = False

    for a in range(4):
        log_exp_sum += np.exp(mat[a] - maxVal)
    
    return mat - (np.log(log_exp_sum) + maxVal)


def set_from_genotype_probs(ind, geno_probs = None, calling_threshold = 0.1, set_genotypes = False, set_dosages = False, set_haplotypes = False) :

    # Check diploid geno_probs; not sure what to do for haploid except assume inbred?
    if geno_probs.shape[0] == 2:

        geno_probs = geno_probs/np.sum(geno_probs, axis = 0)
        called_values = call_genotype_probs(geno_probs, calling_threshold)

        # Assuming the individual is haploid

        if set_dosages:
            if ind.dosages is None:
                ind.dosages = called_values.dosages.copy()
            ind.dosages[:] = 2*called_values.dosages

        if set_genotypes:
            ind.genotypes[:] = 2*called_values.haplotypes
            ind.genotypes[called_values.haplotypes == 9] = 9 # Correctly set missing loci.

        if set_haplotypes:
            ind.haplotypes[0][:] = called_values.haplotypes
            ind.haplotypes[1][:] = called_values.haplotypes

    if geno_probs.shape[0] == 4:
        geno_probs = geno_probs/np.sum(geno_probs, axis = 0)
        called_values = call_genotype_probs(geno_probs, calling_threshold)

        if set_dosages:
            if ind.dosages is None:
                ind.dosages = called_values.dosages.copy()
            ind.dosages[:] = called_values.dosages

        if set_genotypes:
            ind.genotypes[:] = called_values.genotypes

        if set_haplotypes:
            ind.haplotypes[0][:] = called_values.haplotypes[0]
            ind.haplotypes[1][:] = called_values.haplotypes[1]


def call_genotype_probs(geno_probs, calling_threshold = 0.1) :

    if geno_probs.shape[0] == 2:
        # Haploid
        HaploidValues = collections.namedtuple("HaploidValues", ["haplotypes", "dosages"])
        dosages = geno_probs[1,:].copy()
        haplotypes = call_matrix(geno_probs, calling_threshold)

        return HaploidValues(dosages = dosages, haplotypes = haplotypes)

    if geno_probs.shape[0] == 4:
        # Diploid
        DiploidValues = collections.namedtuple("DiploidValues", ["genotypes", "haplotypes", "dosages"])

        dosages = geno_probs[1,:] + geno_probs[2,:] + 2*geno_probs[3,:]

        # Collapse the two heterozygous states into one.
        collapsed_hets = np.array([geno_probs[0,:], geno_probs[1,:] + geno_probs[2,:], geno_probs[3,:]], dtype=np.float32)
        genotypes = call_matrix(collapsed_hets, calling_threshold)

        # aa + aA, Aa + AA
        haplotype_0 = np.array([geno_probs[0,:] + geno_probs[1,:], geno_probs[2,:] + geno_probs[3,:]], dtype=np.float32)
        haplotype_1 = np.array([geno_probs[0,:] + geno_probs[2,:], geno_probs[1,:] + geno_probs[3,:]], dtype=np.float32)
        haplotypes = (call_matrix(haplotype_0, calling_threshold), call_matrix(haplotype_1, calling_threshold))
        
        return DiploidValues(dosages = dosages, haplotypes = haplotypes, genotypes = genotypes)


def call_matrix(matrix, threshold):
    called_genotypes = np.argmax(matrix, axis = 0)
    setMissing(called_genotypes, matrix, threshold)
    return called_genotypes.astype(np.int8)

@jit(nopython=True)
def setMissing(calledGenotypes, matrix, threshold) :
    nLoci = len(calledGenotypes)
    for i in range(nLoci):
        if matrix[calledGenotypes[i],i] < threshold:
            calledGenotypes[i] = 9


@jit(nopython=True)
def setGenoProbsFromGenotypes(genotypes, errorMat, vals):
    nLoci = len(genotypes)
    for i in range(nLoci) :
        if genotypes[i] != 9:
            vals[:, i] = errorMat[genotypes[i], :, i]
def generateErrorMat(error) :
    errorMat = np.array([[1-error, error/2, error/2, error/2], 
                            [error/2, 1-error, 1-error, error/2],
                            [error/2, error/2, error/2, 1-error]], dtype = np.float32)
    errorMat = errorMat/np.sum(errorMat, 1)[:,None]
    return errorMat


def generateSegregationXXChrom(partial=False, e= 1e-06) :
    paternalTransmission = np.array([ [1, 1, 0, 0],[0, 0, 1, 1]])
    maternalTransmission = np.array([ [1, 0, 1, 0],[0, 1, 0, 1]])

    fatherAlleleCoding = np.array([0, 0, 1, 1])
    motherAlleleCoding = np.array([0, 1, 0, 1])

    # !                  fm  fm  fm  fm 
    # !segregationOrder: pp, pm, mp, mm
    
    segregationTensor = np.zeros((4, 4, 4, 4))
    for segregation in range(4):
        #Change so that father always passes on the maternal allele?
        if(segregation == 0) :
            father = maternalTransmission
            mother = paternalTransmission
        if(segregation == 1) :
            father = maternalTransmission
            mother = maternalTransmission
        if(segregation == 2) :
            father = maternalTransmission
            mother = paternalTransmission
        if(segregation == 3) :
            father = maternalTransmission
            mother = maternalTransmission

        # !alleles: aa, aA, Aa, AA
        for allele in range(4) :
            segregationTensor[:, :, allele, segregation] = np.outer(father[fatherAlleleCoding[allele]], mother[motherAlleleCoding[allele]])

    segregationTensor = segregationTensor*(1-e) + e/4 #trace has 4 times as many elements as it should since it has 4 internal reps.
    segregationTensor = segregationTensor.astype(np.float32)

    return(segregationTensor)


def generateSegregationXYChrom(partial=False, e= 1e-06) :
    paternalTransmission = np.array([ [1, 1, 0, 0],[0, 0, 1, 1]])
    maternalTransmission = np.array([ [1, 0, 1, 0],[0, 1, 0, 1]])

    motherAlleleCoding = np.array([0, 1, 0, 1])

    # !                  fm  fm  fm  fm 
    # !segregationOrder: pp, pm, mp, mm
    #They don't get anything from the father -- father is always 0
    segregationTensor = np.zeros((4, 4, 4, 4))
    for segregation in range(4):
        if(segregation == 0) :
            mother = paternalTransmission
        if(segregation == 1) :
            mother = maternalTransmission
        if(segregation == 2) :
            mother = paternalTransmission
        if(segregation == 3) :
            mother = maternalTransmission

        # !alleles: aa, aA, Aa, AA
        for allele in range(4) :
            for fatherAllele in range(4):
                segregationTensor[fatherAllele, :, allele, segregation] = mother[motherAlleleCoding[allele]]

    segregationTensor = segregationTensor*(1-e) + e/4 #trace has 4 times as many elements as it should since it has 4 internal reps.
    segregationTensor = segregationTensor.astype(np.float32)

    return(segregationTensor)



def generateSegregation(partial=False, e= 1e-06) :
    paternalTransmission = np.array([ [1, 1, 0, 0],[0, 0, 1, 1]])
    maternalTransmission = np.array([ [1, 0, 1, 0],[0, 1, 0, 1]])

    fatherAlleleCoding = np.array([0, 0, 1, 1])
    motherAlleleCoding = np.array([0, 1, 0, 1])

    # !                  fm  fm  fm  fm 
    # !segregationOrder: pp, pm, mp, mm
    
    segregationTensor = np.zeros((4, 4, 4, 4))
    for segregation in range(4):
        if(segregation == 0) :
            father = paternalTransmission
            mother = paternalTransmission
        if(segregation == 1) :
            father = paternalTransmission
            mother = maternalTransmission
        if(segregation == 2) :
            father = maternalTransmission
            mother = paternalTransmission
        if(segregation == 3) :
            father = maternalTransmission
            mother = maternalTransmission

        # !alleles: aa, aA, Aa, AA
        for allele in range(4) :
            segregationTensor[:, :, allele, segregation] = np.outer(father[fatherAlleleCoding[allele]], mother[motherAlleleCoding[allele]])

    if partial : segregationTensor = np.mean(segregationTensor, 3)

    segregationTensor = segregationTensor*(1-e) + e/4 #trace has 4 times as many elements as it should since it has 4 internal reps.
    segregationTensor = segregationTensor.astype(np.float32)
    return(segregationTensor)

# def generateErrorMat(error) :
#     # errorMat = np.array([[1-error*3/4, error/4, error/4, error/4], 
#     #                         [error/4, .5-error/4, .5-error/4, error/4],
#     #                         [error/4, error/4, error/4, 1-error*3/4]], dtype = np.float32)
#     errorMat = np.array([[1-error*2/3, error/3, error/3, error/3], 
#                             [error/3, 1-error*2/3, 1-error*2/3, error/3],
#                             [error/3, error/3, error/3, 1-error*2/3]], dtype = np.float32)
#     errorMat = errorMat/np.sum(errorMat, 1)[:,None]
#     return errorMat

## Not sure if below is ever used. 
# def generateTransmission(error) :
#     paternalTransmission = np.array([ [1-error, 1.-error, error, error],
#                                       [error, error, 1-error, 1-error]])
    
#     maternalTransmission = np.array([ [1.-error, error, 1.-error, error],
#                                         [error, 1-error, error, 1-error]] )
#     segregationTransmissionMatrix = np.zeros((4,4))
#     segregationTransmissionMatrix[0,:] = paternalTransmission[0,:]
#     segregationTransmissionMatrix[1,:] = paternalTransmission[0,:]
#     segregationTransmissionMatrix[2,:] = paternalTransmission[1,:]
#     segregationTransmissionMatrix[3,:] = paternalTransmission[1,:]

#     segregationTransmissionMatrix[:,0] = segregationTransmissionMatrix[:,0] * maternalTransmission[0,:]
#     segregationTransmissionMatrix[:,1] = segregationTransmissionMatrix[:,1] * maternalTransmission[1,:]
#     segregationTransmissionMatrix[:,2] = segregationTransmissionMatrix[:,2] * maternalTransmission[0,:]
#     segregationTransmissionMatrix[:,3] = segregationTransmissionMatrix[:,3] * maternalTransmission[1,:]

#     segregationTransmissionMatrix = segregationTransmissionMatrix.astype(np.float32)
#     return(segregationTransmissionMatrix)
