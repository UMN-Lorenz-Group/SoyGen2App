from numba import jit
import numpy as np
from . import NumbaUtils
from . HaplotypeLibrary import haplotype_from_indices

def haploidHMM(individual, source_haplotypes, error, recombination_rate, threshold=0.9, calling_method='dosages'):

    target_haplotype = individual.haplotypes
    n_loci = len(target_haplotype)

    # !!!! May need to cast the source Haplotypes to a matrix. #May also want to handle the probabilistic case.
    if type(source_haplotypes) is list:
        source_haplotypes = np.array(source_haplotypes)

    # Expand error and recombinationRate to arrays as may need to have
    # marker specific error/recombination rates.
    if type(error) is float:
        error = np.full(n_loci, error, dtype=np.float32)
    if type(recombination_rate) is float:
        recombination_rate = np.full(n_loci, recombination_rate, dtype=np.float32)

    # Construct penetrance values (point estimates)
    point_estimates = getHaploidPointEstimates(target_haplotype, source_haplotypes, error)

    # Run forward-backward algorithm on penetrance values
    # Note: don't need these probabilites if using the sample method
    if calling_method != 'sample':
        total_probs = haploidForwardBackward(point_estimates, recombination_rate)

    # Handle the different calling methods
    if calling_method == 'callhaps':
        # Call haplotypes
        called_haps = haploidCallHaps(total_probs, threshold)
        # Call genotypes
        called_genotypes = getHaploidGenotypes(called_haps, source_haplotypes)
        return called_genotypes
    if calling_method == 'dosages':
        dosages = getHaploidDosages(total_probs, source_haplotypes)
        individual.dosages = dosages
    if calling_method == 'sample':
        haplotype = getHaploidSample(point_estimates, recombination_rate, source_haplotypes)
        individual.imputed_haplotypes = haplotype
    if calling_method == 'Viterbi':
        haplotype = get_viterbi(point_estimates, recombination_rate, source_haplotypes)
        individual.imputed_haplotypes = haplotype


@jit(nopython=True, nogil=True)
def getHaploidDosages(hap_est, source_haplotypes):
    """Calculate dosages for a single haplotype"""
    n_loci, n_haps = hap_est.shape
    dosages = np.zeros(n_loci, dtype=np.float32)
    for i in range(n_loci):
        for j in range(n_haps):
            dosages[i] += source_haplotypes[j, i] * hap_est[i, j]
    return dosages


@jit(nopython=True, nogil=True)
def getHaploidSample(point_estimates, recombination_rate, source_haps):
    """Sample a haplotype"""
    forward_probs = haploidForward(point_estimates, recombination_rate)
    haplotype = haploidSampleHaplotype(forward_probs, source_haps, recombination_rate)
    return haplotype


@jit(nopython=True, nogil=True)
def get_viterbi(point_estimates, recombination_rate, haplotype_library):
    """Get most likely haplotype using the Viterbi algorithm"""
    forward_probs = haploidForward(point_estimates, recombination_rate)
    indices = haploid_viterbi(forward_probs, recombination_rate)
    return haplotype_from_indices(indices, haplotype_library)


@jit(nopython=True)
def haploidCallHaps(hapEst, threshold):
    nHaps, nLoci = hapEst.shape
    calledHaps = np.full(nLoci, -1, dtype=np.int64)  # These are haplotype ids. -1 is missing.
    for i in range(nLoci):
        maxVal = -1
        maxLoc = -1
        for j in range(nHaps):
            if hapEst[j, i] > threshold and hapEst[j, i] > maxVal:
                maxLoc = j
                maxVal = hapEst[j, i]
        calledHaps[i] = maxLoc
    return calledHaps


@jit(nopython=True)
def getHaploidGenotypes(calledHaps, sourceHaplotypes):
    nHaps, nLoci = sourceHaplotypes.shape
    calledGenotypes = np.full(nLoci, 9, dtype=np.int8) # These are haplotype ids. -1 is missing.
    for i in range(nLoci):
        if calledHaps[i] != -1:
            calledGenotypes[i] = sourceHaplotypes[calledHaps[i], i]
    return calledGenotypes


@jit(nopython=True, nogil=True)
def getHaploidPointEstimates(targetHaplotype, sourceHaplotypes, error):
    nHaps, nLoci = sourceHaplotypes.shape
    pointMat = np.full((nLoci, nHaps), 1, dtype=np.float32)

    for i in range(nLoci):
        if targetHaplotype[i] != 9:
            for j in range(nHaps):
                if targetHaplotype[i] == sourceHaplotypes[j, i]:
                    pointMat[i, j] = 1 - error[i]
                else:
                    pointMat[i, j] = error[i]
    return pointMat


@jit(nopython=True, nogil=True)
def haploidTransformProbs(previous, new, estimate, point_estimate, recombination_rate):
    """Transforms a probability distribution (over haplotypes, at a single locus)
    to a probability distribution at the next locus by accounting for emission probabilities
    (point_estimates) and transition probabilities (recombination_rate)
    This is a core step in the forward and backward algorithms

    point_estimates     emission probabilities - (1D NumPy array)
    recombination_rate  recombination rate at this locus - (scalar)
    previous            probability distribution over haplotypes (hidden states) at the *previous* locus
    estimate            newly calculated probability distribution over haplotypes at *this* locus
    new                 intermediate probability distribution (passed in to this function for speed)

    Note: previous and estimate are updated by this function
    """
    n_haps = len(previous)

    # Get estimate at this locus and normalize
    new[:] = previous * point_estimate
    new /= np.sum(new)

    # Account for recombination rate
    e = recombination_rate
    e1 = 1-recombination_rate
    for j in range(n_haps):
        new[j] = new[j]*e1 + e/n_haps

    # Update distributions (in place)
    for j in range(n_haps):
        estimate[j] *= new[j]
        previous[j] = new[j]


@jit(nopython=True, nogil=True)
def haploidOneSample(forward_probs, recombination_rate):
    """Sample one haplotype (an individual) from the forward and backward probability distributions
    Returns two arrays:
      sample_indices   array of indices of haplotypes in the haplotype library at each locus
                       e.g. an individual composed of haplotypes 13 and 42 with 8 loci:
                       [42, 42, 42, 42, 42, 13, 13, 13]
    A description of the sampling process would be nice here..."""

    est = forward_probs.copy()  # copy so that forward_probs is not modified
    n_loci, n_haps = forward_probs.shape
    prev = np.ones(n_haps, dtype=np.float32)
    new = np.empty(n_haps, dtype=np.float32)

    # Sampled probability distribution at one locus
    sampled_probs = np.empty(n_haps, dtype=np.float32)
    sample_indices = np.empty(n_loci, dtype=np.int64)

    # Backwards algorithm
    for i in range(n_loci-2, -1, -1): # zero indexed then minus one since we skip the boundary
        # Sample at this locus
        j = NumbaUtils.multinomial_sample(pvals=est[i+1, :])
        sampled_probs[:] = 0
        sampled_probs[j] = 1
        sample_indices[i+1] = j

        # Get estimate at this locus using the *sampled* distribution
        # (instead of the point estimates/emission probabilities)
        haploidTransformProbs(prev, new, est[i, :], sampled_probs, recombination_rate[i+1])
        # No need to normalise at this locus as multinomial_sample()
        # handles un-normalized probabilities

    # Last sample (at the first locus)
    j = NumbaUtils.multinomial_sample(pvals=est[0, :])
    sample_indices[0] = j

    return sample_indices


@jit(nopython=True, nogil=True)
def haploid_viterbi(forward_probs, recombination_rate):
    """Find the most likely haplotype according to the The Viterbi algorithm
    Returns:
      indices   array of indices of haplotypes in the haplotype library at each locus
                e.g. an individual composed of haplotypes 13 and 42 with 8 loci:
                [42, 42, 42, 42, 42, 13, 13, 13]"""

    est = forward_probs.copy()  # copy so that forward_probs is not modified
    n_loci, n_haps = forward_probs.shape
    prev = np.ones(n_haps, dtype=np.float32)
    new = np.empty(n_haps, dtype=np.float32)

    # Most likely probability distribution at one locus
    sampled_probs = np.empty(n_haps, dtype=np.float32)
    indices = np.empty(n_loci, dtype=np.int64)

    # Backwards algorithm
    for i in range(n_loci-2, -1, -1): # zero indexed then minus one since we skip the boundary
        # Choose the most likely state (i.e. max probability) at this locus
        j = np.argmax(est[i+1, :])

        sampled_probs[:] = 0
        sampled_probs[j] = 1
        indices[i+1] = j

        # Get estimate at this locus using the most likely distribution
        # (instead of the point estimates/emission probabilities)
        haploidTransformProbs(prev, new, est[i, :], sampled_probs, recombination_rate[i+1])
        # No need to normalise at this locus as argmax() does not depend on normalisation

    # Most likely state at the first locus
    j = np.argmax(est[0, :])
    indices[0] = j

    return indices


@jit(nopython=True, nogil=True)
def haploidSampleHaplotype(forward_probs, haplotype_library, recombination_rate):
    """Sample one haplotype (an individual) from the forward and backward probability distributions
    Returns: a sampled haploytpe of length n_loci"""
    indices = haploidOneSample(forward_probs, recombination_rate)
    return haplotype_from_indices(indices, haplotype_library)


@jit(nopython=True, nogil=True)
def haploidForward(point_estimate, recombination_rate):
    """Calculate (unnomalized) forward probabilities"""

    n_loci, n_haps = point_estimate.shape
    est = point_estimate.copy()
    prev = np.ones(n_haps, dtype=np.float32)
    new = np.empty(n_haps, dtype=np.float32)

    for i in range(1, n_loci):
        # Update estimates at this locus
        haploidTransformProbs(prev, new, est[i, :], point_estimate[i-1, :], recombination_rate[i])

    return est


@jit(nopython=True)
def haploidBackward(point_estimate, recombination_rate):
    """Calculate (unnomalized) backward probabilities"""

    n_loci, n_haps = point_estimate.shape
    est = np.ones_like(point_estimate, dtype=np.float32)
    prev = np.ones(n_haps, dtype=np.float32)
    new = np.empty(n_haps, dtype=np.float32)

    for i in range(n_loci-2, -1, -1):  # zero indexed then minus one since we skip the boundary
        # Update estimates at this locus
        haploidTransformProbs(prev, new, est[i, :], point_estimate[i+1, :], recombination_rate[i+1])

    return est


@jit(nopython=True)
def haploidForwardBackward(point_estimate, recombination_rate):
    """Calculate normalized state probabilities at each loci using the forward-backward algorithm"""

    est = (haploidForward(point_estimate, recombination_rate) *
           haploidBackward(point_estimate, recombination_rate))

    # Return normalized probabilities
    n_loci = point_estimate.shape[0]
    for i in range(n_loci):
        est[i, :] /= np.sum(est[i, :])

    return est
