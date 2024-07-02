import pickle
import random
import numpy as np
from numba import njit, jit

def profile(x):
    return x

# Helper functions
def slices(start, length, n):
    """Return `n` slices starting at `start` of length `length`"""
    return [slice(i, i+length) for i in range(start, length*n + start, length)]

def bin_slices(l, n):
    """Return a list of slice() objects that split l items into n bins
    The first l%n bins are length l//n+1; the remaining n-l%n bins are length l//n
    Similar to np.array_split()"""
    return slices(0, l//n+1, l%n) + slices((l//n+1)*(l%n), l//n, n-l%n)

def topk_indices(genotype, haplotypes, n_topk):
    """Return top-k haplotype indices with fewest opposite homozygous markers compared to genotype
    The array of number of opposite homozygous loci for each haplotype (`opposite_homozygous`)
    is first shuffled so that tied values are effectively randomly sampled"""

    # Mask of homozygous loci in the genotype
    homozygous_mask = (genotype == 0) | (genotype == 2)  # note: this purposefully ignores missing loci

    # Select just these homozygous loci
    g = genotype[homozygous_mask]
    h = haplotypes[:, homozygous_mask]

    # Number of opposite homozygous loci for all haplotypes
    opposite_homozygous = np.sum(g//2 != h, axis=1)

    # Shuffle indices
    indices = np.arange(len(opposite_homozygous))
    np.random.shuffle(indices)

    # Stable argsort on the shuffled values
    args = np.argsort(opposite_homozygous[indices], kind='stable')

    # Top-k indices (fewest opposite homozygous loci)
    return indices[args[:n_topk]]


class HaplotypeLibrary(object):
    """A library of haplotypes
    Each haplotype can have an identifier (any Python object, but typically a str or int)
    The identifiers are used to select haplotypes for updating or masking
    Haplotypes should be NumPy arrays of dtype np.int8

    Some functions only work on frozen libraries; some only on unfrozen ones.
    Use freeze() and unfreeze() to swap between the two states. Typically a library is
    built with append() and then frozen to enable additional functionality"""

    def __init__(self, n_loci=None):
        self._n_loci = n_loci
        self._frozen = False
        self._haplotypes = []
        self._identifiers = []  # index to identifier mapping
        self.dtype = None


    def __repr__(self):
        return repr(self._identifiers) + '\n' + repr(self._haplotypes)


    def __len__(self):
        """Number of haplotypes in the library"""
        return len(self._haplotypes)


    def __iter__(self):
        """Iterate over tuple of (id, haplotype)"""
        for i in range(len(self)):
            yield self._identifiers[i], self._haplotypes[i]


    def append(self, haplotype, identifier=None):
        """Append a single haplotype to the library.
        Note: a copy of the haplotype is taken"""
        if self._frozen:
            raise RuntimeError('Cannot append to frozen library')

        if self.dtype is None:
            self.dtype = haplotype.dtype

        if self._n_loci is None:
            self._n_loci = len(haplotype)

        self._check_haplotype(haplotype, expected_shape=(self._n_loci,))
        self._identifiers.append(identifier)
        self._haplotypes.append(haplotype.copy())


    def freeze(self):
        """Freeze the library: convert identifier and haplotype lists to NumPy arrays"""
        if self._frozen:
            raise RuntimeError('Cannot freeze an already frozen library')
        self._haplotypes = np.array(self._haplotypes)
        self._identifiers = np.array(self._identifiers)
        self._frozen = True


    def unfreeze(self):
        """Unfreeze the library: convert identifiers and haplotypes to lists"""
        if not self._frozen:
            raise RuntimeError('Cannot unfreeze an unfrozen library')
        self._haplotypes = list(self._haplotypes)
        self._identifiers = list(self._identifiers)
        self._frozen = False


    def update(self, haplotypes, identifier):
        """Update identifier's haplotypes
        'haplotypes' can be a 1d array of loci or a 2d array of shape(#haps, #loci)"""
        if not self._frozen:
            raise RuntimeError('Cannot update an unfrozen library')
        self._check_haplotype_dtype(haplotypes)
        indices = self._indices(identifier)
        # Use Numpy's broadcasting checks to handle mismatch of shape in the following:
        self._haplotypes[indices] = haplotypes


    def exclude_identifiers(self, identifiers):
        """Return a NumPy array of haplotypes excluding specified identifiers
        'identifiers' can be a single identifier or iterable of identifiers"""
        if not self._frozen:
            raise RuntimeError('Cannot exclude from an unfrozen library')
        mask = ~np.isin(self._identifiers, identifiers)
        return self._haplotypes[mask]


    def sample_only_identifiers(self, identifiers):
        """Returns HaplotypeLibrary object of haplotypes only including specified identifiers
        'identifiers' can be a single identifier or iterable of identifiers"""
        if not self._frozen:
            raise RuntimeError('Cannot exclude from an unfrozen library')
        mask = np.isin(self._identifiers, identifiers)
        return self._sampled_library(mask)


    def sample(self, n_haplotypes):
        """Return a NumPy array of randomly sampled haplotypes"""
        if not self._frozen:
            raise RuntimeError('Cannot sample from an unfrozen library')
        if n_haplotypes > len(self):
            n_haplotypes = len(self)
        sampled_indices = np.sort(np.random.choice(len(self), size=n_haplotypes, replace=False))
        return self._haplotypes[sampled_indices]


    def sample_targeted(self, n_haplotypes, genotype, n_bins, exclude_identifiers=None):
        """Sample haplotypes that 'closely match' genotype `genotype`"""
        if not self._frozen:
            raise RuntimeError('Cannot sample from an unfrozen library')
        if n_haplotypes > len(self):
            return self

        # Get top-k in a number of marker bins
        n_topk = n_haplotypes  # unnecessary variable redifinition
        indices = np.empty((n_topk, n_bins), dtype=np.int64)
        for i, s in enumerate(bin_slices(self._n_loci, n_bins)):
            indices[:, i] = topk_indices(genotype[s], self._haplotypes[:, s], n_topk)

        # Get top n_haplotypes across the bins excluding any in exclude_ifentifiers
        sampled_indices = set()
        exclude_indices = set(self._indices(exclude_identifiers))
        for idx in indices.flatten():
            if idx not in exclude_indices:
                sampled_indices.add(idx)
            if len(sampled_indices) >= n_topk:
                break
        sampled_indices = list(sampled_indices)
        # Return HaplotypeLibrary object
        return self._sampled_library(sampled_indices)


    def exclude_identifiers_and_sample(self, identifiers, n_haplotypes):
        """Return a NumPy array of (n_haplotypes) randomly sampled haplotypes
        excluding specified identifiers.
        'identifiers' can be a single identifier or an iterable of identifiers
        Note: A copy of the haplotypes are created because of fancy indexing"""
        # Exclude specified identifiers
        if not self._frozen:
            raise RuntimeError('Cannot sample or exclude from an unfrozen library')
        exclude_mask = ~np.isin(self._identifiers, identifiers)
        n_remaining_haplotypes = exclude_mask.sum()
        # Generate random sample
        if n_haplotypes > n_remaining_haplotypes:
            n_haplotypes = n_remaining_haplotypes
        sampled_indices = np.random.choice(n_remaining_haplotypes, size=n_haplotypes, replace=False)
        sampled_indices.sort()
        # Return HaplotypeLibrary object
        return self._sampled_library(sampled_indices)


    def asMatrix(self):
        """Return the NumPy array - kept for backwards compatibility"""
        if self._frozen:
            return self._haplotypes.copy()
        return np.array(self._haplotypes)


    def removeMissingValues(self):
        """Replace missing values randomly with 0 or 1 with 50 % probability
        kept for backwards compatibility"""
        for hap in self._haplotypes:
            removeMissingValues(hap)


    def get_called_haplotypes(self, threshold = 0.99):
        """Return "called" haplotypes -- these are haplotypes which only contain integer values (0,1,9).
        For haplotypes where there is uncertainty, a threshold is used to determine whether the value is called as a value or is missing. """

        if not self._frozen:
            self.freeze()
        if self.dtype is np.int8:
            return self._haplotypes

        else:
            called_haplotypes = np.full(self._haplotypes.shape, 0, dtype = np.float32)
            for i in range(called_haplotypes.shape[0]):
                called_haplotypes[i,:] = self.call_haplotypes(self._haplotypes[i,:], threshold)
            return called_haplotypes

    @staticmethod            
    @jit(nopython=True)
    def call_haplotypes(hap, threshold):
        nLoci = len(hap)
        output = np.full(nLoci, 9, dtype = np.int8)
        for i in range(nLoci):
            if hap[i] <= 1 :
                if hap[i] > threshold : output[i] = 1
                if hap[i] < 1-threshold : output[i] = 0
        return output


    def get_haplotypes(self):
        if not self._frozen:
            self.freeze()

        return self._haplotypes

    def get_identifiers(self):
        """Get haplotype identifiers"""
        return list(self._identifiers)


    def _indices(self, identifier):
        """Get row indices associated with an identifier. These can be used for fancy indexing"""
        # Return empty list if identifier == None
        if not identifier:
            return list()
        if not self._frozen:
            raise RuntimeError("Cannot get indices from an unfrozen library")
        if identifier not in self._identifiers:
            raise KeyError(f"Identifier '{identifier}' not in library")
        return  np.flatnonzero(self._identifiers == identifier).tolist()

    def _check_haplotype_dtype(self, haplotype):
        """Check the haplotype has expected dtype"""
        if haplotype.dtype != self.dtype:
            raise TypeError('haplotype(s) not equal to library dtype, {self.dtype}')

    def _check_haplotype(self, haplotype, expected_shape):
        """Check haplotype has expected shape and dtype.
        Could extend to check values in {0,1,9}"""
        self._check_haplotype_dtype(haplotype)
        if haplotype.shape != expected_shape:
            raise ValueError('haplotype(s) has unexpected shape')

    def _sampled_library(self, indices):
        """Create a 'duplicated' HaplotypeLibrary consisting of specified indices only"""
        # Create HaplotypeLibrary object to return
        library = HaplotypeLibrary(self._n_loci)
        library.dtype = self.dtype
        library._haplotypes = self._haplotypes[indices]
        library._identifiers = self._identifiers[indices]
        library.freeze()
        return library


@jit(nopython=True)
def haplotype_from_indices(indices, haplotype_library):
    """Helper function that takes an array of indices (for each locus) that 'point' to rows
    in a haplotype library (NumPy array) and extracts the alleles from the corresponding haplotypes
    (in the library)
    Returns: a haplotype array of length n_loci"""

    n_loci = len(indices)
    haplotype = np.empty(n_loci, dtype=np.int8)
    for col_idx in range(n_loci):
        row_idx = indices[col_idx]
        haplotype[col_idx] = haplotype_library[row_idx, col_idx]
    return haplotype


@njit
def removeMissingValues(hap):
    for i in range(len(hap)) :
        if hap[i] == 9:
            hap[i] = random.randint(0, 1)

class ErrorLibrary(object):
    @profile
    def __init__(self, hap, haplotypes):
        self.hap = hap
        self.hapLib = haplotypes
        self.errors = jit_assessErrors(hap, self.hapLib)

    def getWindowValue(self, k):
        return jit_getWindowValue(self.errors, k, self.hap)

@jit(nopython=True)
def jit_getWindowValue(errors, k, hap) :
    window = np.full(errors.shape, 0, dtype = np.int8)
    nHaps, nLoci = errors.shape
    #Let me be silly here.
    for i in range(k+1):
        window[:,0] += errors[:,i]

    for i in range(1, nLoci):
        window[:,i] = window[:,i-1]
        if i > k:
            if hap[i-k-1] != 9:
                window[:,i] -= errors[:,i-k-1] #This is no longer in the window.
        if i < (nLoci-k):
            if hap[i+k] != 9:
                window[:,i] += errors[:,i+k] #This is now included in the window.

    return window

@jit(nopython=True)
def jit_assessErrors(hap, haps):
    errors = np.full(haps.shape, 0, dtype = np.int8)
    nHaps, nLoci = haps.shape
    for i in range(nLoci):
        if hap[i] != 9:
            if hap[i] == 0:
                errors[:, i] = haps[:,i]
            if hap[i] == 1:
                errors[:, i] = 1-haps[:,i]
    return errors

from collections import OrderedDict
class HaplotypeDict(object):
    def __init__(self):
        self.nHaps = 0
        self.haps = []
        self.tree = dict()
    # @profile
    def append(self, haplotype):
        byteVal = haplotype.tobytes()
        if byteVal in self.tree:
            return self.tree[byteVal]
        else:
            self.tree[byteVal] = self.nHaps
            self.haps.append(haplotype)
            self.nHaps += 1
            return self.nHaps -1
        return self.tree[byteVal]

    def get(self, index):
        return self.haps[index]

# hap = np.array([0, 0, 0, 0, 0, 0, 0, 0])

# hapLib = [np.array([1, 0, 0, 0, 0, 0, 1]),
#           np.array([0, 1, 0, 0, 0, 1, 0]),
#           np.array([0, 0, 1, 0, 1, 0, 0]),
#           np.array([0, 0, 0, 1, 0, 0, 0]),
#           np.array([0, 0, 1, 0, 1, 0, 0])]

# aa = ErrorLibrary(hap, hapLib)
# print(aa.errors)
# print(aa.getWindowValue(2))


# node_type = deferred_type()
# jit_randomBinary_spec = OrderedDict()
# jit_randomBinary_spec['array'] = int64[:]
# jit_randomBinary_spec['index'] = int64
# jit_randomBinary_spec['nItems'] = int64

# @jitclass(jit_randomBinary_spec)
# class jit_RandomBinary(object):
#     def __init__(self, nItems):
#         self.index = 0
#         self.nItems = nItems
#         self.array = np.random.randint(2, size = nItems)
#     def next():
#         self.index += 1
#         if self.index == self.nItems:
#             self.array = np.random.randint(2, size = nItems)
#             self.index = 0
#         return self.array[self.index]

# I Don't think this is used any more.
# def getCores(nLoci, lengths, offsets = [0]) :
#     nCores = len(lengths)*len(offsets)
#     startStop = []

#     for length in lengths:
#         for offset in offsets:
#             finished = False
#             if offset > 0:
#                 start = 0
#                 stop = min(offset, nLoci)
#                 startStop.append((start, stop))
#                 if stop == nLoci: finished = True
#             else:
#                 stop = 0
#             while not finished:
#                 start = stop
#                 stop = min(stop + length, nLoci)
#                 startStop.append((start, stop))
#                 if stop == nLoci: finished = True
#     return startStop

# from collections import OrderedDict
# class HaplotypeCount_dict(object):
#     def __init__(self):
#         self.nHaps = 0
#         self.haps = []
#         self.tree = OrderedDict()
#     # @profile
#     @profile
#     def append(self, haplotype, score = 1):
#         byteVal = haplotype.tobytes()
#         if byteVal in self.tree:
#             self.tree[byteVal] += score
#         else:
#             self.tree[byteVal] = score
#             self.haps.append(haplotype)
#             self.nHaps += 1

#     def getLargest(self):
#         #This is sloppy
#         vals = [value for key, value in self.tree.items()]

#         index = np.argmax(vals)
#         return self.haps[index]
