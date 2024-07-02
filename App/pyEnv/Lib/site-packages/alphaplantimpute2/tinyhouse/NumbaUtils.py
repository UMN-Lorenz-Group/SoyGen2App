from numba import jit
import numpy as np


@jit(nopython=True, nogil=True)
def multinomial_sample(pvals):
    """Draw one sample from a multinomial distribution (similar to np.random.multinomial)
    pvals   - array of probabilities (np.float32) (that ideally should sum to one - but doesn't have to)
    return  - 1d index of the sampled value. Use multinomial_sample_2d() to convert to a 2d index"""

    # Cumulative sum of the probabilities, total sum is cumsum[-1]
    cumsum = np.cumsum(pvals)
    # Random float in the half-open interval [0.0, sum(pvals)) - this handles pvals that don't sum to 1
    rand = np.random.random() * cumsum[-1]
    # Find where the cumulative sum exceeds the random float and return that index
    for index, value in enumerate(cumsum):
        if value > rand:
            break
    return index


@jit(nopython=True, nogil=True)
def multinomial_sample_2d(pvals):
    """Draw one sample from a multinomial distribution (similar to np.random.multinomial)
    pvals   - 2d array of probabilities (np.float32) (that ideally should sum to one - but doesn't have to)
    return  - (row_index, col_index) tuple"""

    # Get 1d index
    index = multinomial_sample(pvals)
    # Convert 1d index to 2d and return indices tuple
    ncols = pvals.shape[1]
    return (index//ncols, index%ncols)
