from numba import jit
import numpy as np
from . import NumbaUtils
from . import ProbMath

class HaploidMarkovModel :
    def __init__(self, n_loci, error, recombination_rate = None):

        self.update_paramaters(n_loci, error, recombination_rate)

        self.directional_smoothing = self.create_directional_smoothing()
        self.apply_smoothing = self.create_apply_smoothing()
        self.apply_viterbi = self.create_viterbi_algorithm()
        self.apply_sampling = self.create_sampling_algorithm(NumbaUtils.multinomial_sample)


    def update_paramaters(self, n_loci, error, recombination_rate = None):

        self.n_loci = n_loci

        if type(error) is float:
            self.error = np.full(n_loci, error, dtype=np.float32)
        else:
            self.error = error

        if recombination_rate is None:
            recombination_rate = 1.0/n_loci

        if type(recombination_rate) is float:
            self.recombination_rate = np.full(n_loci, recombination_rate, dtype=np.float32)
        else:
            self.recombination_rate = recombination_rate


    def get_mask(self, called_haplotypes):
        return np.all(called_haplotypes != 9, axis = 0)

    # def get_run_option(default_arg, alternative_arg):
    #     # Return the default arg as true if it is supplied, otherwise return the alternative arg.
    #     if default_arg is not None:
    #         if alternative_arg is not None:
    #             if default_arg and alternative_arg:
    #                 print("Both arguments are true, returning default")
    #             if not default_arg and not alternative_arg:
    #                 print("Both arguments are false, returning default")
    #         return default_arg
    #     else:
    #         if alternative_arg is None:
    #             return True
    #         else:
    #             return not alternative_arg


    def run_HMM(self, point_estimates = None, algorithm = "marginalize", **kwargs):
        
        # return_called_values = get_run_option(return_called_values, return_genotype_probabilities)

        if point_estimates is None:
            point_estimates = self.get_point_estimates(**kwargs)
        
        if algorithm == "marginalize":
            total_probs = self.apply_smoothing(point_estimates, self.recombination_rate)
            genotype_probabilities = self.calculate_genotype_probabilities(total_probs, **kwargs)

        elif algorithm == "viterbi":
            total_probs = self.apply_viterbi(point_estimates, self.recombination_rate)
            genotype_probabilities = self.calculate_genotype_probabilities(total_probs, **kwargs)
        
        elif algorithm == "sample":
            total_probs = self.apply_sampling(point_estimates, self.recombination_rate)
            genotype_probabilities = self.calculate_genotype_probabilities(total_probs, **kwargs)

        else:
            print(f"Valid alrogithm option not given: {alrogithm}")

        return genotype_probabilities

    def call_genotype_probabilities(self, genotype_probabilities, threshold = 0.1):

        return ProbMath.call_genotype_probs(genotype_probabilities, threshold)


    def get_point_estimates(self, individual, haplotype_library, library_calling_threshold = 0.95, **kwargs) :
        called_haplotypes = haplotype_library.get_called_haplotypes(threshold = library_calling_threshold)
        mask = self.get_mask(called_haplotypes)     
        point_estimates = self.njit_get_point_estimates(individual.genotypes, called_haplotypes, self.error, mask)
        return point_estimates

    @staticmethod
    @jit(nopython=True, nogil=True)
    def njit_get_point_estimates(genotypes, haplotypes, error, mask):
        nHap, nLoci = haplotypes.shape
        point_estimates = np.full((nLoci, nHap), 1, dtype = np.float32)
        for i in range(nLoci):
            if genotypes[i] != 9 and mask[i]:
                for j in range(nHap):
                    sourceGeno = haplotypes[j, i]
                    if 2*sourceGeno == genotypes[i]:
                        point_estimates[i, j] = 1-error[i]
                    else:
                        point_estimates[i, j] = error[i]

        return point_estimates

    @staticmethod
    @jit(nopython=True, nogil=True)
    def transmission(cummulative_probabilities, previous_point_probability, recombination_rate, output):
        output[:] = cummulative_probabilities * previous_point_probability
        normalize(output)
        output[:] *= (1-recombination_rate)
        output[:] += recombination_rate


    def create_directional_smoothing(self) :
        transmission = self.transmission

        @jit(nopython=True, nogil=True)
        def directional_smoothing(point_estimate, recombination_rate, forward = False, backward = False):
            
            output = np.full(point_estimate.shape, 1, dtype = np.float32)
            n_loci = point_estimate.shape[0]

            if forward:
                start = 1
                stop = n_loci
                step = 1

            if backward:
                start = n_loci - 2
                stop = -1
                step = -1

            for i in range(start, stop, step):
                transmission(output[i-step,:], point_estimate[i - step,:], recombination_rate[i], output[i,:])

            return output

        return directional_smoothing

    def create_apply_smoothing(self):
        directional_smoothing = self.directional_smoothing

        @jit(nopython=True, nogil=True)
        def apply_smoothing(point_estimate, recombination_rate):
            """Calculate normalized state probabilities at each loci using the forward-backward algorithm"""

            est = ( point_estimate * 
                    directional_smoothing(point_estimate, recombination_rate, forward = True) *
                    directional_smoothing(point_estimate, recombination_rate, backward = True) )

            # Return normalized probabilities

            normalize_along_first_axis(est)
            return est

        return apply_smoothing

    def create_sampling_algorithm(self, selection_function):
        directional_smoothing = self.directional_smoothing
        transmission = self.transmission

        @jit(nopython=True, nogil=True)
        def sample_path(point_estimate, recombination_rate):
            """Calculate normalized state probabilities at each loci using the forward-backward algorithm"""

            # Right now using a matrix output; will improve later.
            n_loci = point_estimate.shape[0]

            output = np.full(point_estimate.shape, 0, dtype = np.float32)

            forward_and_point_estimate = point_estimate * directional_smoothing(point_estimate, recombination_rate, forward = True)

            # First index.
            selected_index = selection_function(forward_and_point_estimate[-1].ravel())
            output[- 1].ravel()[selected_index] = 1 # Set the output value at the selected_index to 1.

            # Always sample backward (for tradition mostly).
            locus_estimate = np.full(point_estimate[0].shape, 0, dtype = np.float32)
            matrix_ones = np.full(point_estimate[0].shape, 1, dtype = np.float32)

            start = n_loci - 2
            stop = -1
            step = -1

            for i in range(start, stop, step):
                # Pass along sampled value at the locus.
                transmission(output[i-step,:], matrix_ones, recombination_rate[i], locus_estimate)
                # Combine forward_estimate with backward_estimate
                locus_estimate *= forward_and_point_estimate[i,:]
                selected_index = selection_function(locus_estimate.ravel())
                output[i].ravel()[selected_index] = 1 # Set the output value at the selected_index to 1.

            # Return probabilities
            return output

        return sample_path


    def create_viterbi_algorithm(self):

        maximum_likelihood_step = self.maximum_likelihood_step
        @jit(nopython=True, nogil=True)
        def viterbi_path(point_estimate, recombination_rate):
            """Calculate normalized state probabilities at each loci using the forward-backward algorithm"""

            # Right now using a matrix output; will improve later.
            n_loci = point_estimate.shape[0]

            path_score = np.full(point_estimate.shape, 0, dtype = np.float32)
            previous_index = np.full(point_estimate.shape, 0, dtype = np.int64)

            output = np.full(point_estimate.shape, 0, dtype = np.float32)

            path_score[0] = point_estimate[0]
            start = 1; stop = n_loci; step = 1
            for i in range(start, stop, step):
                # Pass along sampled value at the locus.
                maximum_likelihood_step(path_score[i-step], recombination_rate[i], point_estimate[i], path_score[i], previous_index[i])

            # Traceback
            start_index = np.argmax(path_score[-1])
            output[n_loci-1].ravel()[start_index] = 1

            index = start_index
            start = n_loci-2; stop = -1; step = -1
            for i in range(start, stop, step):
                index = previous_index[i-step].ravel()[index]
                output[i].ravel()[index] = 1
            return output

        return viterbi_path

    @staticmethod
    @jit(nopython=True, nogil=True)
    def maximum_likelihood_step(previous_path_score, recombination_rate, point_estimate, output_path_score, output_index):

        best_index = np.argmax(previous_path_score)
        best_score = previous_path_score[best_index]

        n_hap = previous_path_score.shape[0]
        for i in range(n_hap):

            no_rec_score = (1-recombination_rate)*previous_path_score[i]
            rec_score = best_score*recombination_rate

            if no_rec_score > rec_score:
                # No recombination
                output_path_score[i] = no_rec_score*point_estimate[i]
                output_index[i] = i
            else:
                # Recombination
                output_path_score[i] = rec_score/n_hap*point_estimate[i]
                output_index[i] = best_index
        
        output_path_score /= np.sum(output_path_score)

    def calculate_genotype_probabilities(self, total_probs, haplotype_library, **kwargs):
        haplotype_dosages = haplotype_library.get_haplotypes()
        return self.njit_calculate_genotype_probabilities(total_probs, haplotype_dosages)


    @staticmethod        
    @jit(nopython=True, nogil=True)
    def njit_calculate_genotype_probabilities(total_probs, reference_haplotypes) :
        n_hap, n_loci = reference_haplotypes.shape
        geno_probs = np.full((2, n_loci), 0.0000001, dtype = np.float32) # Adding a very small value as a prior incase all of the values are missing.

        for i in range(n_loci):
            for j in range(n_hap):
                hap_value = reference_haplotypes[j, i]
                prob_value = total_probs[i,j]
                if hap_value != 9:
                    # Add in a sum of total_probs values. 
                    geno_probs[0, i] += prob_value * (1-hap_value)
                    geno_probs[1, i] += prob_value * hap_value

        geno_probs = geno_probs/np.sum(geno_probs, axis = 0)
        return geno_probs


class DiploidMarkovModel(HaploidMarkovModel) :
    def __init__(self, n_loci, error, recombination_rate = None):
        HaploidMarkovModel.__init__(self, n_loci, error, recombination_rate)


    def extract_reference_panels(self, haplotype_library = None, paternal_haplotype_library = None, maternal_haplotype_library = None) :
        if maternal_haplotype_library is not None and paternal_haplotype_library is not None:
            seperate_reference_panels = True
            return paternal_haplotype_library, maternal_haplotype_library, seperate_reference_panels

        else:
            seperate_reference_panels = False
            return haplotype_library, haplotype_library, seperate_reference_panels


    def get_point_estimates(self, individual, library_calling_threshold= 0.95, **kwargs):
        paternal_haplotype_library, maternal_haplotype_library, seperate_reference_panels = self.extract_reference_panels(**kwargs)

        paternal_called_haplotypes = paternal_haplotype_library.get_called_haplotypes(threshold = library_calling_threshold)
        maternal_called_haplotypes = maternal_haplotype_library.get_called_haplotypes(threshold = library_calling_threshold)

        mask = self.get_mask(paternal_called_haplotypes) & self.get_mask(maternal_called_haplotypes) 
        return self.njit_get_point_estimates(individual.genotypes, paternal_called_haplotypes, maternal_called_haplotypes, self.error, mask)


    @staticmethod
    @jit(nopython=True, nogil=True)
    def njit_get_point_estimates(indGeno, paternalHaplotypes, maternalHaplotypes, error, mask):
        nPat, nLoci = paternalHaplotypes.shape
        nMat, nLoci = maternalHaplotypes.shape

        point_estimates = np.full((nLoci, nPat, nMat), 1, dtype = np.float32)

        for i in range(nLoci):
            if indGeno[i] != 9 and mask[i]:
                for j in range(nPat):
                    for k in range(nMat):
                        sourceGeno = paternalHaplotypes[j, i] + maternalHaplotypes[k, i]
                        if sourceGeno == indGeno[i]:
                            point_estimates[i, j, k] = 1-error[i]
                        else:
                            point_estimates[i, j, k] = error[i]

        return point_estimates


    def calculate_genotype_probabilities(self, total_probs, haplotype_library = None, paternal_haplotype_library= None, maternal_haplotype_library= None, **kwargs):
        paternal_haplotype_library, maternal_haplotype_library, seperate_reference_panels = self.extract_reference_panels(haplotype_library, paternal_haplotype_library, maternal_haplotype_library)
        return self.njit_calculate_genotype_probabilities(total_probs, paternal_haplotype_library.get_haplotypes(), maternal_haplotype_library.get_haplotypes(), seperate_reference_panels)

    @staticmethod
    @jit(nopython=True, nogil=True)
    def njit_calculate_genotype_probabilities(total_probs, paternal_haplotypes, maternal_haplotypes, seperate_reference_panels) :
        n_pat, n_loci = paternal_haplotypes.shape
        n_mat, n_loci = maternal_haplotypes.shape
        geno_probs = np.full((4, n_loci), 0.00001, dtype = np.float32)

        for i in range(n_loci):
            for j in range(n_pat):
                for k in range(n_mat):
                    # diploid case where the markers are assumed independent.
                    if seperate_reference_panels or j != k: 
                        pat_value = paternal_haplotypes[j, i]
                        mat_value = maternal_haplotypes[k, i]
                        prob_value = total_probs[i,j,k]
                        if pat_value != 9 and mat_value != 9:
                            # Add in a sum of total_probs values. 
                            geno_probs[0, i] += prob_value * (1-pat_value)*(1-mat_value)   #aa
                            geno_probs[1, i] += prob_value * (1-pat_value)*mat_value     #aA
                            geno_probs[2, i] += prob_value * pat_value*(1-mat_value)     #Aa
                            geno_probs[3, i] += prob_value * pat_value*mat_value       #AA

                    # Haploid case where the markers are not independent
                    else:
                        hap_value = paternal_haplotypes[j, i]
                        prob_value = total_probs[i,j,k]
                        if hap_value != 9:
                            geno_probs[0, i] += prob_value * (1-hap_value)
                            geno_probs[1, i] += 0
                            geno_probs[2, i] += 0
                            geno_probs[3, i] += prob_value * hap_value

        geno_probs = geno_probs/np.sum(geno_probs, axis = 0)
        return geno_probs


    
    
    @staticmethod
    @jit(nopython=True, nogil=True)
    def transmission(cummulative_probabilities, previous_point_probability, recombination_rate, output):
        output[:] = cummulative_probabilities * previous_point_probability
        normalize(output)

        row_sums = np.sum(output, 0)
        col_sums = np.sum(output, 1)

        output[:] *= (1 - recombination_rate)**2 # No recombination on either chromosome.
        output[:] += np.expand_dims(row_sums, 0)/output.shape[0]*recombination_rate*(1-recombination_rate) # recombination on the maternal (second) chromosome)
        output[:] += np.expand_dims(col_sums, 1)/output.shape[1]*recombination_rate*(1-recombination_rate) # recombination on the paternal (first) chromosome)
        output[:] += recombination_rate**2/output.size # double recombination
    
    @staticmethod
    @jit(nopython=True, nogil=True)
    def maximum_likelihood_step(previous_path_score, recombination_rate, point_estimate, output_path_score, output_index):

        n_pat = previous_path_score.shape[0]
        n_mat = previous_path_score.shape[1]

        combined_max_index = np.argmax(previous_path_score)
        combined_max_score = previous_path_score.ravel()[combined_max_index] * recombination_rate**2/(n_mat*n_pat)

        paternal_max_index = np.full(n_pat, 0, dtype = np.int64)
        paternal_max_value = np.full(n_pat, 0, dtype = np.float32)

        maternal_max_index = np.full(n_mat, 0, dtype = np.int64)
        maternal_max_value = np.full(n_mat, 0, dtype = np.float32)

        # Recombination on the maternal side, paternal side is fixed
        for i in range(n_pat):
            index = np.argmax(previous_path_score[i,:])
            paternal_max_value[i] = previous_path_score[i, index] * (1-recombination_rate)*recombination_rate/n_mat
            paternal_max_index[i] = i*n_mat + index


        # Recombination on the paternal side, maternal side is fixed
        for j in range(n_mat):
            index = np.argmax(previous_path_score[:, j])
            maternal_max_value[j] = previous_path_score[index, j] * (1-recombination_rate)*recombination_rate/n_pat
            maternal_max_index[j] = index*n_mat + j

        for i in range(n_pat):
            for j in range(n_mat):

                best_score = (1-recombination_rate)**2*previous_path_score[i,j]
                best_index = i*n_mat + j

                # Paternal recombination 
                if paternal_max_value[i] > best_score:
                    best_score = paternal_max_value[i]
                    best_index = paternal_max_index[i]

                if maternal_max_value[j] > best_score:
                    best_score = maternal_max_value[j]
                    best_index = maternal_max_index[j]
                
                if combined_max_score > best_score:
                    best_score = combined_max_score
                    best_index = combined_max_index

                output_path_score[i,j] = best_score*point_estimate[i,j]
                output_index[i,j] = best_index
        
 
        output_path_score /= np.sum(output_path_score)

class JointMarkovModel(HaploidMarkovModel) :
    def __init__(self, n_loci, error, recombination_rate = None):
        HaploidMarkovModel.__init__(self, n_loci, error, recombination_rate)


    @staticmethod
    @jit(nopython=True, nogil=True)
    def njit_get_point_estimates(indGeno, haplotypes, error, mask):
        n_hap, n_loci = haplotypes.shape

        point_estimates = np.full((n_loci, n_hap, n_hap + 1), 1, dtype = np.float32)

        diploid_section = point_estimates[:,:,0:-1]
        haploid_section = point_estimates[:,:,-1]

        # Diploid point estimates/Emission probabilities

        for i in range(n_loci):
            if indGeno[i] != 9 and mask[i]:
                for j in range(n_hap):
                    for k in range(n_hap):
                        sourceGeno = haplotypes[j, i] + haplotypes[k, i]
                        if sourceGeno == indGeno[i]:
                            diploid_section[i, j, k] = 1-error[i]
                        else:
                            diploid_section[i, j, k] = error[i]

        # Diploid point estimates/Emission probabilities
        for i in range(n_loci):
            if indGeno[i] != 9 and mask[i]:
                for j in range(n_hap):
                        sourceGeno = 2*haplotypes[j, i]
                        if sourceGeno == indGeno[i]:
                            haploid_section[i, j] = 1-error[i]
                        else:
                            haploid_section[i, j] = error[i]

        return point_estimates


    @staticmethod
    @jit(nopython=True, nogil=True)
    def njit_calculate_genotype_probabilities(total_probs, reference_haplotypes) :
        n_hap, n_loci = reference_haplotypes.shape
        geno_probs = np.full((4, n_loci), 0.00001, dtype = np.float32)

        diploid_section = total_probs[:,:,0:-1]
        haploid_section = total_probs[:,:,-1]

        for i in range(n_loci):
            for j in range(n_hap):
                for k in range(n_hap):
                    # diploid case where the markers are assumed independent.
                    if j != k: 
                        pat_value = reference_haplotypes[j, i]
                        mat_value = reference_haplotypes[k, i]
                        prob_value = diploid_section[i,j,k]
                        if pat_value != 9 and mat_value != 9:
                            # Add in a sum of total_probs values. 
                            geno_probs[0, i] += prob_value * (1-pat_value)*(1-mat_value)   #aa
                            geno_probs[1, i] += prob_value * (1-pat_value)*mat_value     #aA
                            geno_probs[2, i] += prob_value * pat_value*(1-mat_value)     #Aa
                            geno_probs[3, i] += prob_value * pat_value*mat_value       #AA

                    # markers are not independent
                    else:
                        hap_value = reference_haplotypes[j, i]
                        prob_value = diploid_section[i,j,k]
                        if hap_value != 9:
                            geno_probs[0, i] += prob_value * (1-hap_value)
                            geno_probs[1, i] += 0
                            geno_probs[2, i] += 0
                            geno_probs[3, i] += prob_value * hap_value

        for i in range(n_loci):
            for j in range(n_hap):
                hap_value = reference_haplotypes[j, i]
                prob_value = haploid_section[i,j]
                if hap_value != 9:
                    geno_probs[0, i] += prob_value * (1-hap_value)
                    geno_probs[1, i] += 0
                    geno_probs[2, i] += 0
                    geno_probs[3, i] += prob_value * hap_value

        geno_probs = geno_probs/np.sum(geno_probs, axis = 0)
        return geno_probs


    
    
    @staticmethod
    @jit(nopython=True, nogil=True)
    def transmission(cummulative_probabilities, previous_point_probability, recombination_rate, output):
        
        output[:] = cummulative_probabilities * previous_point_probability
        normalize(output)

        diploid_section = output[:,0:-1]
        haploid_section = output[:,-1]

        diploid_weight = np.sum(diploid_section)
        haploid_weight = np.sum(haploid_section)

        row_sums = np.sum(diploid_section, 0)
        col_sums = np.sum(diploid_section, 1)

        diploid_section[:] *= (1 - recombination_rate)**2
        diploid_section[:] += np.expand_dims(row_sums, 0)/diploid_section.shape[0]*recombination_rate*(1-recombination_rate) # recombination on the maternal (second) chromosome)
        diploid_section[:] += np.expand_dims(col_sums, 1)/diploid_section.shape[1]*recombination_rate*(1-recombination_rate) # recombination on the paternal (first) chromosome)
        diploid_section[:] += diploid_weight*recombination_rate**2/diploid_section.size # double recombination

        haploid_section[:] *= (1 - recombination_rate)
        haploid_section[:] += haploid_weight*recombination_rate/haploid_section.size

        # loose the recombination to the haploid section and add the haploid recombination to diploid
        diploid_section[:] *= (1 - recombination_rate)
        diploid_section[:] += recombination_rate * haploid_weight/diploid_section.size

        # loose the recombination to the haploid section and add the haploid recombination to diploid
        haploid_section[:] *= (1 - recombination_rate)
        haploid_section[:] += recombination_rate * diploid_weight/haploid_section.size

    
    # @staticmethod
    # @jit(nopython=True, nogil=True)
    # def maximum_likelihood_step(previous_path_score, recombination_rate, point_estimate, output_path_score, output_index):

    #     n_pat = previous_path_score.shape[0]
    #     n_mat = previous_path_score.shape[1]

    #     combined_max_index = np.argmax(previous_path_score)
    #     combined_max_score = previous_path_score.ravel()[combined_max_index] * recombination_rate**2/(n_mat*n_pat)

    #     paternal_max_index = np.full(n_pat, 0, dtype = np.int64)
    #     paternal_max_value = np.full(n_pat, 0, dtype = np.float32)

    #     maternal_max_index = np.full(n_mat, 0, dtype = np.int64)
    #     maternal_max_value = np.full(n_mat, 0, dtype = np.float32)

    #     # Recombination on the maternal side, paternal side is fixed
    #     for i in range(n_pat):
    #         index = np.argmax(previous_path_score[i,:])
    #         paternal_max_value[i] = previous_path_score[i, index] * (1-recombination_rate)*recombination_rate/n_mat
    #         paternal_max_index[i] = i*n_mat + index


    #     # Recombination on the paternal side, maternal side is fixed
    #     for j in range(n_mat):
    #         index = np.argmax(previous_path_score[:, j])
    #         maternal_max_value[j] = previous_path_score[index, j] * (1-recombination_rate)*recombination_rate/n_pat
    #         maternal_max_index[j] = index*n_mat + j

    #     for i in range(n_pat):
    #         for j in range(n_mat):

    #             best_score = (1-recombination_rate)**2*previous_path_score[i,j]
    #             best_index = i*n_mat + j

    #             # Paternal recombination 
    #             if paternal_max_value[i] > best_score:
    #                 best_score = paternal_max_value[i]
    #                 best_index = paternal_max_index[i]

    #             if maternal_max_value[j] > best_score:
    #                 best_score = maternal_max_value[j]
    #                 best_index = maternal_max_index[j]
                
    #             if combined_max_score > best_score:
    #                 best_score = combined_max_score
    #                 best_index = combined_max_index

    #             output_path_score[i,j] = best_score*point_estimate[i,j]
    #             output_index[i,j] = best_index
        
 
    #     output_path_score /= np.sum(output_path_score)


@jit(nopython=True, nogil=True)
def normalize(mat):
    mat[:] /= np.sum(mat)    

@jit(nopython=True, nogil=True)
def normalize_along_first_axis(mat):
    for i in range(mat.shape[0]):
        normalize(mat[i,:])


