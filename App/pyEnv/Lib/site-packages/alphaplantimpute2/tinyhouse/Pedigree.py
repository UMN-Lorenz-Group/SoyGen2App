import numpy as np
import numba
import sys

try:
    from numba.experimental import jitclass
except ModuleNotFoundError:
    from numba import jitclass

from numba import jit, int8, int64, boolean, deferred_type, optional, float32, double
from collections import OrderedDict
from . import InputOutput
from . import ProbMath
from . import MultiThreadIO




class Family(object):
    """Family is a container for fullsib families"""
    def __init__(self, idn, sire, dam, offspring):
        self.idn = idn
        self.sire = sire
        self.dam = dam
        self.offspring = offspring
        self.generation = max(sire.generation, dam.generation) + 1

        # Add this family to both the sire and dam's family.
        self.sire.families.append(self)
        self.dam.families.append(self)

    def addChild(self, child) :
        self.offspring.append(child)

    def toJit(self):
        """Returns a just in time version of itself with Individuals replaced by id numbers"""
        offspring = np.array([child.idn for child in self.offspring])
        return jit_Family(self.idn, self.sire.idn, self.dam.idn, offspring)


spec = OrderedDict()
spec['idn'] = int64
spec['sire'] = int64
spec['dam'] = int64
spec['offspring'] = int64[:]

@jitclass(spec)
class jit_Family(object):
    def __init__(self, idn, sire, dam, offspring):
        self.idn = idn
        self.sire = sire
        self.dam = dam
        self.offspring = offspring

class Individual(object):
    
    def __init__(self, idx, idn) :

        self.genotypes = None
        self.haplotypes = None
        self.dosages = None

        self.imputed_genotypes = None
        self.imputed_haplotypes = None

        self.reads = None
        self.longReads = []

        # Do we use this?
        self.genotypeDosages = None
        self.haplotypeDosages = None
        self.hapOfOrigin = None

        # Info is here to provide other software to add in additional information to an Individual.
        self.info = None

        #For plant impute. Inbred is either DH or heavily selfed. Ancestors is historical source of the cross (may be more than 2 way so can't handle via pedigree).
        self.inbred = False
        self.imputationAncestors = [] #This is a list of lists. Either length 0, length 1 or length 2.
        self.selfingGeneration = None


        self.sire = None
        self.dam = None
        self.idx = idx # User inputed string identifier
        self.idn = idn # ID number assigned by the pedigree
        self.fileIndex = dict() # Position of an animal in each file when reading in. Used to make sure Input and Output order are the same.
        self.fileIndex["id"] = idx

        self.dummy = False

        self.offspring = []
        self.families = []

        self.sex = -1
        self.generation = None

        self.initHD = False
        # used in pythonHMM, but how should this best be coded when already have initHD? Should probably be set when data is read in,
        # but need high_density_threshold to be set in the pedigree first
        self.is_high_density = False

        self.genotypedFounderStatus = None #?
    
    def __eq__(self, other):
        return self is other

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.idn)

    def subset(self, start, stop):
        # Take an individual and create an individual that just contains information on those markers.
        # It's okay if this wipes other information.

        # print(self.__class__)
        new_ind = self.__class__(self.idx, self.idn)
        if self.genotypes is not None:
            new_ind.genotypes = self.genotypes[start:stop].copy() # Maybe could get away with not doing copies... doing them just to be safe.

        if self.haplotypes is not None:
            new_ind.haplotypes = (self.haplotypes[0][start:stop].copy(), self.haplotypes[1][start:stop].copy())

        if self.reads is not None:
            new_ind.reads = (self.reads[0][start:stop].copy(), self.reads[1][start:stop].copy())
        return new_ind

    def getPercentMissing(self):
        return np.mean(self.genotypes == 9)

    def getGeneration(self):

        if self.generation is not None : return self.generation

        if self.dam is None: 
            damGen = -1
        else:
            damGen = self.dam.getGeneration()
        if self.sire is None: 
            sireGen = -1
        else:
            sireGen = self.sire.getGeneration()

        self.generation = max(sireGen, damGen) + 1
        return self.generation


    def constructInfo(self, nLoci, genotypes = True,  haps = False, reads = False) :
        if genotypes and self.genotypes is None:
            self.genotypes = np.full(nLoci, 9, dtype = np.int8)
        
        if  haps and self.haplotypes is None:
            self.haplotypes = (np.full(nLoci, 9, dtype = np.int8), np.full(nLoci, 9, dtype = np.int8))

        if reads and self.reads is None:
            self.reads = (np.full(nLoci, 0, dtype = np.int64), np.full(nLoci, 0, dtype = np.int64))
    
    def isFounder(self):
        return (self.sire is None) and (self.dam is None)

    def getGenotypedFounderStatus(self):
        # options: 1: "GenotypedFounder", 0:"ChildOfNonGenotyped", 2:"ChildOfGenotyped"
        if self.genotypedFounderStatus is None:
            if self.isFounder() :
                if self.genotypes is None or np.all(self.genotypes == 9):
                    self.genotypedFounderStatus = 0 
                else:
                    self.genotypedFounderStatus = 1
            else:
                parentStatus = max(self.sire.getGenotypedFounderStatus(), self.sire.getGenotypedFounderStatus())
                if parentStatus > 0:
                    self.genotypedFounderStatus = 2
                else:
                    if self.genotypes is None or np.all(self.genotypes == 9):
                        self.genotypedFounderStatus = 0 
                    else:
                        self.genotypedFounderStatus = 1

        return self.genotypedFounderStatus
    def isGenotypedFounder(self):
        return (self.getGenotypedFounderStatus() == 1)


class PlantImputeIndividual(Individual):
    """Simple derived class for AlphaPlantImpute2
    with some extra member variables"""
    def __init__(self, idx, idn):
        super().__init__(idx, idn)
        self.founders = []
        self.descendants = []


# Not sure of the code source: https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/
# Slightly modified.
import re 
def sorted_nicely( l , key): 
    """ Sort the given iterable in the way that humans expect.""" 
    return sorted(l, key = lambda k: alphanum_key(key(k)))

def alphanum_key(k):
    convert = lambda text: int(text) if text.isdigit() else text 
    return [ convert(c) for c in re.split('([0-9]+)', str(k)) ] 


class Generation(object):

    def __init__(self, number):

        self.number = number

        self.families = []

        self.individuals = []
        self.sires = set()
        self.dams = set()
        self.parents = set()

    def add_individual(self, ind):
        self.individuals.append(ind)

    def add_family(self, fam):
        self.families.append(fam)
        self.sires.add(fam.sire)
        self.dams.add(fam.dam)
        self.parents.add(fam.sire)
        self.parents.add(fam.dam)
        # Note: Individuals are added seperately.



class Pedigree(object):
 
    def __init__(self, fileName = None, constructor = Individual):

        self.maxIdn = 0
        self.maxFam = 0

        self.individuals = dict()
        self.families = None
        self.constructor = constructor
        self.nGenerations = 0
        self.generations = None #List of lists

        self.truePed = None
        self.nLoci = 0
        
        self.startsnp = 0
        self.endsnp = self.nLoci

        self.referencePanel = [] #This should be an array of haplotypes. Or a dictionary?

        self.maf=None #Maf is the frequency of 2s.

        # Threshold that determines if an individual is high-density genotyped
        self.high_density_threshold = 0.9

        if fileName is not None:
            self.readInPedigree(fileName)

        self.args = None
        self.writeOrderList = None

        self.allele_coding = None

    def reset_families(self):

        for ind in self.individuals.values():
            ind.families = []
            ind.generation = None

        self.nGenerations = 0
        self.generations = None
        self.families = None
        for ind in self:
            ind.families = []

        self.setupFamilies()

    def subset(self, start, stop):
        new_pedigree = Pedigree(constructor = self.constructor)
        new_pedigree.nLoci = stop - start
        # Add all of the individuals.

        for ind in self:
            # Note: ind.subset strips away all of the family information.
            new_pedigree[ind.idx] = ind.subset(start, stop)

        for ind in self:
            if ind.sire is not None:
                new_ind = new_pedigree[ind.idx]
                new_sire = new_pedigree[ind.sire.idx]


                # Add individuals
                new_ind.sire = new_sire
                new_sire.offspring.append(new_ind)

            if ind.dam is not None:
                new_ind = new_pedigree[ind.idx]
                new_dam = new_pedigree[ind.dam.idx]

                # Add individuals
                new_ind.dam = new_dam
                new_dam.offspring.append(new_ind)
        return new_pedigree

    def merge(self, new_pedigree, start, stop):
        # This just merged genotype, haplotype, and dosage information (if availible).
        
        for ind in self:
            new_ind = new_pedigree[ind.idx]
            if new_ind.genotypes is not None:
                if ind.genotypes is None:
                    ind.genotypes = np.full(self.nLoci, 9, dtype = np.int8)
                ind.genotypes[start:stop] = new_ind.genotypes
            
            if new_ind.dosages is not None:
                if ind.dosages is None:
                    ind.dosages = np.full(self.nLoci, -1, dtype = np.float32)

                ind.dosages[start:stop] = new_ind.dosages
            
            if new_ind.haplotypes is not None:
                if ind.haplotypes is None:
                    ind.haplotypes = (np.full(self.nLoci, 9, dtype = np.int8), np.full(self.nLoci, 9, dtype = np.int8))

                ind.haplotypes[0][start:stop] = new_ind.haplotypes[0]
                ind.haplotypes[1][start:stop] = new_ind.haplotypes[1]



    def __len__(self):
        return len(self.individuals)


    def writeOrder(self):
        if self.writeOrderList is None:
            inds = [ind for ind in self if (not ind.dummy) and (self.args.writekey in ind.fileIndex)]
            self.writeOrderList = sorted_nicely(inds, key = lambda ind: ind.fileIndex[self.args.writekey])
            
            if not self.args.onlykeyed:
                indsNoDummyNoFileIndex = [ind for ind in self if (not ind.dummy) and (not self.args.writekey in ind.fileIndex)]
                self.writeOrderList.extend(sorted_nicely(indsNoDummyNoFileIndex, key = lambda ind: ind.idx))
                
                dummys = [ind for ind in self if ind.dummy]
                self.writeOrderList.extend(sorted_nicely(dummys, key = lambda ind: ind.idx))

        for ind in self.writeOrderList :
            yield (ind.idx, ind)

    def setMaf(self):
        """Calculate minor allele frequencies at each locus"""

        # The default values of 1 (maf) and 2 (counts) provide a sensible prior
        # For example, a locus where all individuals are missing, the MAF will be 0.5
        maf = np.full(self.nLoci, 1, dtype = np.float32)
        counts = np.full(self.nLoci, 2, dtype = np.float32)
        for ind in self.individuals.values():
            if ind.genotypes is not None:
                addIfNotMissing(maf, counts, ind.genotypes)
        
        self.maf = maf/counts

    def getMissingness(self):
        missingness = np.full(self.nLoci, 1, dtype = np.float32)

        counts = 0
        for ind in self.individuals.values():
            if ind.genotypes is not None:
                counts += 1
                addIfMissing(missingness, ind.genotypes)
        return missingness/counts


    def set_high_density(self):
        """Set whether each individual is high-density"""
        for individual in self:#.individuals.values():
            is_high_density = np.mean(individual.genotypes != 9) >= self.high_density_threshold
            if is_high_density:
                individual.is_high_density = True


    def fillIn(self, genotypes = True, haps = False, reads = False):

        for individual in self:
            individual.constructInfo(self.nLoci, genotypes = True, haps = haps, reads = reads)


    def __getitem__(self, key) :
        return self.individuals[key]

    def __setitem__(self, key, value):
        self.individuals[key] = value

    def __iter__(self) :
        if self.generations is None:
            self.setUpGenerations()
        for gen in self.generations:
            for ind in gen.individuals:
                yield ind

    def __reversed__(self) :
        if self.generations is None:
            self.setUpGenerations()
        for gen in reversed(self.generations):
            for ind in gen.individuals:
                yield ind

    def sort_individuals(self, individuals):

        return {k:v for k, v in sorted(individuals.items(), key = lambda pair: alphanum_key(pair[0]))}

    # Generation code

    def setUpGenerations(self) :
        # Try and make pedigree independent.
        self.individuals = self.sort_individuals(self.individuals)

        self.nGenerations = 0
        #We can't use a simple iterator over self here, becuase __iter__ calls this function.
        for idx, ind in self.individuals.items():
            gen = ind.getGeneration()
            self.nGenerations = max(gen, self.nGenerations)
        self.nGenerations += 1 #To account for generation 0. 

        self.generations = [Generation(i) for i in range(self.nGenerations)]

        for idx, ind in self.individuals.items():
            gen = ind.getGeneration()
            self.generations[gen].add_individual(ind)

    #This is really sloppy, but probably not important.
    def setupFamilies(self) :

        if self.generations is None:
            self.setUpGenerations()

        self.families = dict()
        for ind in self:
            if not ind.isFounder():
                parents = (ind.sire.idx, ind.dam.idx)
                if parents in self.families :
                    self.families[parents].addChild(ind)
                else:
                    self.families[parents] = Family(self.maxFam, ind.sire, ind.dam, [ind])
                    self.maxFam += 1

        
        for family in self.families.values():
            self.generations[family.generation].add_family(family)
    

    def getFamilies(self, rev = False) :
        if self.generations is None:
            self.setUpGenerations()
        if self.families is None:
            self.setupFamilies()

        gens = range(0, len(self.generations))
        if rev: gens = reversed(gens)

        for i in gens:
            for family in self.generations[i].families:
                yield family
 


    def getIndividual(self, idx) :
        if idx not in self.individuals:
            self.individuals[idx] = self.constructor(idx, self.maxIdn)
            self.maxIdn += 1
            self.generations = None
        return self.individuals[idx]

    def readInPedigree(self, fileName):
        with open(fileName) as f:
            lines = f.readlines()
        pedList = [line.split() for line in lines]
        self.readInPedigreeFromList(pedList)

    def readInPlantInfo(self, fileName):
        with open(fileName) as f:
            lines = f.readlines()

        for line in lines:
            parts = line.split()
            idx = parts[0]; 

            if idx not in self.individuals:
                self.individuals[idx] = self.constructor(idx, self.maxIdn)
                self.maxIdn += 1

            ind = self.individuals[idx]
            if len(parts) > 1:
                if parts[1] == "DH" or parts[1] == "INBRED":
                    ind.inbred = True
                elif parts[1][0] == "S" :
                    ind.inbred = False
                    ind.selfingGeneration = int(parts[1][1:])
                else:
                    ind.inbred = False

            if len(parts) > 2:
                if "|" in line:
                    first, second = line.split("|")
                    self.addAncestors(ind, first.split()[2:])
                    self.addAncestors(ind, second.split())
                else:
                    self.addAncestors(ind, parts[2:])


    def addAncestors(self, ind, parts):
        ancestors = []
        for idx in parts:
            if idx not in self.individuals:
                self.individuals[idx] = self.constructor(idx, self.maxIdn)
                self.maxIdn += 1
            ancestor = self.individuals[idx]
            ancestors.append(ancestor)
        ind.imputationAncestors.append(ancestors)


    def readInPedigreeFromList(self, pedList):
        index = 0
        for parts in pedList :
            idx = parts[0]
            self.individuals[idx] = self.constructor(idx, self.maxIdn)
            self.maxIdn += 1
            self.individuals[idx].fileIndex['pedigree'] = index; index += 1

        for parts in pedList :
            idx = parts[0]
            if parts[1] == "0": parts[1] = None
            if parts[2] == "0": parts[2] = None
            
            if parts[1] is not None and parts[2] is None:
                parts[2] = "MotherOf"+parts[0]
            if parts[2] is not None and parts[1] is None:
                parts[1] = "FatherOf"+parts[0] 

            ind = self.individuals[parts[0]]
            
            if parts[1] is not None:
                if parts[1] not in self.individuals:
                    self.individuals[parts[1]] = self.constructor(parts[1], self.maxIdn)
                    self.maxIdn += 1
                    self.individuals[parts[1]].fileIndex['pedigree'] = index; index += 1
                    self.individuals[parts[1]].dummy=True

                sire = self.individuals[parts[1]]
                ind.sire = sire
                sire.offspring.append(ind)
                sire.sex = 0

            if parts[2] is not None:
                if parts[2] not in self.individuals:
                    self.individuals[parts[2]] = self.constructor(parts[2], self.maxIdn)
                    self.maxIdn += 1
                    self.individuals[parts[2]].fileIndex['pedigree'] = index; index += 1
                    self.individuals[parts[1]].dummy=True

                dam = self.individuals[parts[2]]
                ind.dam = dam
                dam.offspring.append(ind)
                dam.sex = 1

            # Optional fourth column contains sex OR inbred/outbred status
            if len(parts) > 3:
                male, female = {'m', '0', 'xy'}, {'f', '1', 'xx'}
                inbred, outbred = {'dh', 'inbred'}, {'outbred'}
                expected_entries = male | female | inbred | outbred
                if parts[3].lower() not in expected_entries:
                    print(f"ERROR: unexpected entry in pedigree file, fourth field: '{parts[3]}'\nExiting...")
                    sys.exit(2)
                # Sex
                if parts[3].lower() in male:
                    ind.sex = 0
                elif parts[3].lower() in female:
                    ind.sex = 1
                # Inbred/DH
                if parts[3].lower() in inbred:
                    ind.inbred = True
                elif parts[3].lower() in outbred:
                    ind.inbred = False


    def readInFromPlink(self, idList, pedList, bed, externalPedigree = False):
        index = 0

        if not externalPedigree:
            self.readInPedigreeFromList(pedList)
    
        for i, idx in enumerate(idList):
            genotypes=bed[:, i].copy() ##I think this is the right order. Doing the copy to be safe.
            nLoci = len(genotypes)
            if self.nLoci == 0:
                self.nLoci = nLoci
            if self.nLoci != nLoci:
                print(f"ERROR: incorrect number of loci when reading in plink file. Expected {self.nLoci} got {nLoci}.\nExiting...")
                sys.exit(2)
            if idx not in self.individuals:
                self.individuals[idx] = self.constructor(idx, self.maxIdn)
                self.maxIdn += 1

            ind = self.individuals[idx]
            ind.constructInfo(nLoci, genotypes=True)
            ind.genotypes = genotypes

            ind.fileIndex['plink'] = index; index += 1

            if np.mean(genotypes == 9) < .1 :
                ind.initHD = True


    
    def check_line(self, id_data, fileName, idxExpected=None, ncol=None, getInd=True, even_cols=False):
        idx, data = id_data

        if idxExpected is not None and idx != idxExpected:
            print(f"ERROR: Expected individual {idxExpected} but got individual {idx}.\nExiting...")
            sys.exit(2)
        if ncol is None:
            ncol = len(data)
        if ncol != len(data):
            print(f"ERROR: incorrect number of columns in {fileName}. Expected {ncol} values but got {len(data)} for individual {idx}.\nExiting...")
            sys.exit(2)
        if even_cols and ncol % 2 != 0:
            print(f"ERROR: file {fileName} doesn't contain an even number of allele columns for individual {idx}.\nExiting...")
            sys.exit(2)
        
        nLoci = len(data)
        if self.nLoci == 0:
            self.nLoci = nLoci
        if self.nLoci != nLoci:
            print(f"ERROR: inconsistent number of markers or alleles in {fileName}. Expected {self.nLoci} got {nLoci}.")
            sys.exit(2)

        ind = None
        if getInd :
            ind = self.getIndividual(idx)

        return ind, data, ncol 


    def update_allele_coding(self, alleles):
        """Update allele codings with new alleles
        self.allele_coding      - array of shape (2, n_loci) such that:
        self.allele_coding[0]   - array of alleles that are coded as 0
                                  (these are set to the first alleles 'seen')
        self.allele_coding[1]   - array of alleles that are coded as 1
                                  (these are alleles that are different from
                                  self.allele_coding[0])
        alleles  - alleles as read in from PLINK file
                array of dtype np.bytes_, b'0' is 'missing'

        This function is much like finding unique entries in a list:
        only add a new item if it is different from those seen before
        In this case, only record the first two uniques found,
        but also check there are only 2 alleles in total"""

        # If allele coding is complete then skip the update
        if self.allele_coding_complete():
            return

        # Update any missing entries in self.allele_coding[0]
        mask = self.allele_coding[0] == b'0'
        self.allele_coding[0, mask] = alleles[mask]

        # Update entries in self.allele_coding[1] if:
        # - the alleles have not already been seen in self.allele_coding[0]
        # - and the entry (in self.allele_coding[1]) is missing
        mask = self.allele_coding[1] == b'0'
        mask &= self.allele_coding[0] != alleles
        self.allele_coding[1, mask] = alleles[mask]

        # Check for > 2 alleles at any loci. These are:
        # - alleles that are not missing
        # - alleles that are not in either of self.allele_coding[0] or self.allele_coding[1]
        mask = alleles != b'0'
        mask &= ((alleles != self.allele_coding[0]) & (alleles != self.allele_coding[1]))  # poss speedup alleles != self.allele_coding[0] done above
        if np.sum(mask) > 0:
            print(f'ERROR: more than two alleles found in input file(s) at loci {np.flatnonzero(mask)}\nExiting...')
            sys.exit(2)


    def allele_coding_complete(self):
        """Check whether the allele coding is complete (contains no missing values)"""
        if self.allele_coding is None:
            return False
        else:
            return np.sum(self.allele_coding == b'0') == 0


    def decode_alleles(self, alleles):
        """Decode PLINK plain text alleles to AlphaGenes genotypes or haplotypes
        handles single individuals - alleles has shape (n_loci*2, )
        or multiple individuals    - alleles has shape (n_individuals, n_loci*2)"""
        # 'Double' self.allele_coding as there are two allele columns at each locus
        coding = np.repeat(self.allele_coding, 2, axis=1)

        decoded = np.full_like(alleles, b'0', dtype=np.int8)
        decoded[alleles == coding[0]] = 0  # set alleles coded as 0
        decoded[alleles == coding[1]] = 1  # set alleles coded as 1
        decoded[alleles == b'0'] = 9       # convert missing (b'0' -> 9)

        # Extract haplotypes
        decoded = np.atleast_2d(decoded)
        n_haps = decoded.shape[0] * 2
        n_loci = decoded.shape[1] // 2
        haplotypes = np.full((n_haps, n_loci), 9, dtype=np.int8)
        haplotypes[::2] = decoded[:, ::2]
        haplotypes[1::2] = decoded[:, 1::2]
        
        # Convert to genotypes
        genotypes = decoded[:, ::2] + decoded[:, 1::2]
        genotypes[genotypes > 9] = 9  # reset missing values

        return genotypes.squeeze(), haplotypes


    def encode_alleles(self, haplotypes):
        """Encode haplotypes as PLINK plain text
        handles any even number of haplotypes with shape (n_individuals*2, n_loci)"""
        assert len(haplotypes) % 2 == 0
        assert len(haplotypes[0])== self.nLoci

        # 'Double' self.allele_coding as there are two allele columns at each locus in PLINK format
        coding = np.repeat(self.allele_coding, 2, axis=1)

        # Encoded array is 'reshaped' - one individual per line, each locus is a pair of alleles
        encoded = np.full((len(haplotypes)//2, self.nLoci*2), b'0', dtype=np.bytes_)
        # 'Splice' haplotypes into (adjacent) pairs of alleles
        encoded[:, ::2] = haplotypes[::2]
        encoded[:, 1::2] = haplotypes[1::2]

        # Encode
        mask0 = encoded == b'0'  # major alleles (0)
        mask1 = encoded == b'1'  # minor alleles (1)
        mask9 = encoded == b'9'  # missing (9)
        encoded[mask0] = np.broadcast_to(coding[0], encoded.shape)[mask0]
        encoded[mask1] = np.broadcast_to(coding[1], encoded.shape)[mask1]
        encoded[mask9] = b'0'

        return encoded.squeeze()


    def check_allele_coding(self, filename):
        """Check coding is sensible"""
        # Monoallelic loci: 
        #   allele_coding[0] filled, but allele_coding[1] unfilled, i.e. coding[1] == b'0'
        n_monoallelic = (self.allele_coding[1] == b'0').sum()
        # Unusual letters
        unusual = ~np.isin(self.allele_coding, [b'A', b'C', b'G', b'T', b'0'])  # unexpected letters
        if np.sum(unusual) > 0:
            letters = ' '.join(np.unique(self.allele_coding[unusual].astype(str)))
            print(f'ERROR: unexpected values found in {filename}: [{letters}].\n'
                  f'Please check the file is in PLINK .ped format\nExiting...')
            sys.exit(2)
        elif n_monoallelic > 0:
            print(f'WARNING: allele coding from {filename} has {n_monoallelic} monoallelic loci')
        else:
            # Reassuring message if tests pass
            print('Allele coding OK')


    def readInPed(self, filename, startsnp=None, stopsnp=None, haps=False, update_coding=False):
        """Read in genotypes, and optionally haplotypes, from a PLINK plain text formated file, usually .ped
        If update_coding is True, the allele coding is interpreted from the .ped file and any coding
        in self.allele_coding is updated (if the coding is incomplete).
        Note: to force a full read of the allele coding set self.allele_codine = None first"""

        print(f'Reading in PLINK .ped format: {filename}')
        # Check the allele coding is to be got from file or is provided via self.allele_coding
        if not update_coding and self.allele_coding is None:
            raise ValueError('readInPed () called with no allele coding')

        data_list = MultiThreadIO.readLinesPlinkPlainTxt(filename, startsnp=startsnp, stopsnp=stopsnp, dtype=np.bytes_)
        index = 0
        ncol = None
        if self.nLoci != 0:
            # Temporarilly double nLoci while reading in PLINK plain text formats (two fields per locus)
            # otherwise reading of multiple PLINK files results in an 'Incorrect number of values'
            # error in check_line()
            self.nLoci = self.nLoci * 2

        if not self.allele_coding_complete():
            if self.allele_coding is None:
                print(f'Interpreting allele coding from {filename}')
            else:
                print(f'Updating allele coding with coding from {filename}')

        for value in data_list:
            ind, alleles, ncol = self.check_line(value, filename, idxExpected=None, ncol=ncol, even_cols=True)
            ind.constructInfo(self.nLoci, genotypes=True)
            ind.fileIndex['plink'] = index; index += 1

            if update_coding:
                # Initialise allele coding array if none exists
                # read_or_create = 'Reading' if self.allele_coding is None else 'Updating'
                if self.allele_coding is None:
                    self.allele_coding = np.full((2, self.nLoci//2), b'0', dtype=np.bytes_)

                # Update allele codes
                # print(f'{read_or_create} allele coding with coding from {filename}')
                self.update_allele_coding(alleles[::2])   # first allele in each pair
                self.update_allele_coding(alleles[1::2])  # second allele

            # Decode haplotypes and genotypes
            ind.genotypes, haplotypes = self.decode_alleles(alleles)
            if haps:
                ind.haplotypes = haplotypes

            if np.mean(ind.genotypes == 9) < .1:
                ind.initHD = True

        # Check allele coding
        self.check_allele_coding(filename)

        # Reset nLoci back to its undoubled state
        self.nLoci = self.nLoci//2


    def readInGenotypes(self, fileName, startsnp=None, stopsnp = None):

        print("Reading in AlphaGenes format:", fileName)
        index = 0
        ncol = None

        data_list = MultiThreadIO.readLines(fileName, startsnp = startsnp, stopsnp = stopsnp, dtype = np.int8)

        for value in data_list:
            ind, genotypes, ncol = self.check_line(value, fileName, idxExpected = None, ncol = ncol)

            ind.constructInfo(self.nLoci, genotypes=True)
            ind.genotypes = genotypes

            ind.fileIndex['genotypes'] = index; index += 1

            if np.mean(genotypes == 9) < .1 :
                ind.initHD = True

    def readInReferencePanel(self, fileName, startsnp=None, stopsnp = None):

        print("Reading in reference panel:", fileName)
        index = 0
        ncol = None

        data_list = MultiThreadIO.readLines(fileName, startsnp = startsnp, stopsnp = stopsnp, dtype = np.int8)

        for value in data_list:
            ind, haplotype, ncol = self.check_line(value, fileName, idxExpected = None, ncol = ncol, getInd=False)
            self.referencePanel.append(haplotype)

    def readInPhase(self, fileName, startsnp=None, stopsnp = None):
        index = 0
        ncol = None

        data_list = MultiThreadIO.readLines(fileName, startsnp = startsnp, stopsnp = stopsnp, dtype = np.int8)

        e = 0
        currentInd = None
        for value in data_list:

            if e == 0: 
                idxExpected = None
            else:
                idxExpected = currentInd.idx

            ind, haplotype, ncol = self.check_line(value, fileName, idxExpected = idxExpected, ncol = ncol)
            currentInd = ind

            ind.constructInfo(self.nLoci, haps=True)
            ind.haplotypes[e][:] = haplotype
            e = 1-e

            ind.fileIndex['phase'] = index; index += 1

        
    def readInSequence(self, fileName, startsnp=None, stopsnp = None):
        index = 0
        ncol = None

        print("Reading in sequence data :", fileName)
        
        data_list = MultiThreadIO.readLines(fileName, startsnp = startsnp, stopsnp = stopsnp, dtype = np.int64)
        e = 0
        currentInd = None

        for value in data_list:
            if e == 0: 
                idxExpected = None
            else:
                idxExpected = currentInd.idx

            ind, reads, ncol = self.check_line(value, fileName, idxExpected = idxExpected, ncol = ncol)
            currentInd = ind

            ind.constructInfo(self.nLoci, reads=True)
            ind.fileIndex['sequence'] = index; index += 1

            ind.reads[e][:] = reads
            e = 1-e


    def callGenotypes(self, threshold):
        for idx, ind in self.writeOrder():
            matrix = ProbMath.getGenotypeProbabilities_ind(ind, InputOutput.args)
            
            matrixCollapsedHets = np.array([matrix[0,:], matrix[1,:] + matrix[2,:], matrix[3,:]], dtype=np.float32)
            calledGenotypes = np.argmax(matrixCollapsedHets, axis = 0)
            setMissing(calledGenotypes, matrixCollapsedHets, threshold)
            if InputOutput.args.sexchrom and ind.sex == 0:
                doubleIfNotMissing(calledGenotypes)
            ind.genotypes = calledGenotypes


    def writePedigree(self, outputFile):
        with open(outputFile, 'w+') as f:
            for ind in self:
                sire = "0"
                if ind.sire is not None:
                    sire = ind.sire.idx
                dam = "0"
                if ind.dam is not None:
                    dam = ind.dam.idx
                f.write(ind.idx + ' ' + sire + ' ' + dam + '\n')


    def writeGenotypes(self, outputFile):

        data_list = []
        for ind in self :
            data_list.append( (ind.idx, ind.genotypes) )

        MultiThreadIO.writeLines(outputFile, data_list, str)


    def writePhase(self, outputFile):
        data_list = []
        for ind in self:
            if ind.haplotypes.ndim == 2:  # diploid
                data_list.append((ind.idx, ind.haplotypes[0]))
                data_list.append((ind.idx, ind.haplotypes[1]))
            elif ind.haplotypes.ndim == 1:  # haploid
                data_list.append((ind.idx, ind.haplotypes))
                data_list.append((ind.idx, ind.haplotypes))

        MultiThreadIO.writeLines(outputFile, data_list, str)


    def writeDosages(self, outputFile):
        data_list = []
        for ind in self :
            if ind.dosages is not None:
                data_list.append( (ind.idx, ind.dosages) )
            else:
                dosages = ind.genotypes.copy()
                dosages[dosages == 9] = 1
                data_list.append( (ind.idx, dosages) )

        MultiThreadIO.writeLines(outputFile, data_list, "{:.4f}".format)


    def writeGenotypes_prefil(self, outputFile):
        # print("Output is using filled genotypes. Filling missing with a value of 1")
        # fillValues = np.full(1, self.nLoci)

        print("Output is using filled genotypes. Filling missing with rounded allele frequency")
        self.setMaf()
        fillValues = np.round(self.maf)

        with open(outputFile, 'w+') as f:
            for idx, ind in self.individuals.items():
                fill(ind.genotypes, fillValues)
                self.writeLine(f, ind.idx, ind.genotypes, str)


    def writeGenotypesPed(self, outputFile):
        """Write genotypes in PLINK plain text format"""
        data_list = []
        for ind in self:
            # Split genotypes into 'pseudo' haplotypes such that
            # the first allele/haplotype of a heterozygous locus is always 0
            missing = ind.genotypes == 9
            h1 = ind.genotypes//2
            h2 = ind.genotypes - h1
            h1[missing] = 9
            h2[missing] = 9
            alleles = self.encode_alleles(np.vstack([h1, h2]))
            data_list.append( (ind.idx, alleles) )
        MultiThreadIO.writeLinesPlinkPlainTxt(outputFile, data_list)


    def writePhasePed(self, outputFile):
        """Write phased data (i.e. haplotypes) in PLINK plain text .ped format"""
        data_list = []
        for ind in self:
            if ind.haplotypes.ndim == 2:  # diploid
                alleles = self.encode_alleles(ind.haplotypes)
            elif ind.haplotypes.ndim == 1:  # haploid
                diploid = np.vstack([ind.haplotypes, ind.haplotypes])
                alleles = self.encode_alleles(diploid)
            data_list.append( (ind.idx, alleles) )

        MultiThreadIO.writeLinesPlinkPlainTxt(outputFile, data_list)


    def writeLine(self, f, idx, data, func) :
        f.write(idx + ' ' + ' '.join(map(func, data)) + '\n')


@jit(nopython=True)
def fill(genotypes, fillValue):
    for i in range(len(genotypes)):
        if genotypes[i] == 9:
            genotypes[i] = fillValue[i]

@jit(nopython=True)
def addIfNotMissing(array1, counts, array2):
    for i in range(len(array1)):
        if array2[i] != 9:
            array1[i] += array2[i]
            counts[i] += 2


@jit(nopython=True)
def addIfMissing(array1, array2):
    for i in range(len(array1)):
        if array2[i] == 9:
            array1[i] += 1

@jit(nopython=True)
def doubleIfNotMissing(calledGenotypes):
    nLoci = len(calledGenotypes)
    for i in range(nLoci):
        if calledGenotypes[i] == 1:
            calledGenotypes[i] = 2

@jit(nopython=True)
def setMissing(calledGenotypes, matrix, thresh) :
    nLoci = len(calledGenotypes)
    for i in range(nLoci):
        if matrix[calledGenotypes[i],i] < thresh:
            calledGenotypes[i] = 9















