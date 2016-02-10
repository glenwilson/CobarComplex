from dual_st_alg import Monomial, Polynomial, TensorMonomial, TensorPolynomial, all_indices, tau_deg, xi_deg, tau_wt, xi_wt
from options import opts
from cohomology import Cohomology
from mod_lin_alg import ModVector, ModMatrix
import random
import copy
import cPickle as pickle

class CobarMonomial(object):
    """
    initialized with a list of Monomial objects, and a coefficient
    """
    @staticmethod
    def random(x, filt):
        return CobarMonomial([Monomial.random(x) for i in
                              range(filt)], filt,
                             random.randrange(opts.prime))

    @staticmethod
    def unit(filt):
        return CobarMonomial([Monomial.unit() for i in range(filt)],
                             filt, 1)

    @staticmethod
    def null(filt):
        return CobarMonomial([Monomial.unit() for i in range(filt)],
                              filt, 0)

    def __init__(self, factors, filt, coeff):
        self.factors = factors
        self.filt = filt
        if len(self.factors) != self.filt:
            raise TypeError
        self.coeff = coeff

    def get_coeff(self):
        return self.coeff

    def get_factors(self):
        return self.factors
    
    def get_degree(self):
        deg = 0
        wt = 0
        for mon in self.factors:
            deg += mon.get_degree()[0]
            wt += mon.get_degree()[1]
        if self.coeff:
            return (deg, wt)
        else:
            return (0,0)

    def __str__(self):
        if not self.coeff:
            return "0"
        string = str(self.coeff) + "*"
        for mon in self.factors:
            string += str(mon) + " | "
        return string[:-3]

    def simplify(self):
        """Warning! Mutatior!

        """
        for mon in self.factors:
            mon.simplify()
        newcoeff = self.coeff
        for mon in self.factors:
            newcoeff = (newcoeff * mon.coeff) % opts.prime
            mon.coeff = 1
        self.coeff = newcoeff
        if not self.coeff:
            self.factors = CobarMonomial.null(self.filt).factors

    def __eq__(self, other):
        """
        mutator
        """
        self.simplify()
        other.simplify()
        if self.filt != other.filt or self.coeff != other.coeff:
            return False
        else:
            for i in xrange(self.filt):
                if not self.factors[i] == other.factors[i]:
                    return False
        return True

    def factors_equal(self, other):
        """warning, mutator! will check if two CObarMOnomials have the same
        factors, ignoring the coefficient unless one is 0

        """
        self.simplify()
        other.simplify()
        if self.filt != other.filt:
            return False
        else:
            for i in xrange(self.filt):
                if not self.factors[i] == other.factors[i]:
                    return False
        return True
        
    # def __mul__(self, other):
    #     ### THIS IS THE WRONG MUL! THIS IS CONCATENATION
    #     if isinstance(other, CobarMonomial):
    #         newcoeff = (self.coeff * other.coeff) %opts.prime
    #         newfactors = self.factors + other.factors
    #         return CobarMonomial(newfactors, len(newfactors), newcoeff)
    #     elif isinstance(other, TensorMonomial):
    #         newcoeff = (self.coeff * other.coeff) % opts.rime
    #         newfactors = self.factors + other.tuple
    #         return CobarMonomial(newfactors, len(newfactors), newcoeff)
    #     elif isinstance(other, Monomial):
    #         newfactors = self.factors[:]
    #         newfactors.append(other)
    #         return CobarMonomial(newfactors, len(newfactors), self.coeff)
    #     else:
    #         raise TypeError

    def __rmul__(self, n):
        newcoeff = (n * self.coeff) % opts.prime
        newfactors = self.factors[:]
        xx = CobarMonomial(newfactors, len(newfactors), newcoeff)
        xx.simplify()
        return xx
        
    def __and__(self, other):
        """
        this defines the concatenation product 

        this is invoked by self & other
        """
        if isinstance(other, CobarMonomial):
            out = self.factors + other.factors
            newcoeff = (self.coeff * other.coeff) % opts.prime
            return CobarMonomial(out, len(out), newcoeff)
        elif isinstance(other, CobarPolynomial) \
             or isinstance(other, TensorPolynomial):
            newsummands = []
            for mon2 in other.get_summands():
                newsummands.append(self & mon2)
            return CobarPolynomial(newsummands)
        elif isinstance(other, TensorMonomial):
            out = self.factors + other.pair
            newcoeff = (self.coeff * other.coeff) % opts.prime
            return CobarMonomial(out, len(out), newcoeff)
        else:
            print "self", self, type(self), self.factors
            print "other", other, type(other)
            raise TypeError

    def insert(self, i, monomial):
        """Returns a new CobarMonomial with monomial inserted into the ith
        spot in the list of factors

        """
        newfactors = self.factors[:]
        newfactors.insert(i, monomial)
        return CobarMonomial(newfactors, len(newfactors), self.coeff)

    def _map(self):
        out = CobarPolynomial([])
        self.simplify()
        print self
        ###
        out = out + (CobarMonomial.unit(1) & self)
        print out
        out = out + ((-1)**(self.filt+1) % opts.prime) * (self & CobarMonomial.unit(1))
        print out
        ###
        
        for i in range(self.filt):
            mon = self.factors[i]
            #change this to .reduced_coproduct() if you figure out
            #correct deal for reduced thing!
            coprod = mon.coproduct()
#            print "mon", mon
#            print "coprod", coprod
            left = CobarMonomial(self.factors[:i], i, 1)
            right = CobarMonomial(self.factors[i+1:], \
                                  len(self.factors[i+1:]), 1)
            out = out + self.coeff * (((-1)**(i+1) % opts.prime) \
                         * ((left & coprod) & right))
            print out
        out.simplify()
        return out

    def _map_reduced(self):
        out = CobarPolynomial([])
        self.simplify()
        for i in range(self.filt):
            mon = self.factors[i]
            coprod = mon.reduced_coproduct()
            left = CobarMonomial(self.factors[:i], i, 1)
            right = CobarMonomial(self.factors[i+1:], \
                                  len(self.factors[i+1:]), 1)
            out = out + (self.coeff * (((-1)**(i+1) % opts.prime )) * ((left & coprod) & right))
            #for summand in out.summands:
            #    coeff_fix = 1
            #    for term in summand.factors[:i+1]:
            #        coeff_fix *= (-1)**(mon.get_degree()[0]+1)
            #    coeff_fix = coeff_fix % opts.prime
            #    summand.coeff = summand.coeff * coeff_fix    
            #print out
        out.simplify()
        return out

        
            
class CobarPolynomial(object):
    """
    initialized with a list of CobarMonomial objects.
    """
    @staticmethod
    def random(x, y, filt):
        summands = [ CobarMonomial.random(x, filt) for i in range(y)]
        return CobarPolynomial(summands)

    @staticmethod
    def unit(filt):
        return CobarPolynomial([CobarMonomial.unit(filt)])

    @staticmethod
    def null():
        return CobarPolynomial([])

    def __init__(self, summands):
        self.summands = summands
        if not summands:
            self.filt = 0
        else:
            self.filt = summands[0].filt
            for term in summands:
                if term.filt != self.filt:
                    print "not homogeneous!"

    def get_summands(self):
        return self.summands

    def __str__(self):
        string = ""
        if not self.summands:
            return "0"
        for mon in self.get_summands():
            string += str(mon) + " + "
        return string[:-3]

    def __add__(self,other):
        if isinstance(other, CobarPolynomial):
            out = self.summands[:] + other.summands
        elif isinstance(other, CobarMonomial):
            out = self.summands[:]
            out.append(other)
        else:
            raise TypeError
        return CobarPolynomial(out)

    def stupid_simplify(self):
        for mon in self.summands[:]:
            if mon == CobarMonomial.null(mon.filt):
                self.summands.remove(mon)

    def simplify(self):
        self.stupid_simplify()
        stack = self.get_summands()[:]
        #this simplify may be unnecessary
        for mon in stack:
            mon.simplify()
        outsum = []
        while stack:
            mon = stack.pop()
            similar = [ x for x in stack if mon.factors_equal(x) ]
            coefficient = mon.get_coeff()
            for x in similar:
                coefficient = (coefficient + x.get_coeff()) \
                              % opts.prime
                stack.remove(x)
            if coefficient:
                outsum.append(CobarMonomial(mon.get_factors(),
                                            mon.filt, coefficient))
        self.summands = outsum

    def get_degree(self):
        self.simplify()
        if not self.summands:
            return (0,0)
        deg = self.summands[0].get_degree()[0]
        wt = self.summands[0].get_degree()[1]
        for mon in self.summands:
            if (deg, wt) != mon.get_degree():
                print (deg, wt)
                print mon.get_degree()
                print self
                return False
        return (deg, wt)
            
        
    # def __mul__(self, other):
    #     # THIS IS THE WRONG MUL!!!
    #     newsummands = []
    #     for mon1 in self.get_summands():
    #         for mon2 in other.get_summands():
    #             newsummands.append(mon1*mon2)
    #     return CobarPolynomial(newsummands)

    def copy(self):
        return copy.deepcopy(self)
    
    def __rmul__(self, n):
        """ scalar multiplication by n"""
        out = self.copy()
        for mon in out.summands:
            mon.coeff = (n * mon.coeff) % opts.prime
        return out
            
    def __and__(self, other):
        """
        This defines the concatenation product 

        invoked with self & other
        """
        if isinstance(other, CobarPolynomial) \
           or isinstance(other, TensorPolynomial) :
            newsummands = []
            for mon1 in self.get_summands():
                for mon2 in other.get_summands():
                    newsummands.append(mon1 & mon2)
            return CobarPolynomial(newsummands)
        elif isinstance(other, CobarMonomial) \
             or isinstance(other, TensorMonomial):
            newsummands = []
            for mon1 in self.get_summands():
                newsummands.append(mon1 & other)
            return CobarPolynomial(newsummands)

    def is_summand(self, monomial):
        """checks if the monomial appears in the polynomial with non-zero
        coefficient

        """
        self.simplify()
        for mon in self.summands:
            if mon.factors_equal(monomial):
                return True
        return False

    def monomial_coefficient(self, monomial):
        self.simplify()
        for mon in self.summands:
            if mon.factors_equal(monomial):
                return mon.coeff
        return 0

        
class CobarModule(object):

    def __init__(self):
        # keys are bidegrees (deg, wt)
        # value is a list consisting of
        # CobarMonomials which form a basis
        self._dict = {}
        # this retains the outgoing maps as matrices
        self._maps = {}
        #check flag to see if dictionary has been
        #populated with bases
        self.flag = False

#    def elt_to_vector(self, elt):

#    def vector_to_elt(self, vector, bidegree):
    
        
        
class CobarComplex(object):
    """
    This is a wrapper class
    """
    def __init__(self, length = opts.bounds):
        self.length = length
        self.cplx = [CobarModule() for i in range(self.length + 1)]
        self.module_flag = False
        self.map_flag = False

    def generate_modules(self):
        """
        This will produce the first opts.bounds +1 modules 
        in the cobar complex

        check to make sure only simplified monomials appear in basis!
        """
        if self.module_flag:
          return  
        #this sets up 0th module
        self.cplx[0]._dict[(0,0)] = [CobarMonomial([], 0, 1)]
        for f in range(self.length):
             self.cplx[f+1]._dict[(0,0)] = []
        #this will add all generators fo first module
        first = self.cplx[1]._dict
        Indices = all_indices()
        #print Indices
        for pair in Indices:
            elt = CobarMonomial([Monomial.tau_xi_list(pair)], 1, 1)
            deg = elt.get_degree()
            #print elt
            #print deg
            try:
                first[deg] += [elt]
            except KeyError:
                first[deg] = [elt]
        # Now proceed to fill in remaining modules
        for f in range(2,self.length+1):
            prev = self.cplx[f-1]._dict
            curr = self.cplx[f]._dict
            # for each existing bidegree
            for bideg in prev.keys():
                #try to insert a monomial 
                for pair in Indices:
                    deg = bideg[0] + tau_deg(pair[0]) + xi_deg(pair[1])
                    wt = bideg[1] + tau_wt(pair[0]) + xi_wt(pair[1])
                    # check degree works
                    if deg <= opts.bounds:
                        try:
                            curr[(deg, wt)]
                        except KeyError:
                            curr[(deg, wt)] = []
                        #if it does, insert
                        for monomial in prev[bideg]:
                            for k in range(f):
                                elt = monomial.insert(k, Monomial.tau_xi_list(pair))
                                flag = False
                                for other in curr[(deg, wt)]:
                                    if other == elt:
                                        flag = True
                                        break
                                if not flag:
                                    curr[(deg, wt)] += [elt]
                    else:
                        break
        self.module_flag = True

    def extend_bounds(self, bounds):
        """This will take an existing complex and extend the modules which are
        present to be complete up to the given bounds

        """
        pass
        
    def extend_modules(self, length):
        """This will pickup an existing complex of the given length and add
        on the additional modules to the complex.

        """
        pass

    def product_structure(self):
        pass

    def massey_products(self):
        pass
    
    def get_cplx(self):
        """Get this to implement unpickling. check to see if there is a
        pickled file? Or will that be too time consuming?

        """
        if not self.module_flag:
            self.generate_modules()
        return self.cplx

    def get_maps(self):
        if not self.map_flag:
            self.make_maps()
        return self.cplx
    
    def make_map(self, i):
        """
        make matrices in all given bidegrees

        this needs a map defined on CobarMonomial objects

        does this handle 0 dim vector spaces? 
        """
        #print "in make map", i
        dom = self.get_cplx()[i]
        rng = self.get_cplx()[i+1]
        for bideg in dom._dict.keys():
            out = []
            #print "at bideg", bideg
            dom_basis = dom._dict[bideg]
            cols = len(dom_basis)
            try:
                rng_basis = rng._dict[bideg]
                rows = len(rng_basis)
            except KeyError:
                rng._dict[bideg] = []
                rows = 0
                if cols == 0:
                    # we must account for this later!
                    cols = 1
                #print "bideg in make maps",bideg
                dom._maps[bideg] = ModMatrix.null(1,cols)
                continue
            if cols == 0:
                cols = 1
                if rows == 0:
                    rows = 1
                #print "bideg in make maps. somethings 0", bideg
                dom._maps[bideg] = ModMatrix.null(rows, cols)
                continue
            for mon in dom_basis:
                #print "filtration", i, "at bidegree", bideg
                #print "mon", mon
                value = mon._map_reduced()
                #print "val", value
                vect = self.vector_from_element(value, i+1, bideg[0],
                                                bideg[1])
                #print "vector to append", vect
                out.append(vect)
            matrix = ModMatrix(out).get_transpose()
            #print matrix
            dom._maps[bideg] = matrix 
            
    def make_maps(self):
        if self.map_flag:
            return
        for i in range(self.length):
            self.make_map(i)
        self.map_flag = True

    def get_cohomology(self, filt, deg, wt):
        """
        returns a cohomology object corr. to filt, deg, wt
        #Returns 0 cohomology object if deg, wt doesn't appera

        THIS MUTATES ._maps and ._dict IF NECESSARY!
        """
        cplx = self.get_maps()
        if filt == 0:
            return "you know what this is!"
        try:
            #this checks to make sure cohom is non-zero here
            B = cplx[filt]._maps[(deg, wt)]
            B_basis = cplx[filt]._dict[(deg, wt)]
        except KeyError:
            if deg <= opts.bounds:
                return Cohomology(ModMatrix.identity(1), \
                                  ModMatrix.null(1,1))
            else:
                print filt, deg, wt
                raise TypeError("We aren't computing that far!")
        try:
            #this checks to see if incoming map is 0 or not
            A_basis = cplx[filt - 1]._dict[(deg, wt)]
            A = cplx[filt - 1]._maps[(deg, wt)]
        except KeyError:
            A = ModMatrix.null(len(B_basis), 1)
            A_basis = []
            cplx[filt - 1]._maps[(deg, wt)] = A
            cplx[filt - 1]._dict[(deg, wt)] = []
            # We must now check dimensions of A_basis and B_basis
            # to account for 0 diml vspaces
#        print "Abasis", len(A_basis)
#        print "Bbasis", len(B_basis)
        if len(B_basis) == 0:
            #this means cohom is 0 diml
            #resulting cohom object will not match up with vect/elt
            #translation
            #print "B basis null"
            A = ModMatrix.null(1,1)
            B = ModMatrix.identity(1)
        elif not B.row_count:
            #print "B has no rows"
            B = ModMatrix.null(1, len(B_basis))
        cohom = Cohomology(B, A)
#        print "bideg key error", filt, deg, wt
#        print "filt - 1",  cplx[filt - 1]._dict.keys()
#        print "filt", cplx[filt]._dict.keys()
#        print "filt - 1 ", cplx[filt - 1]._maps.keys()
#        print "filt", cplx[filt]._maps.keys()
        return cohom

    def vector_from_element(self, elt, filt, deg, wt):
        """returns vector in standard basis as determined by the
        generate_modules algorithm corr. to elt.

        may need to add in a case to handle 0 elements
        """
        elt.simplify()
        try:
            basis = self.get_cplx()[filt]._dict[(deg, wt)]
        except KeyError:
            print elt
            print filt, (deg, wt)
            raise KeyError("Elt does not appear in complex!")
        vect = ModVector.null(len(basis))
        for i in range(len(basis)):
            mon = basis[i]
            vect[i] = elt.monomial_coefficient(mon)
        return vect
            

    def element_from_vector(self, vector, filt, deg, wt):
        """
        returns a CobarPolynomial corresponding to vector from the 
        standard basis in filt, deg, wt.
        """
        elt = CobarPolynomial([])
        basis = self.get_cplx()[filt]._dict[(deg,wt)]
        if len(vector.vector) != len(basis):
            print "vector", vector
            print "basis"
            for thing in basis:
                print thing
            raise TypeError("This vector doesn't seem to belong here")
        for i in range(len(basis)):
            elt = elt + ( vector[i] * basis[i])
        elt.simplify()
        return elt

    def pickle_cplx(self):
        filename = "complex" + str(opts.prime) + \
        "-" + str(opts.bounds) + "-" + str(self.length) + ".pickle"
        str(self.length) +"-" + str(opts.bounds) + ".pickle"
        pickle_file = open(filename, "wb")
        pickle.dump(self.cplx, pickle_file, -1)
        pickle_file.close()

    def get_pickled_cplx(self):
        filename = "complex" + str(opts.prime) + \
        "-" + str(opts.bounds) + "-" + str(self.length) + ".pickle"
        try:
            pickle_file = open(filename, "rb")
            self.cplx = pickle.load(pickle_file)
            pickle_file.close()
            self.module_flag = True
            self.map_flag = True
        except IOError:
            print "No pickled file"
        
        
