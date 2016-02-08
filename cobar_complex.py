from dual_st_alg import Monomial, Polynomial, TensorMonomial, TensorPolynomial, all_indices, tau_deg, xi_deg, tau_wt, xi_wt
from options import opts
from cohomology import Cohomology
from mod_lin_alg import ModVector, ModMatrix
import random

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
        
    def __mul__(self, other):
        if isinstance(other, CobarMonomial):
            newcoeff = (self.coeff * other.coeff) %opts.prime
            newfactors = self.factors + other.factors
            return CobarMonomial(newfactors, len(newfactors), newcoeff)
        elif isinstance(other, TensorMonomial):
            newcoeff = (self.coeff * other.coeff) % opts.rime
            newfactors = self.factors + other.tuple
            return CobarMonomial(newfactors, len(newfactors), newcoeff)
        elif isinstance(other, Monomial):
            newfactors = self.factors[:]
            newfactors.append(other)
            return CobarMonomial(newfactors, len(newfactors), self.coeff)
        else:
            raise TypeError
        
    def __and__(self, other):
        """
        this defines the concatenation product 

        this is invoked by self & other
        """
        if isinstance(other, CobarMonomial):
            out = self.factors + other.factors
            newcoeff = (self.coeff * other.coeff) % opts.prime
            return CobarMonomial(out, len(out), newcoeff)
        elif isinstance(other, CobarPolynomial):
            newsummands = []
            for mon2 in other.get_summands():
                newsummands.append(self & mon2)
            return CobarPolynomial(newsummands)
        elif isinstance(other, TensorMonomial):
            out = self.factors + other.pair
            newcoeff = (self.coeff * other.coeff) % opts.prime
            return CobarMonomial(out, len(out), newcoeff)
        else:
            raise TypeError

    def insert(self, i, monomial):
        """Returns a new CobarMonomial with monomial inserted into the ith
        spot in the list of factors

        """
        newfactors = self.factors[:]
        newfactors.insert(i, monomial)
        return CobarMonomial(newfactors, len(newfactors), self.coeff)
        
            
class CobarPolynomial(object):
    """
    initialized with a list of CobarMonomial objects.
    """
    @staticmethod
    def random(x, y, filt):
        summands = [ CobarMonomial.random(y, filt) for i in range(x)]
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
        out = self.summands + other.summands
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

    def __mul__(self, other):
        newsummands = []
        for mon1 in self.get_summands():
            for mon2 in other.get_summands():
                newsummands.append(mon1*mon2)
        return CobarPolynomial(newsummands)

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

class CobarModule(object):

    def __init__(self):
        # keys are bidegrees (deg, wt)
        # value is a list consisting of
        # CobarMonomials which form a basis
        self._dict = {}
        #check flag to see if dictionary has been
        #populated with bases
        self.flag = False

    #def elt_to_vector(self, elt):

    #def vector_to_elt(self, vector, bidegree):
    
        
        
class CobarComplex(object):
    """
    This is a wrapper class
    """
    def __init__(self):
        self.length = opts.bounds
        self.cplx = [CobarModule() for i in range(self.length + 1)]


    def generate_modules(self):
        """
        This will produce the first opts.bounds +1 modules 
        in the cobar complex
        
        """
        #this sets up 0th module
        self.cplx[0]._dict[(0,0)] = [CobarMonomial([], 0, 1)]
        #this will add all generators fo first module
        first = self.cplx[1]._dict
        Indices = all_indices()
        for pair in Indices:
            elt = CobarMonomial([Monomial.tau_xi_list(pair)], 1, 1)
            deg = elt.get_degree()
            try:
                first[deg] += [elt]
            except KeyError:
                first[deg] = [elt]
        # Now proceed to fill in remaining modules
        for f in range(2,opts.bounds+1):
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
        
                                


        
            
        
