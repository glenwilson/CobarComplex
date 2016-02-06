
class Options(object):
    def __init__(self, prime, degree_bounds, numpcs):
        self.prime = prime
        self.bounds = degree_bounds
        self.numpcs = numpcs

    def get_prime(self):
        return self.prime

    def get_bounds(self):
        return self.bounds

    def get_numpcs(self):
        return self.numpcs

opts = Options(17, 10, 2)

import random 
import copy
class Monomial(object):
    """terms is a list of tuples of the form (y, i) where y is a string
    "t" or "x" and i is a natural number. We assume terms of the form
    ("x",0) do not appear. Such terms are 1 in the dual steenrod
    algebra. An empty list is treated as (coeff * 1). 

    """
    @staticmethod
    def random(length):
        terms = []
        for i in range(length):
            randint = random.randrange(6)
            parity = random.randrange(2)
            if parity:
                terms.append(("t", randint))
            else:
                terms.append(("x", randint+1))

        return Monomial( terms,  random.randrange(1,opts.get_prime()))
        
    def __init__(self, terms, coeff):
        self.terms = terms
        self.coeff = coeff

    def __getstate__(self):
        return {'terms' : self.terms,
                'coeff' : self.coeff}
        
    def __setstate__(self, _dict):
        self.terms = _dict['terms']
        self.coeff = _dict['coeff']

    def get_terms(self):
        return self.terms

    def get_coeff(self):
        return self.coeff % opts.get_prime()
    
    def get_degree(self):
        """
        Possibly time consuming, simplifies! 
        """
        self.simplify()
        deg = 0
        wt = 0
        for term in self.terms:
            if term[0] == "x":
                deg += 2*(opts.prime**term[1] - 1) 
                wt += opts.prime**term[1] - 1
            elif term[0] == "t":
                deg += 2*opts.prime**term[1] - 1
                wt += opts.prime**term[1] - 1
        return (deg, wt)

    def copy(self):
        return copy.deepcopy(self)
                
    def __str__(self):
        mystr = ""
        mystr += str(self.get_coeff()) 
        if not self.terms:
            return mystr
        mystr += "."
        for term in self.terms:
            mystr += term[0] + str(term[1])
        return mystr
    
    def simplify(self):
        """
        Warning! Mutator!
        """
        newterms = sorted([ term for term in self.terms if term[0] == "x" and term[1] > 0], key=lambda term: term[1])
        # tauterms = [ term for term in self.terms if term[0] == "t"]
        tauindices = [ term[1] for term in self.terms if term[0] == "t"]
        newcoeff = self.coeff
        if repeats(tauindices):
            self.terms = []
            self.coeff = 0
        elif not self.get_coeff():
            self.terms = []
            self.coeff = 0
        else:
            tausorted = tauindices[:]
            tausorted.sort()
            perm = arePermsEqualParity(tauindices, tausorted)
            if not perm:
                self.coeff = -self.coeff
            newterms += [("t", x) for x in tausorted]
            self.terms = newterms

    def __mul__(self,other):
        """
        Does not mutate self or other. 
        Returns the product
        """
        prodterms = self.get_terms() + other.get_terms()
        prodcoeff = self.get_coeff() * other.get_coeff()
        prod = Monomial(prodterms, prodcoeff)
        prod.simplify()
        return prod

    def __eq__(self, other):
        """
        This should maybe mutate monomials to simplified form?
        Currently uses costly copies
        """
        selfcopy = self.copy()
        othercopy = other.copy()
        selfcopy.simplify()
        othercopy.simplify()
        if selfcopy.coeff != othercopy.coeff:
            return False
        if selfcopy.terms != othercopy.terms:
            return False
        return True
        
class Polynomial(object):
    """
    summands is a list of monomial objects
    """
    @staticmethod
    def random(x, y):
        """
        x is the number of monomials
        y is the length of each monomial
        """
        summands = [ Monomial.random(y) for i in range(x) ]
        return Polynomial(summands)

    def __init__(self, summands):
        self.summands = summands

    def __str__(self):
        string = ""
        if not self.get_summands():
            return "0"
        for monomial in self.get_summands():
            string += str(monomial) + " + "
        return string[:-3]
        
    def get_summands(self):
        return self.summands

    def homogeneous(self):
        """
        Does NOT simplify!!!
        """
        degs = { mon.get_degree() for mon in self.get_summands() }
        if len(degs) > 1:
            return False
        else:
            return True

    def get_degree(self):
        if self.homogeneous():
            return self.get_summands()[0].get_degree()
        else:
            return False

    def __add__(self, other):
        return Polynomial(self.get_summands() + other.get_summands())

    def __mul__(self, other):
        newsummands = []
        for mon1 in self.get_summands():
            for mon2 in other.get_summands():
                newsummands.append(mon1*mon2)
        return Polynomial(newsummands)

    def stupid_simplify(self):
        for monomial in self.get_summands()[:]:
            if monomial == Monomial([],0):
                self.summands.remove(monomial)            
    
    def simplify(self):
        """
        Mutator!
        """
        self.stupid_simplify()
        stack = self.get_summands()[:]
        # might need case to rule out 0
        for mon in stack:
            mon.simplify()
        outsum = []
        while stack:
            mon = stack.pop()
            similar = [ x for x in stack if str(x)[2:] == str(mon)[2:]]
            coefficient = mon.get_coeff()
            for x in similar:
                coefficient += x.get_coeff() % opts.prime
                stack.remove(x)
            if coefficient:
                outsum.append(Monomial(mon.get_terms(), coefficient))
        self.summands = outsum

    def copy(self):
        return copy.deepcopy(self)
        
    def __eq__(self, other):
        """
        This should maybe mutate monomials to simplified form?
        Currently uses costly copies
        """
        selfcopy = self.copy()
        othercopy = other.copy()
        selfcopy.simplify()
        othercopy.simplify()
        stack1 = selfcopy.get_summands()[:]
        stack2 = othercopy.get_summands()[:]
        for x in selfcopy.get_summands():
            if x in stack2:
                stack1.remove(x)
                stack2.remove(x)
            else:
                return False
        if len(stack1) or len(stack2):
            return False
        else:
            return True
                

class TensorMonomial(object):
    """
    pair is a list [x,y] with x and y Monomials, 
    corresponds to x \otimes y
    """
    def __init__(self, pair, coeff):
        self.pair  = pair
        self.coeff = coeff

    def get_pair(self):
        return self.get_pair

    def get_coeff(self):
        return self.coeff
    
    def get_degree(self):
        adeg = pair[0].get_degree()
        bdeg = pair[1].get_degree()
        if coeff:
            return (adeg[0] + bdeg[0], adeg[1] + bdeg[1])
        else:
            return (0,0)

    def __str__(self):
        if not self.coeff:
            return "0"
        return str(self.coeff) + "*" + str(self.pair[0]) + " | " + str(self.pair[1])

    def simplify(self):
        """
        Warning! Mutator!
        """
        self.pair[0].simplify()
        self.pair[1].simplify()
        newcoeff = self.pair[0].coeff * self.pair[1].coeff * self.coeff % opts.prime
        self.pair[0].coeff = 1
        self.pair[1].coeff = 1
        self.coeff = newcoeff
        if not self.coeff:
            self.pair[0] = Monomial([],0)
            self.pair[1] = Monomial([],0)

    def __mul__(self, other):
        newpair = [self.pair[0] * other.pair[0], self.pair[1] * other.pair[1] ]
        newcoeff = self.coeff * other.coeff % opts.prime
        return TensorMonomial(newpair, newcoeff)

    def __eq__(self, other):
        """
        Warning! Mutator...
        """
        self.simplify()
        other.simplify()
        if self.coeff != other.coeff:
            return False
        elif self.pair[0] != other.pair[0]:
            return False
        elif self.pair[1] != other.pair[1]:
            return False
        return True

    def monomials_equal(self, other):
        """Warning! mutator!  This will check if the monomials are equal,
        ignoring the coefficients, unless they are 0.

        """
        self.simplify()
        other.simplify()
        if self.pair[0] != other.pair[0]:
            return False
        elif self.pair[1] != other.pair[1]:
            return False
        return True

    
class TensorPolynomial(object):
    """
    summands is a list of TensorMonomial objects
    """
    def __init__(self, summands):
        self.summands = summands

    def get_summands(self):
        return self.summands

    def stupid_simplify(self):
        for mon in self.summands[:]:
            if mon == TensorMonomial([Monomial([],0),
                                      Monomial([],0)],0):
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
            similar = [ x for x in stack if mon.monomials_equal(x) ]
            coefficient = mon.get_coeff()
        for x in similar:
            coefficient += x.get_coeff() % opts.prime
            stack.remove(x)
        if coefficient:
            outsum.append(TensorMonomial(mon.get_pair(), coefficient))
        self.summands = outsum

    def __str__(self):
        string = ""
        if not self.summands:
            return "0"
        for mon in self.get_summands():
            string += str(mon) + " + "
        return string[:-3]

    def __add__(self, other):
        return TensorPolynomial(self.get_summands() \
                                + other.get_summands())

    def __mul__(self, other):
        newsummands = []
        for mon1 in self.get_summands():
            for mon2 in other.get_summands():
                newsummands.append(mon1*mon2)
        return TensorPolynomial(newsummands)



    #multiply
    
def repeats(alist):
    for i in alist:
        if alist.count(i) > 1:
            return True
    return False
            
def arePermsEqualParity(perm0, perm1):
    """Check if 2 permutations are of equal parity.

    Assume that both permutation lists are of equal length and have
    the same elements. No need to check for these conditions.

    Thanks to Stackoverflow "how-to-check-if-permutations-have-equal-parity"

    """
    perm1 = perm1[:] ## copy this list so we don't mutate the
                             ## original

    transCount = 0
    for loc in range(len(perm0) - 1): # Do (len - 1) transpositions
        p0 = perm0[loc]
        p1 = perm1[loc]
        if p0 != p1:
            sloc = perm1[loc:].index(p0)+loc          # Find position in perm1
            perm1[loc], perm1[sloc] = p0, p1          # Swap in perm1
            transCount += 1
                            
                            # Even number of transpositions means equal parity
    if (transCount % 2) == 0:
        return True
    else:
        return False
