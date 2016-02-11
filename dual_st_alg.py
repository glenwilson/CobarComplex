from options import opts
from math import log
import random 
import copy

class Monomial(object):
    """terms is a list of tuples of the form (y, i) where y is a string
    "t" or "x" and i is a natural number. We assume terms of the form
    ("x",0) do not appear. Such terms are 1 in the dual steenrod
    algebra. An empty list is treated as (coeff * 1). 

    Should implement a flag system for checking simplification. 
    This will cut down on repeated simplifications.
    """
    @staticmethod
    def random(length):
        terms = []
        for i in range(length):
            randint = random.randrange(4)
            parity = random.randrange(2)
            if parity:
                terms.append(("t", randint))
            else:
                terms.append(("x", randint+1))

        return Monomial( terms,  random.randrange(opts.get_prime()))

    @staticmethod
    def unit():
        return Monomial([],1)
    @staticmethod
    def null():
        return Monomial([],0)
    
    @staticmethod
    def term_coproduct(term):
        """
        term must be a tuple of the form ("x", i) or ("t", i)
        """
        if term[0] == "x":
            if term[1] == 0:
                return TensorPolynomial.unit()
            out = TensorPolynomial.null()
            out += TensorPolynomial([TensorMonomial([Monomial.unit(), Monomial([term], 1)], 1)])
            out += TensorPolynomial([TensorMonomial([Monomial([term], 1), Monomial.unit()], 1)])
            for j in range(1, term[1]):
                out += TensorPolynomial([TensorMonomial([Monomial([("x", term[1]-j)], 1)**(opts.prime**j),\
                                                         Monomial([("x", j)],1)],1)])
            return out
        elif term[0] == "t":
            out = TensorPolynomial.null()
            out += TensorPolynomial([TensorMonomial([Monomial.unit(), Monomial([term], 1)],1)])
            out += TensorPolynomial([TensorMonomial([Monomial([term], 1), Monomial.unit()],1)])
            for j in range(0, term[1]):
                out += TensorPolynomial([TensorMonomial([Monomial([("x", term[1]-j)],1)**(opts.prime**j),\
                                                         Monomial([("t", j)],1)],1)])
            return out

    @staticmethod
    def tau_list(indices):
        out = []
        for i in range(len(indices)):
            for j in range(indices[i]):
                out.append(("t",i))
        return Monomial(out, 1)

    @staticmethod
    def xi_list(indices):
        out = []
        for i in range(len(indices)):
            for j in range(indices[i]):
                out.append(("x", i+1))
        return Monomial(out, 1)

    @staticmethod
    def tau_xi_list(pair):
        xx = Monomial.tau_list(pair[0])
        yy = Monomial.xi_list(pair[1])
        return xx * yy
        
    def __init__(self, terms, coeff):
        self.terms = terms
        self.coeff = coeff
        self.simplify_flag = False

    def __getstate__(self):
        return {'terms' : self.terms,
                'coeff' : self.coeff,
                'simplify_flag' : self.simplify_flag}
        
    def __setstate__(self, _dict):
        self.terms = _dict['terms']
        self.coeff = _dict['coeff']
        self.simplify_flag = _dict['simplify_flag']

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
        try:
            if self.simplify_flag:
                return
        except AttributeError:
            print self
        newterms = sorted([ term for term in self.terms if term[0] == "x" and term[1] > 0], key=lambda term: term[1])
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
        self.simplify_flag = True

    def __mul__(self,other):
        """
        Does not mutate self or other. 
        Returns the product
        """
        prodterms = self.get_terms() + other.get_terms()
        prodcoeff = (self.get_coeff() * other.get_coeff() % opts.prime)
        prod = Monomial(prodterms, prodcoeff)
        prod.simplify()
        return prod

    def __pow__(self, num):
        out = Monomial([],1)
        for i in xrange(num):
            out = out * self
        out.simplify()
        return out
    
    def __eq__(self, other):
        """
        This should maybe mutate monomials to simplified form?
        Currently uses costly copies
        """
        self.simplify()
        other.simplify()
        if self.coeff != other.coeff:
            return False
        if self.terms != other.terms:
            return False
        return True

    def coproduct(self):
        self.simplify()
        out = TensorPolynomial([TensorMonomial([Monomial.unit(),
                                                Monomial.unit()],
                                               self.get_coeff())])
        for term in self.get_terms():
            out = out * Monomial.term_coproduct(term)
            out.simplify()
        out.simplify()
        return out
            
    def reduced_coproduct(self):
        xx = self.coproduct()
        xx._reduce()
        return xx 
    
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

    @staticmethod
    def unit():
        return Polynomial([Monomial.unit()])

    @staticmethod
    def null():
        return Polynomial([])
    
    def __init__(self, summands):
        self.summands = summands

    def __str__(self):
        string = ""
        if not self.get_summands():
            return "0"
        for monomial in self.get_summands():
            string += str(monomial) + " + "
        return string[:-3]

    def __setstate__(self, _dict):
        self.summands = _dict['summands']
    
    def __getstate__(self):
        return {'summands' : self.summands }
               
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
                coefficient = (coefficient + x.get_coeff()) \
                              % opts.prime
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

    This needs to allow for terms to have exponentials!
    """
    @staticmethod
    def random(x):
        return TensorMonomial([Monomial.random(x), Monomial.random(x)], random.randrange(opts.prime))

    @staticmethod
    def unit():
        return TensorMonomial([Monomial.unit(), Monomial.unit()], 1)

    @staticmethod
    def null():
        return TensorMonomial([Monomial.unit(), Monomial.unit()], 0)
    
    def __init__(self, pair, coeff):
        self.pair  = pair
        self.coeff = coeff
        self.simplify_flag = False

    def __getstate__(self):
        return {'pair' : self.pair, 'coeff' : self.coeff}

    def __setstate__(self, _dict):
        self.pair = _dict['pair']
        self.coeff = _dict['coeff']
    
    def get_pair(self):
        return self.pair

    def get_coeff(self):
        return self.coeff
    
    def get_degree(self):
        adeg = self.pair[0].get_degree()
        bdeg = self.pair[1].get_degree()
        if self.coeff:
            return (adeg[0] + bdeg[0], adeg[1] + bdeg[1])
        else:
            return (0,0)

    def __str__(self):
        if not self.coeff:
            return "0"
        return str(self.coeff) + "*" \
            + str(self.pair[0]) + " | " \
            + str(self.pair[1])

    def simplify(self):
        """
        Warning! Mutator!
        """
        if self.simplify_flag:
            return
        self.pair[0].simplify()
        self.pair[1].simplify()
        newcoeff = (self.pair[0].coeff * self.pair[1].coeff * self.coeff) % opts.prime
        self.pair[0].coeff = 1
        self.pair[1].coeff = 1
        self.coeff = newcoeff
        if not self.coeff:
            self.pair[0] = Monomial([],1)
            self.pair[1] = Monomial([],1)
        self.simplify_flag = True

    def __mul__(self, other):
        newpair = [self.pair[0] * other.pair[0], self.pair[1] * other.pair[1] ]
        sign_fix = self.pair[1].get_degree()[0] *other.pair[0].get_degree()[0] 
        newcoeff = ((-1)**sign_fix * self.coeff * other.coeff) % opts.prime
        out = TensorMonomial(newpair, newcoeff)
        out.simplify()
        return out

    def __eq__(self, other):
        """
        Warning! Mutator...
        """
        self.simplify()
        other.simplify()
        if self.coeff != other.coeff:
#            print "coeffs diff"
            return False
        elif not self.pair[0] == other.pair[0]:
#            print "first terms diff"
#            print self.pair[0].terms, self.pair[0].coeff
#            print other.pair[0].terms, self.pair[0].coeff
            return False
        elif not self.pair[1] == other.pair[1]:
#            print "second terms diff"
            return False
        return True

    def monomials_equal(self, other):
        """Warning! mutator!  This will check if the monomials are equal,
        ignoring the coefficients, unless they are 0.

        """
        self.simplify()
        other.simplify()
        if not self.pair[0] == other.pair[0]:
            return False
        elif not self.pair[1] == other.pair[1]:
            return False
        return True

    def _reduce(self):
        self.simplify()
        for mon in self.get_pair():
            if mon == Monomial.unit():
                self.coeff = 0
            
    
class TensorPolynomial(object):
    """
    summands is a list of TensorMonomial objects
    """
    @staticmethod
    def random(x, y):
        summands = [ TensorMonomial.random(y) for i in range(x)]
        return TensorPolynomial(summands)

    @staticmethod
    def unit():
        return TensorPolynomial([TensorMonomial.unit()])

    @staticmethod
    def null():
        return TensorPolynomial([])
    
    def __init__(self, summands):
        self.summands = summands
        self.simplify_flag = False

    def __getstate__(self):
        return { 'summands' : self.summands}

    def __setstate__(self, _dict):
        self.summands = _dict['summands']

    def get_summands(self):
        return self.summands

    def stupid_simplify(self):
        for mon in self.summands[:]:
            
            if mon == TensorMonomial([Monomial([],0),
                                      Monomial([],0)],0):
                self.summands.remove(mon)

    def simplify(self):
        if self.simplify_flag:
            return
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
                coefficient = (coefficient + x.get_coeff()) \
                              % opts.prime
                stack.remove(x)
            if coefficient:
                outsum.append(TensorMonomial(mon.get_pair(), \
                                             coefficient))
        self.summands = outsum
        self.simplify_flag = True

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

    def _reduce(self):
        for tmon in self.summands:
            tmon._reduce()
        self.simplify()

#######################################################
        
# class CobarMonomial(object):
#     """
#     This is a generalized TensorMonomial object 
#     which allows for more tensor factors
#     """
#     def __init__(self, factors, coeff):
#         self.factors = factors
#         self.coeff = coeff
        
# class CobarPolynomial(object):
#     """
#     This is a generalized TensorPolynomial object
#     which allows for more tensor factors.

#     As we may use reduced cobar construction, should 
#     can perform reduced calculations.
#     """
#     def __init__(self, summands):
#         self.summands = summands

    
############################################################
    
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
            sloc = perm1[loc:].index(p0)+loc # Find position in perm1
            perm1[loc], perm1[sloc] = p0, p1  # Swap in perm1
            transCount += 1
                            
    # Even number of transpositions means equal parity
    if (transCount % 2) == 0:
        return True
    else:
        return False

def xi_deg(J):
    """
    J is a list of natural numbers
    The first term in list J corresponds to xi_1
    """
    out = 0
    for i in range(len(J)):
        out += J[i]*(2*(opts.prime**(i+1))-2)
    return out

def xi_wt(J):
    """
    J is a list of natural numbers
    returns weight of xi(J)
    """
    out = 0
    for i in range(len(J)):
        out += J[i]*(opts.prime**(i+1)-1)
    return out

def tau_deg(I):
    """
    J is a list of natural numbers
    First term in list I corresponds to tau_0
    """
    out = 0
    for i in range(len(I)):
        out += I[i]*(2*(opts.prime**(i))-1)
    return out

def tau_wt(I):
    out = 0
    for i in range(len(I)):
        out += I[i]*(opts.prime**i - 1)
    return out
    
def xi_indices():
    """
    Returns a list of all allowable indices J for which x(J) has 
    degree less than degree bounds
    """
    indices = []
    maximum = int(log((opts.bounds + 3)/2.0, opts.prime))
    for i in range(maximum):
        index = [0]*maximum
        index[i] = 1
        if xi_deg(index) <= opts.bounds:
            indices.append(index)
    j = maximum - 1
    while j >= 0:
        for x in indices[:]:
            y = x[:]
            y[j] += 1
            while xi_deg(y) <= opts.bounds and not (y in indices):
                indices.append(y[:])
                y[j] += 1
        j += (-1)
    sort_ind = sorted(indices, key=lambda y: xi_deg(y))
    return sort_ind

def tau_indices():
    """Returns a list of all allowable indices I for which t(I) has
    degree less than degree bounds

    """
    indices = []
    maximum = int(log((opts.bounds + 2)/2.0, opts.prime))
    for i in range(maximum + 1):
        index = [0]*(maximum+1)
        index[i] = 1
        if tau_deg(index) <= opts.bounds:
            indices.append(index)
    j = maximum
    while j >= 0:
        for x in indices[:]:
            y = x[:]
            if y[j] == 0:
                y[j] = 1
                if tau_deg(y) <= opts.bounds and not( y in indices):
                    indices.append(y[:])
        j += (-1)
    sort_ind = sorted(indices, key=lambda y: tau_deg(y))
    return sort_ind
            
def all_indices():
    """This returns a list of pairs (I,J) where t(I)x(J) is within degree
    bounds.

    """
    out = []
    all_I = tau_indices()
    #add in 0 list
    xx = [0]*(len(all_I[0]))
    all_I.insert(0, xx)
    all_J = xi_indices()
    yy = [0]*(len(all_J[0]))
    all_J.insert(0,yy)
    #Now go through in order and try to pair a tau index I with a xi
    #index J, until it fails. Then move on to next tau index I.
    for I in all_I:
        for J in all_J:
            if tau_deg(I) + xi_deg(J) <= opts.bounds:
                out.append((I,J))
            else:
                break
    sort_out = sorted(out, key=lambda pair: \
                      tau_deg(pair[0]) + xi_deg(pair[1]))
    #this needs to go back for reduced list
    sort_out.remove((xx, yy))
    return sort_out

    
