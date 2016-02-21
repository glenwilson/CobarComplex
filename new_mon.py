from options import opts
from math import log
import random 
import copy

class Monomial(object):
    """Declare a monomial with two lists of natural numbers I, J, and an
    additional tuple (e, i) corresponding to gamma^e zeta^i.

    """
    @staticmethod
    def random(length):
        randtaus = [ random.randrange(2) for i in
                     range(length)]
        randxis = [ random.randrange(3) for i in
                     range(length)]
        randcohom = [random.randrange(2), random.randrange(4) ]
        randcoeff = random.randrange(opts.prime)
        return Monomial(randtaus, randxis, randcohom, randcoeff)

    @staticmethod
    def random_no_cohom(length):
        randtaus = [ random.randrange(2) for i in
                     range(length)]
        randxis = [ random.randrange(3) for i in
                     range(length)]
        randcohom = [0,0]
        randcoeff = random.randrange(opts.prime)
        return Monomial(randtaus, randxis, randcohom, randcoeff)

    @staticmethod
    def unit():
        return Monomial([], [], [0,0], 1)

    @staticmethod
    def null():
        return Monomial([], [], [0,0], 0)
    
    @staticmethod
    def term_coproduct(term):
        """
        term must be a tuple of the form ("x", i) or ("t", i)
        """
        if term[0] == "x":
            if term[1] == 0:
                return TensorPolynomial.unit()
            out = TensorPolynomial.null()
            out += TensorPolynomial([TensorMonomial([Monomial.unit(), Monomial.term(term)], [0,0], 1)])
            out += TensorPolynomial([TensorMonomial([Monomial.term(term), Monomial.unit()], [0,0], 1)])
            for j in range(1, term[1]):
                out += TensorPolynomial([TensorMonomial([Monomial.term(("x", term[1]-j))**(opts.prime**j),\
                                                         Monomial.term(("x", j))], [0,0], 1)])
            return out
        elif term[0] == "t":
            out = TensorPolynomial.null()
            out += TensorPolynomial([TensorMonomial([Monomial.unit(), Monomial.term(term)],[0,0], 1)])
            out += TensorPolynomial([TensorMonomial([Monomial.term(term), Monomial.unit()],[0,0], 1)])
            for j in range(0, term[1]):
                out += TensorPolynomial([TensorMonomial([Monomial.term(("x", term[1]-j))**(opts.prime**j),\
                                                         Monomial.term(("t", j))], [0,0], 1)])
            return out

    @staticmethod
    def tau_list(indices):
        return Monomial(indices, [], [0,0], 1)

    @staticmethod
    def xi_list(indices):
        return Monomial([], indices, [0,0], 1)

    @staticmethod
    def tau_xi_list(pair):
        return Monomial(pair[0], pair[1], [0,0], 1)

    @staticmethod
    def term(term):
        """Term is either (t, i) or (x, j)  for i>=0 and j>= 1

        """
        if term[0] == "t":
            I = [0]*(term[1]+1)
            I[term[1]]=1
            return Monomial(I, [], [0,0], 1)
        elif term[0] == "x":
            J = [0]*(term[1])
            J[term[1]-1] = 1
            return Monomial([], J, [0,0], 1)
        else:
            raise TypeError
    
    @staticmethod
    def cohom_list(indices):
        return Monomial([], [], indices, 1)
    
    def __init__(self, taus, xis, cohom, coeff):
        """
        Taus is a list I = [i0, i1, ..., is] corr. to the monomial tau(I)

        Xis is a list J = [j1, j2, ..., jk] corr. to the monomial xi(J)

        cohom is a list [e, i] corr to gamma^e zeta^i

        coeff is a residue mod opts.prime
        """
        
        self.taus = taus
        self.xis = xis
        self.cohom = cohom
        self.coeff = coeff % opts.prime
        self.simplify_flag = False

    # def __getstate__(self):
    #     return {'terms' : self.terms,
    #             'coeff' : self.coeff,
    #             'simplify_flag' : self.simplify_flag}
        
    # def __setstate__(self, _dict):
    #     self.terms = _dict['terms']
    #     self.coeff = _dict['coeff']
    #     self.simplify_flag = _dict['simplify_flag']

    # def get_terms(self):
    #     return self.terms

    def get_coeff(self):
        return self.coeff 

    def get_degree(self):
        """
        """
        self.simplify()
        deg = 0
        wt = 0
        deg += xi_deg(self.xis)
        wt += xi_wt(self.xis)
        deg += tau_deg(self.taus)
        wt+= tau_wt(self.taus)
        deg += self.cohom[0]*(-1)
        wt += self.cohom[0]*(-1)*opts.ind
        wt += self.cohom[1]*(-1)*opts.ind
        return (deg, wt)

    def copy(self):
        newtaus = self.taus[:]
        newxis = self.xis[:]
        newcohom = self.cohom[:]
        newcoeff = self.coeff
        return Monomial(newtaus, newxis, newcohom, newcoeff)
                
    def __str__(self):
        mystr = ""
        mystr += str(self.get_coeff()) 
        if self.cohom[0]:
            mystr += ".g"
        if self.cohom[1]:
            mystr += ".z^" + str(self.cohom[1]) 
        for i in range(len(self.taus)):
            if self.taus[i]:
                mystr += ".t" + str(i) + "^" + str(self.taus[i])
        for i in range(len(self.xis)):
            if self.xis[i]:
                mystr += ".x" + str(i+1) + "^" + str(self.xis[i])
        return mystr

    def simplify(self):
        if self.simplify_flag:
            return
        if (not self.coeff) or (self.cohom[0] >= 2):
            self.taus = []
            self.xis = []
            self.cohom = [0,0]
            self.coeff = 0
        self.simplify_flag = True

    def __rmul__(self, n):
        if isinstance(n, int):
            out = self.copy()
            out.coeff = (n * out.coeff) % opts.prime
            return out
        else:
            raise TypeError("Expecting an integer")
        
    def __mul__(self, other):
        """
        Does not mutate self or other. 
        Returns the product
        """
        newcoeff = (self.coeff * other.coeff) % opts.prime
        if not newcoeff:
            return Monomial.null()
        newcohom = []
        if self.cohom[0] + other.cohom[0] >= 2:
            return Monomial.null()
        newcohom.append(self.cohom[0] + other.cohom[0])
        newcohom.append(self.cohom[1] + other.cohom[1])
        maxtaus = max(len(self.taus), len(other.taus))
        newtaus = []
        for i in xrange(maxtaus):
            try:
                n = self.taus[i]
            except IndexError:
                n = 0
            try:
                m = other.taus[i]
            except IndexError:
                m = 0
            if n + m == 2:
                return Monomial.null()
            newtaus.append(n+m)
        maxxis = max(len(self.xis), len(other.xis))
        newxis = []
        for i in xrange(maxxis):
            try:
                n = self.xis[i]
            except IndexError:
                n = 0
            try:
                m = other.xis[i]
            except IndexError:
                m = 0
            newxis.append(n+m)
        selftauindex = [i for i in range(len(self.taus)) if
                        self.taus[i] == 1]
        othertauindex = [i for i in range(len(other.taus)) if
                         other.taus[i] == 1]
        perm1 = selftauindex + othertauindex
        perm2 = sorted(perm1)
        if not arePermsEqualParity(perm1, perm2):
            newcoeff = (-newcoeff) % opts.prime
        if other.cohom[0] == 1:
            newcoeff = ((-1)**(sum(self.taus)) * newcoeff) % opts.prime
        return Monomial(newtaus, newxis, newcohom, newcoeff)
        
    def __pow__(self, num):
        out = Monomial.unit()
        for i in xrange(num):
            out = out * self
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
        if self.cohom != other.cohom:
            return False
        maxtaus = max(len(self.taus), len(other.taus))
        for i in xrange(maxtaus):
            try:
                n = self.taus[i]
            except IndexError:
                n = 0
            try:
                m = other.taus[i]
            except IndexError:
                m = 0
            if n != m:
                return False
        maxxis = max(len(self.xis), len(other.xis))
        for i in xrange(maxxis):
            try:
                n = self.xis[i]
            except IndexError:
                n = 0
            try:
                m = other.xis[i]
            except IndexError:
                m = 0
            if n != m:
                return False
        return True

    def __neq__(self, other):
        return not (self == other)

    def monomials_equal(self, other):
        self.simplify()
        other.simplify()
        if self.cohom != other.cohom:
            return False
        maxtaus = max(len(self.taus), len(other.taus))
        for i in xrange(maxtaus):
            try:
                n = self.taus[i]
            except IndexError:
                n = 0
            try:
                m = other.taus[i]
            except IndexError:
                m = 0
            if n != m:
                return False
        maxxis = max(len(self.xis), len(other.xis))
        for i in xrange(maxxis):
            try:
                n = self.xis[i]
            except IndexError:
                n = 0
            try:
                m = other.xis[i]
            except IndexError:
                m = 0
            if n != m:
                return False
        return True
    
    def coproduct(self):
        coeff = self.coeff
        cohom = self.cohom[:]
        out = TensorPolynomial([ TensorMonomial([Monomial.unit(),
                                                 Monomial.unit()],
                                                cohom, coeff)])
        for i in range(len(self.taus)):
            if self.taus[i]:
                out = out * Monomial.term_coproduct(("t", i))
        for i in range(len(self.xis)):
            if self.xis[i]:
                out = out * (Monomial.term_coproduct(("x", i+1))**(self.xis[i]))
        out.simplify()
        return out
            
    # def reduced_coproduct(self):
    #     xx = self.coproduct()
    #     xx._reduce()
    #     return xx 
    
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
        self.simplify_flag = False

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

    def __rmul(self, n):
        if isinstance(n, int):
            newsummands = []
            for mon in self.get_summands():
                newsummands.append(n*mon)
            return Polynomial(newsummands)
        else:
            raise TypeError("Expecting integer")
                
    
    def __mul__(self, other):
        newsummands = []
        for mon1 in self.get_summands():
            for mon2 in other.get_summands():
                newsummands.append(mon1*mon2)
        return Polynomial(newsummands)
            
    def stupid_simplify(self):
        for monomial in self.get_summands():
            monomial.simplify()
        for monomial in self.get_summands()[:]:
            if monomial == Monomial.null():
                self.summands.remove(monomial)            
    
    def simplify(self):
        """
        Mutator!
        """
        if self.simplify_flag:
            return
        print "poly simplify"
        self.stupid_simplify()
        stack = self.get_summands()[:]
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
                moncopy = mon.copy()
                moncopy.coeff = coefficient
                outsum.append(moncopy)
        self.summands = outsum
        self.simplify_flag = True

    def copy(self):
        return Polynomial([mon.copy() for mon in self.summands])

        
    def __eq__(self, other):
        """
        This should maybe mutate monomials to simplified form?
        Currently uses costly copies
        """
        self.simplify()
        other.simplify()
        stack1 = self.get_summands()[:]
        stack2 = other.get_summands()[:]
        for x in self.get_summands()[:]:
            if x in stack2:
                stack1.remove(x)
                stack2.remove(x)
            else:
                return False
        if len(stack1) or len(stack2):
            return False
        else:
            return True

    def __neq__(self, other):
        return not (self == other)

class TensorMonomial(object):
    """
    pair is a list [x,y] with x and y Monomials, 
    corresponds to x \otimes y

    This needs to allow for terms to have exponentials!
    """
    @staticmethod
    def random(x):
        randcoeff = random.randrange(opts.prime)
        randcohom = [random.randrange(2), random.randrange(4) ]
        xx = Monomial.random(x)
        yy = Monomial.random(x)
        xx.cohom = [0,0]
        yy.cohom = [0,0]
        return TensorMonomial([xx, yy], randcohom, randcoeff)

    @staticmethod
    def unit():
        return TensorMonomial([Monomial.unit(), Monomial.unit()], [0,0], 1)

    @staticmethod
    def null():
        """
        Should the tensor factors be MOnomial.null()? Careful! 
        """
        return TensorMonomial([Monomial.unit(), Monomial.unit()], [0,0], 0)
    
    def __init__(self, pair, cohom, coeff):
        """Pair is a list consisting of two Monomial objects. The monomial
objects should not have any terms from mot. coh., i.e., cohom = [0,0].

        cohom = [e, i] describes the coefficient of motivic cohomology
        (mult. on left)

        coeff is the coefficient

        """
        self.pair  = pair
        self.cohom = cohom                      
        self.coeff = coeff
        self.simplify_flag = False

    # def __getstate__(self):
    #     return {'pair' : self.pair, 'coeff' : self.coeff}

    # def __setstate__(self, _dict):
    #     self.pair = _dict['pair']
    #     self.coeff = _dict['coeff']
    
    def get_pair(self):
        return self.pair

    def get_coeff(self):
        return self.coeff

    def get_cohom(self):
        return self.cohom

    def get_degree(self):
        adeg = self.pair[0].get_degree()
        bdeg = self.pair[1].get_degree()
        cohomdeg = cohom_bideg(self.cohom)
        if self.coeff:
            return (adeg[0] + bdeg[0] + cohomdeg[0], adeg[1] + bdeg[1] + cohomdeg[1])
        else:
            return (0,0)

    def __str__(self):
        if not self.coeff:
            return "0"
        out = str(self.coeff) 
        if self.cohom[0]:
            out += ".g"
        if self.cohom[1]:
            out += ".z^" + str(self.cohom[1])
        out += "." + str(self.pair[0]) + " | " + str(self.pair[1])
        return out

    def copy(self):
        newpair = [mon.copy() for mon in self.pair]
        newcohom = self.cohom[:]
        newcoeff = self.coeff
        return TensorMonomial(newpair, newcohom, newcoeff)
    
    def simplify(self):
        """This will become more complicated. Will need to return a tensor
        polynomial object to be entirely general.

        """
        if self.simplify_flag:
            return
        self.pair[0].simplify()
        self.pair[1].simplify()
        newcoeff = (self.pair[0].coeff * self.pair[1].coeff * self.coeff) % opts.prime
        self.pair[0].coeff = 1
        self.pair[1].coeff = 1
        self.coeff = newcoeff
        if self.pair[0].cohom != [0,0] or self.pair[1].cohom != [0,0]:
            raise TypeError("Can't handle cohom in pair yet")
        if not self.coeff:
            self.pair[0] = Monomial.unit()
            self.pair[1] = Monomial.unit()
        self.simplify_flag = True

    def __rmul__(self, n):
        if isinstance(n, int):
            out = self.copy()
            out.coeff = (n * out.coeff) % opts.prime
            return out
        else:
            raise TypeError("Expecting an integer")

        
    def __mul__(self, other):
        newpair = [self.pair[0] * other.pair[0], self.pair[1] * other.pair[1] ]
        sign_fix = self.pair[1].get_degree()[0] *(other.pair[0].get_degree()[0] + cohom_deg(other.cohom))
        newcoeff = ((-1)**sign_fix * self.coeff * other.coeff) % opts.prime
        newcohom = [self.cohom[0] + other.cohom[0], self.cohom[1] +
                    other.cohom[1]]
        
        out = TensorMonomial(newpair, newcohom, newcoeff)
        out.simplify()
        return out

    def __eq__(self, other):
        """
        Warning! Mutator...
        """
        self.simplify()
        other.simplify()
        if self.coeff != other.coeff:
            return False
        elif not self.pair[0] == other.pair[0]:
            return False
        elif not self.pair[1] == other.pair[1]:
            return False
        elif self.cohom != other.cohom:
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
        elif not self.cohom == other.cohom:
            return False
        return True

    # def _reduce(self):
    #     self.simplify()
    #     for mon in self.get_pair():
    #         if mon == Monomial.unit():
    #             self.coeff = 0
            
    
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

    def copy(self):
        return TensorPolynomial([mon.copy() for mon in self.summands])
        
    def get_summands(self):
        return self.summands

    def stupid_simplify(self):
        for mon in self.summands:
            mon.simplify()
        for mon in self.summands[:]:
            if mon == TensorMonomial.null():
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
                outsum.append(TensorMonomial(mon.get_pair(),
                                             mon.get_cohom(),
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

    def __rmul__(self, n):
        if isinstance(n, int):
            newsummands = []
            for mon in self.get_summands():
                newsummands.append(n*mon)
            return Polynomial(newsummands)
        else:
            raise TypeError("Expecting integer")

    def __mul__(self, other):
        newsummands = []
        for mon1 in self.get_summands():
            for mon2 in other.get_summands():
                newsummands.append(mon1*mon2)
        return TensorPolynomial(newsummands)

    def __pow__(self, n):
        out = TensorPolynomial.unit()
        for i in xrange(n):
            out = out * self
        return out

        
    # def _reduce(self):
    #     for tmon in self.summands:
    #         tmon._reduce()
    #     self.simplify()

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

def cohom_deg(cohom):
    return -cohom[0]

def cohom_wt(cohom):
    return -(cohom[0] + cohom[1])*opts.ind    

def cohom_bideg(cohom):
    return (cohom_deg(cohom), cohom_wt(cohom))

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
    #sort_out.remove((xx, yy))
    return sort_out

    
