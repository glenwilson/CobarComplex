from scratch import Options, Monomial, Polynomial, opts, TensorMonomial
import random

opts.prime =7

def simplify_test(length, n):
    """
    length is length of monomials
    n dictates how many tests are done
    """
    for i in range(n):
        xx = Monomial.random(length)
        print xx
        xx.simplify()
        print xx
        print "\n"
    
#simplify_test(10,5)


def product_test(length, n):
    for i in range(n):
        xx = Monomial.random(length)
        yy = Monomial.random(length)
        zz = xx*yy
        tt = yy*xx
        xx.simplify()
        yy.simplify()
        print xx
        print yy
        print zz
        print tt
        print "\n"

#product_test(4, 10)


def equality_test(length, n):
    for i in range(n):
        xx = Monomial.random(length)
        yy = Monomial.random(length)
        print xx, yy, xx==yy

#equality_test(2, 100)

def degree_test(length, n):
    for i in range(n):
        xx = Monomial.random(length)
        print xx
        print xx.get_degree()

#print opts.prime
#degree_test(2,10)

def homogeneous_test(x,y,n):
    for i in range(n):
        xx = Polynomial.random(x,y)
        xx.simplify()
        print xx, " ... ", str(xx.homogeneous())

#homogeneous_test(2,2,100)
        
def poly_simplify_test(x,y,n):
    for i in range(n):
        xx = Polynomial.random(x,y)
        print xx
        xx.stupid_simplify()
        print str(xx) + "\n"

#print opts.prime
#poly_simplify_test(4,2,10)

def poly_simplify_test(x,y,n):
    for i in range(n):
        xx = Polynomial.random(x,y)
        yy = xx + xx
        print yy
        yy.simplify()
        print str(yy) + "\n"

#poly_simplify_test(2,3,20)

def poly_test_eq(x,y,n):
    for i in range(n):
        xx = Polynomial.random(x,y)
        yy = Polynomial.random(x,y)
        print xx, "...", yy, "...",  xx==yy

#poly_test_eq(2,2,1000)

def tensmon_simplify(x,y,n):
    for i in range(n):
        xx = Monomial.random(x)
        yy = Monomial.random(y)
        zz = TensorMonomial([xx,yy], random.randrange(opts.prime))
        print zz
        zz.simplify()
        print zz

tensmon_simplify(3,3,100)
