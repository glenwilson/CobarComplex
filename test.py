from options import opts
#from dual_st_alg import Monomial, Polynomial, opts, TensorMonomial, TensorPolynomial
from dual_st_alg import *
from mod_lin_alg import ModVector, ModMatrix
from cohomology import Cohomology
from cobar_complex import CobarMonomial, CobarPolynomial, CobarModule, CobarComplex
import random

#opts.prime = 3
opts.prime = 3
#opts.prime = 13
#opts.prime = 541

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

#tensmon_simplify(3,3,100)

def tenspoly_str(x,y,n):
    for i in range(n):
        xx = TensorPolynomial.random(x,y)
        print xx

#tenspoly_str(2,3,10)

def tenspoly_sum(x,y,n):
    for i in range(n):
        xx = TensorPolynomial.random(x,y)
        yy = TensorPolynomial.random(x,y)
        zz = xx + yy
        print xx
        print yy
        print str(zz) + "\n"

#tenspoly_sum(3,3,10)

def tensmon_eq(x,n):
    for i in range(n):
        xx = TensorMonomial.random(x)
        print xx
        yy =TensorMonomial.random(x)
        print yy
        print xx == yy 

        
#tensmon_eq(3,100)

def tenspoly_simplify(x,y,n):
    for i in range(n):
        xx = TensorPolynomial.random(x,y)
        print xx
        xx.simplify()
        print xx

#tenspoly_simplify(3,3,100)

def tenspoly_add(x,y,n):
    for i in range(n):
        xx = TensorPolynomial.random(x,y)
        yy = TensorPolynomial.random(x,y)
        print xx
        print yy
        print xx + yy + xx
        zz = xx + yy + xx
        zz.simplify()
        print zz, "\n"

#tenspoly_add(1,1,100)

def tenspoly_mult(x,y,n):
    for i in range(n):
        xx = TensorPolynomial.random(x,y)
        yy = TensorPolynomial.random(x,y)
        print xx
        print yy
        zz = xx * yy
        zz.simplify()
        print zz, "\n"

#tenspoly_mult(3,3,100)

def coprod_term(n):
    for i in range(1,n):
        print Monomial.term_coproduct(("x", i))
    for i in range(n):
        print Monomial.term_coproduct(("t", i))

#coprod_term(4)

def coprod_test(x,n):
    for i in range(n):
        xx = Monomial.random(x)
        print xx
        print xx.coproduct(), "\n"
    
#coprod_test(2, 5)

def reduced_coprod_test(x,n):
    for i in range(n):
        xx = Monomial.random(x)
        print xx
        print "coprod",  xx.coproduct(), "\n"
        print "reduced", xx.reduced_coproduct(), "\n"

#reduced_coprod_test(2,5)


#####################################################

def vector_add(x, n):
    for i in range(n):
        xx = ModVector.random(x)
        yy = ModVector.random(x)
        print xx
        print yy
        print xx+yy, "\n"

# opts.prime = 13
# print opts.prime
# vector_add(5, 10)

def vector_mul(x,n):
    for i in range(n):
        xx = ModVector.random(x)
        yy = ModVector.random(x)
        print xx
        print yy
        print xx * yy, "\n"

# opts.prime = 3
# print opts.prime
# vector_mul(5, 10)

def scalar_mul(x,n):
    for i in range(n):
        xx = ModVector.random(x)
        aa = random.randrange(opts.prime)
        print aa,  xx
        print aa * xx

# opts.prime = 5
# scalar_mul(3,10)
# vector_mul(3,10)

def leading_index(x,n):
    for i in range(n):
        xx = ModVector.random(x)
        print xx
        print xx.get_leading_index()
        print len(xx)

#leading_index(2,10)

# print ModMatrix.random(4,5)
# print ModMatrix.null(4, 4)
# print ModMatrix.identity(7)

def matrix_add(i, j, n):
    for a in range(n):
        xx = ModMatrix.random(i, j)
        yy = ModMatrix.random(i, j)
        print xx, "\n", yy, "\n",  xx + yy, "\n"
    xx = ModMatrix.random(i, j)
    yy = ModMatrix.random(i+1, j+1)
    try:
        print xx + yy
    except TypeError:
        print xx, "\n", "wrong size", "\n",  yy

#matrix_add(3,4,2)

def matrix_mul(i,j,k,n):
    for a in range(n):
        xx = ModMatrix.random(i,j)
        yy = ModMatrix.random(j,k)
        print xx, "\n\n", yy, "\n\n", xx*yy, "\n"

#matrix_mul(10,23,15,2)

def row_op_test(x,y,n):
    for i in range(n):
        xx = ModMatrix.random(x,y)
        print xx
        xx.el_row_op( random.randrange(opts.prime), 0, 1) 
        print xx
        
#row_op_test(2,4,5)

def rref_test(x,y,n):
    for i in range(n):
        xx = ModMatrix.random(x,y)
        xx.compute_rref()
        print xx, "\n"
        print xx.rref, "\n"
        print xx.basis_change, "\n"
        pp = xx.basis_change
        print pp * xx
        rr = pp * xx
        print rr == xx.rref

#rref_test(7,4,1)

def solve_test(x,y,n):
    for i in range(n):
        AA = ModMatrix.random(x,y)
        bb = ModVector.random(x)
        print AA.can_solve(bb)
        xx = AA.solve(bb)
        print "A\n",  AA, "\n"
        print "rrefA\n", AA.get_rref(), "\n"
        print "b\n", bb, "\n"
        print "x\n", xx, "\n"
        print "Ax\n", AA * xx
        print AA * xx == ModMatrix([bb]).get_transpose()

#solve_test(12,28,1)
# AA = ModMatrix([ModVector([0,1]), ModVector([2,2]), ModVector([1,2])])
# bb = ModVector([1,0,1])
# print AA, "\n"
# print bb, "\n"
# print "row weights\n", AA.get_row_weights(), "\n"
# print "rref\n", AA.get_rref(), "\n"
# print "basis change\n", AA.get_basis_change(), "\n"
# print AA.can_solve(bb)

def ker_test(x,y,n):
    for i in range(n):
        AA = ModMatrix.random(x,y)
        print AA.get_rref(), "\n"
        ker = AA.get_kernel()
        for xx in ker:
            print "soln\n", xx, "\n"
            print "result\n", AA * xx

#ker_test(3,4,1)

def rank_test(x,y,n):
    for i in range(n):
        AA = ModMatrix.random(x,y)
        print "AA\n", AA, "\n"
        print AA.get_rank()

#rank_test(3,3,5)

def inv_test(x,n):
    for i in range(n):
        AA = ModMatrix.random(x,x)
        print "AA\n", AA, "\n"
        print "inv\n", AA.get_inverse(), "\n"
        if AA.get_inverse():
            print "Ix\n", AA * AA.get_inverse(), "\n"

#inv_test(2,2)
#inv_test(3,2)

def append_test(x,y,n):
    for i in range(n):
        AA = ModMatrix.random(x,y)
        zz = ModVector.random(x)
        print "AA\n", AA, "\n"
        print "zz\n", zz, "\n"
        print "append\n", AA.get_append_columns([zz,zz + zz]), "\n"

#append_test(2,3,5)

def cohom_test(x,y,z,n):
    for i in range(n):
        AA = ModMatrix.random(y,z)
        BB = ModMatrix.random(x,y)

        if (BB * AA).is_zero():
            print"BB\n", BB, "\n"
            print "AA\n", AA, "\n"
            coh = Cohomology(BB, AA)
            print coh.get_cohomology()

#cohom_test(0,0,0,2)


def cobar_mon(x,f,n):
    for i in range(n):
        xx = CobarMonomial.random(x,f)
        print xx, "\n"
        print xx.get_degree(), "\n"

#cobar_mon(2,4,3)

def cobar_simp(x,f,n):
    for i in range(n):
        xx = CobarMonomial.random(x,f)
        print xx
        xx.simplify()
        print xx
        for term in xx.factors:
            print term

#cobar_simp(2,3,2)

def cobar_poly_simp(x,y,f,n):
    for i in range(n):
        xx = CobarPolynomial.random(x,y,f)
        yy = CobarPolynomial.random(x,y,f)
        zz = xx + yy
        print zz
        zz.simplify()
        print zz
        zz = zz + zz
        zz.simplify()
        print zz
        
#cobar_poly_simp(2,3,3,2)

def cobar_concat(x,y,f,n):
    for i in range(n):
        xx = CobarPolynomial.random(x,y,f)
        yy = CobarPolynomial.random(x,y,f)
        xx.simplify()
        yy.simplify()
        zz = xx & yy
        print xx
        print yy
        print zz, "\n"

#cobar_concat(2,2,2,2)

# opts.prime = 3
# opts.bounds = 20
# print opts.prime, opts.bounds
# print xi_indices()
# print [xi_deg(y) for y in xi_indices()]
# print tau_indices()
# print [tau_deg(y) for y in tau_indices()]
# new = sorted(tau_indices(), key=lambda y: tau_deg(y))
# print new
# print [tau_deg(y) for y in new]
# print "\n"
# print all_indices()
# print [tau_deg(p[0]) + xi_deg(p[1]) for p in all_indices()]

def cplx_test(length):
    C = CobarComplex(length)
    C.generate_modules()
    for i in range(C.length):
        print "module", i
        for deg in C.cplx[i]._dict:
            print "deg", deg
            for thing in C.cplx[i]._dict[deg]:
                print thing

#opts.prime = 3
#opts.bounds = 12
#cplx_test(7)

def map_test(f,x,n):
    for i in range(n):
        xx = CobarMonomial.random(x,f)
        print xx
        print xx._map_reduced()

#map_test(1,1,1)

def summand_test(x,y,f,n):
    for i in range(n):
        xx = CobarMonomial.random(x,f)
        yy = CobarPolynomial.random(x,y,f)
        xx.simplify()
        yy.simplify()
        print xx
        print yy
        print yy.is_summand(xx), "\n"

#summand_test(2,5,2,500)

def vfe_test(x,f,n):
    C= CobarComplex(7)
    print "calculating complex"
    cplx = C.get_cplx()
    for i in range(n):
        xx = CobarMonomial.random(x,f)
        print "finding good xx"
        while xx.get_degree()[0] > opts.bounds:
            xx = CobarMonomial.random(x,f)
        yy = xx._map()
        vect = C.vector_from_element(yy)
        deg = yy.get_degree()[0]
        wt = yy.get_degree()[1]
        for thing in C.get_cplx()[f+1]._dict[(deg,wt)]:
            print thing
        print "\n"
        print xx
        print yy
        print vect

#opts.prime=5
#opts.bounds=15
#vfe_test(1,1,10)

def cohom_test(f):
    C = CobarComplex(f+2)
    C.get_pickled_cplx()
    #C.make_maps()
    for i in range(1,f):
        for bideg in C.get_cplx()[i]._dict.keys():
            d = bideg[0]
            w = bideg[1]
            cohom = C.get_cohomology(i,d,w)
            print "Filt, deg, wt", i, d, w
            cc = cohom.get_cohomology()
            print cc
            if cc.get_basis():
                for vect in cc.get_basis():
                    print C.element_from_vector(vect, i,d,w)
    C.pickle_cplx()

    
opts.prime = 7
opts.bounds = 12
#cplx_test(3)
cohom_test(5)

#C = CobarComplex(7)
#C.make_map(1)
#C.make_map(2)
