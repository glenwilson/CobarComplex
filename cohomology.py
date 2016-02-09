from mod_lin_alg import ModMatrix, ModVector
from vect_space import ModVectorSpace

class Cohomology(object):
    """This class defines a cohomology object. It is initialized with two
    ModMatrix objects B, A such that B*A is defined. We assume that
    both matrices are expressed in terms of a common basis in
    codomain(A) = domain(B).

    This class has methods which will determine a basis for Ker(B), a
    basis for Im(A), and for Ker(B)/Im(A). There is a method to
    determine if two cocycles are equivalent in cohomology.

    Now will have an attribute to store a basis for cohomology in
    terms of product structure.

    """
    
    def __init__(self, B, A):
        self.B = B
        self.A = A
        self.kernel = None
        self.kernel_flag = False
        self.image = None
        self.image_flag = False
        self.cohomology = None
        self.cohomology_flag = False
        self.product = ModVectorSpace([])
        # if self.A.rows == []:
        #     self.A = ModMatrix.null(self.B.get_size()[1], 1)
        # elif self.B.rows == []:
        #     raise TypeError
        if self.B.get_size()[1] != self.A.get_size()[0]:
            raise TypeError("wrong sizes B: " + str(self.B) + "A: " +
                            str(self.A))
        self.extended_ker_basis = None
        self.extended_ker_basis_flag = False
    
    def get_A(self):
        return self.A

    def get_B(self):
        return self.B

    def get_kernel_flag(self):
        return self.kernel_flag

    def get_image_flag(self):
        return self.image_flag

    def get_cohomology_flag(self):
        return self.cohomology_flag

    def get_zero_vector(self):
        dim = self.get_A().get_size()[0]
        return ModVector.null(dim)

    def im_in_ker(self):
        C = self.get_B() * self.get_A()
        if not C.is_zero():
            return False
        else:
            return True
    
    def get_kernel(self):
        if not self.get_kernel_flag():
            self.kernel = ModVectorSpace(self.get_B().get_kernel_vout())
            self.kernel_flag = True
        return self.kernel

    def in_kernel(self, vect):
        """
        Is vect a ModMatrix or a ModVector??? Fix this
        """
        C = self.get_B() * vect
        return C.is_zero()

    def compute_image(self):
        image = ModVectorSpace([])
        for row in self.get_A().get_rref():
            if not row.is_zero():
                index = row.get_leading_index()
                image.add_basis_element([ self.get_A().get_column(index) ])
        self.image = image
        self.image_flag = True
    
    def get_image(self):
        if not self.get_image_flag():
            self.compute_image()
        return self.image

    def are_cohomologous(self, vect_1, vect_2):
        """
        inputs are bitvectors, not bitmatrix object representing col. vector
        """
        vector_sum = vect_1 + ( (-1) * vect_2)
        return self.get_A().can_solve(vector_sum)

    def compute_cohomology(self):
        cohomology = ModVectorSpace([])
        index_list = self.get_kernel().get_basis()[:]
        basis_list = []
        for ker_element in index_list:
            if not self.in_cohomology_subspace(basis_list, ker_element):
                basis_list.append(ker_element)
        cohomology.add_basis_element(basis_list)
        self.cohomology = cohomology
        if (self.get_B().col_count - self.get_B().get_rank() 
            - self.get_A().get_rank() != len(cohomology.get_basis())):
            print "A"
            print self.get_A()
            print "B"
            print self.get_B()
            print "rank A"
            print self.get_A().get_rank()
            print "nullity B"
            print self.get_B().col_count - self.get_B().get_rank() 
            print "basis for cohomology"
            for thing in cohomology.get_basis():
                print thing
            raise TypeError("cohomology basis computation incorrect!")
        self.cohomology_flag = True

    def get_cohomology(self):
        if not self.get_cohomology_flag():
            if not self.im_in_ker():
                print "matrix for image"
                print str(self.A)
                print "matrix for kernel"
                print str(self.B)
                raise TypeError("Not exact here!")
            self.compute_cohomology()
        return self.cohomology

    def get_product(self):
        return self.product

    def add_to_product_basis(self, a_list):
        """
        input is a list of bitvectors
        """
        self.product.add_basis_element(a_list)

    def in_cohomology_subspace(self, subspace_basis, vector):
        """
        subspace_basis should be given as bitvectors
        vector is also a bit vector
        all should be in terms of standard basis for 
        codomain(A) = domain(B)
        """
        B = self.get_A().get_append_columns(subspace_basis)
#        print "checking if\n", vector, "\nin subspace\n", B
#        print "result", B.can_solve(vector)
        return B.can_solve(vector)
                
    # def coset(self, vector):
    #     """
    #     given a vector V, returns the set 
    #     V + Im(A)
    #     """
    #     out_set = set()
    #     im = self.get_image().get_basis()
    #     coeffs = BitVector.vector_space_set(len(im))
    #     for cc in coeffs:
    #         new = vector.copy()
    #         for i in range(0,len(im)):
    #             if cc[i]:
    #                 new += im[i]
    #         out_set.add(new)
    #     return out_set

    # def short_representative(self, vector):
    #     coset = self.coset(vector)
    #     min_cpt = vector.number_of_components()
    #     out = vector
    #     for w in coset:
    #         if w.number_of_components() < min_cpt:
    #             min_cpt = w.number_of_components()
    #             out = w
    #     return out

    def extend_ker_basis(self):
        ker_basis = self.get_kernel().get_basis()
        #assuming vectors are bit vectors
        if not ker_basis:
            C = ModMatrix.null(1, self.B.col_count)
        else:
            C = ModMatrix(ker_basis)
        comp = C.complement_row_space()
        tot = ker_basis[:] + comp[:]
        self.extended_ker_basis = tot
        self.extended_ker_basis_flag = True

    def get_extended_ker_basis(self):
        if not self.extended_ker_basis_flag:
            self.extend_ker_basis()
        return self.extended_ker_basis
            
    def basis_to_ker_basis(self):
        C = ModMatrix(self.get_extended_ker_basis())
        C.transpose()
        inv = C.get_inverse()
        if not inv:
            raise TypeError("the matrix is not invertible! Something went very wrong")
        return inv
