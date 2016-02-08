from options import opts
import copy
import random

class ModVector(object):
    """
    Vectors are initialized with all entries 0. 
    """
    @staticmethod
    def random(length):
        x = ModVector([ random.randrange(opts.prime) for i in xrange(length) ])
        return x

    @staticmethod
    def null(length):
        return ModVector([0 for i in range(length)])
    
    def __str__(self):
        string = ""
        num_len = len(str(opts.prime))
        for i in self.vector:
            string += str(i).zfill(num_len) + " "
        return string
    
    def __init__(self, vector):
        """
        vector is a list of residues mod opts.prime
        """
        self.vector = vector
        self.length = len(vector)

    def __getitem__(self, position):
        return self.vector[position]

    def __setitem__(self, position, value):
        self.vector[position] = value
        
    def __add__(self, other):
        if self.length != other.length:
            raise TypeError("Wrong length")
        out = ModVector([(self.vector[i] + other.vector[i])\
                         % opts.prime for \
                         i in range(self.length)])
        return out

    def __mul__(self, other):
        """
        this is a dot product
        """
        if self.length != other.length:
            raise TypeError("WRong length")
        out = 0 
        for i in xrange(self.length):
            out = (out + self[i] * other[i]) % opts.prime
        return out

    def copy(self):
        return copy.deepcopy(self)
    
    def __rmul__(self, n):
        """
        this is scalar multiplication by n
        """
        out = self.copy()
        for i in xrange(self.length):
            out[i] = (n * self[i]) % opts.prime
        return out

    def get_leading_index(self):
        """
        Returns index of first non-zero entry.
        If the vector is null, then returns -1
        """
        for i in xrange(self.length):
            if (self[i] % opts.prime):
                return i
        return -1

    def __eq__(self, other):
        if self.length != other.length:
            return False
        for i in xrange(self.length):
            if (self[i] - other[i]) % opts.prime:
                return False
        return True

    def __len__(self):
        return len(self.vector)
    
    def is_zero(self):
        return self == ModVector.null(self.length)

    def __getstate__(self):
        return { 'vector' : self.vector,
                 'length' : self.length}

    def __setstate__(self, _dict):
        self.vector = _dict['vector']
        self.length = _dict['length']

class ModMatrix(object):
    @staticmethod
    def random(m, n):
        return ModMatrix([ModVector.random(n) for i in range(m)])

    @staticmethod
    def null(m, n):
        return ModMatrix([ModVector.null(n) for i in range(m)])

    @staticmethod
    def identity(n):
        out = ModMatrix.null(n,n)
        for i in range(n):
            out.rows[i][i] = 1
        return out
            
    def __init__(self, rows):
        self.rows = rows
        self.row_count = len(rows)
        self.col_count = 0
        if rows:
            self.col_count = len(rows[0])
            for row in rows:
                if len(row) != self.col_count:
                    print self.col_count
                    print len(row)
                    raise TypeError("length of rows inconsistent")
        self.basis_change = []
        self.basis_change_flag = False
        self.rref = []
        self.rref_flag = False
        self.row_weights = []
        self.row_weights_flag = False

    def __eq__(self, other):
        if self.row_count != other.row_count:
            return False
        for i in xrange(self.row_count):
            if not (self.rows[i] == other.rows[i]):
                return False
        return True
            
    def is_zero(self):
        return self == ModMatrix.null(self.row_count, self.col_count)
    
    def __str__(self):
        string = ""
        for row in self.rows:
            string += str(row) + "\n"
        return string[:-1]

    def get_size(self):
        return (self.row_count, self.col_count)
    
    def get_row(self, i):
        return self.rows[i].copy()

    def get_entry(self, i, j):
        return self.rows[i][j]

    def get_column(self, j):
        col = ModVector([self.get_entry(i, j) for \
                         i in xrange(self.row_count)])
        return col

    def copy(self):
        return copy.deepcopy(self)

    def __add__(self, other):
        if (self.row_count != other.row_count) or (self.col_count != other.col_count):
            raise TypeError
        out = []
        for i in range(self.row_count):
            out.append(self.get_row(i) + other.get_row(i))
        return ModMatrix(out)
            
    def __mul__(self, other):
        if self.col_count != other.row_count:
            raise TypeError("not compatible for mult")
        out = []
        for i in xrange(self.row_count):
            out.append(ModVector([self.get_row(i) \
                                  * other.get_column(j) \
                                  for j in \
                                  xrange(other.col_count) ]))
        return ModMatrix(out)

    def get_transpose(self):
        new_rows = []
        for i in xrange(self.col_count):
            new_rows.append(self.get_column(i))
        return ModMatrix(new_rows)

    def get_append_columns(self, vect_list):
        """vect_list is a list of vectors

        each vector is appended to the matrix as a column

        """
        newmatrix = []
        for i in xrange(self.row_count):
            newrow = self.rows[i][:]
            for vect in vect_list:
                newrow.append(vect[i])
            newmatrix.append(ModVector(newrow))
        return ModMatrix(newmatrix)

    # Methods for calculating RREF and Change Of Basis matrix
    
    def swap_rows(self, i, j):
        temp = self.rows[i]
        self.rows[i] = self.rows[j]
        self.rows[j] = temp

    def el_row_op(self, const, i, j):
        """
        This adds const * row[i] to row[j]
        """
        const = const % opts.prime
        self.rows[j] = const * self.rows[i] + self.rows[j]

    def scale_row(self, const, i):
        """
        scales row[i] by const
        """
        self.rows[i] = (const % opts.prime) * self.rows[i]

    # Methods to calculate RREF and ChangeOfBasis

    def compute_rref(self):
        """Denote the matrix self by A. This calculates a matrix P for which
        PA = R, where R is in RREF.

        P is stored as self.basis_change and R is stored as self.rref
        """
        if self.rref_flag and self.basis_change_flag:
            return
        self.rref = ModMatrix(self.rows[:])
        rows = self.rref.rows
        self.basis_change = ModMatrix.identity(self.row_count)
        curr_i = 0
        # forward pass of gaussian elimination
        # curr_i denotes the pivot row
        # i denotes the column under consideration
        for i in xrange(self.col_count):
            if curr_i >= self.row_count:
                break
            curr_row = rows[curr_i]
            # if the current rows entry is 0, swap rows so it is
            # non-zero
            if curr_row[i] == 0:
                has_swapped = False
                for j in xrange(curr_i + 1, self.row_count):
                    if rows[j][i] != 0:
                        self.rref.swap_rows(curr_i, j)
                        self.basis_change.swap_rows(curr_i, j)
                        has_swapped = True
                        break
                # if we did not swap, this column consists of
                # all zeros, so move on to next column
                if not has_swapped:
                    continue
            # We now scale the current row to have leading term 1
            constant = pow(rows[curr_i][i], opts.prime - 2, \
                           opts.prime)
            self.rref.scale_row(constant, curr_i)
            self.basis_change.scale_row(constant, curr_i)
            # We now kill off non-zero entries below our current
            # pivot
            for j in xrange(curr_i + 1, self.row_count):
                if rows[j][i] != 0 :
                    const = (-rows[j][i] % opts.prime)
                    self.rref.el_row_op(const, curr_i, j)
                    self.basis_change.el_row_op(const, curr_i, j)
            curr_i += 1
        # Forward pass complete.
        # Now backwards pass
        for j in range(1, curr_i):
            #find pivot in row curr_i - j
            i = rows[curr_i - j].get_leading_index()
            # k is the row whose entry will be killed
            for k in range(curr_i - j):
                if rows[k][i] != 0:
                    c = (-rows[k][i] % opts.prime)
                    self.rref.el_row_op(c, curr_i -j, k)
                    self.basis_change.el_row_op(c, curr_i -j, k)
        self.rref_flag = True
        self.basis_change_flag = True

    def get_rref(self):
        if not self.rref_flag:
            self.compute_rref()
        return self.rref

    def get_basis_change(self):
        if not self.basis_change_flag:
            self.compute_rref()
        return self.basis_change

    def calculate_row_weights(self):
        if self.row_weights_flag:
            return
        self.row_weights = ModVector.null(self.row_count)
        for index in xrange(self.row_count):
            if self.get_rref().rows[index].get_leading_index() >= 0:
                self.row_weights[index] = 1
        self.row_weights_flag = True

    def get_row_weights(self):
        if not self.row_weights_flag:
            self.calculate_row_weights()
        return self.row_weights

    def can_solve(self, vector):
        """
        input is a ModVector, not a ModMatrix
        
        Output is True/False
        """
        if self.row_count != vector.length:
            raise TypeError("size mismatch")
        yy = ModMatrix([vector]).get_transpose()
        yy = self.get_basis_change() * yy
        for i in xrange(self.row_count):
            if not self.get_row_weights()[i] and yy.rows[i][0]:
                return False
        return True

    def solve(self, vector):
        """Input is a ModVector of correct length

        returns False if no solution else returns a column vector (as a
        ModMatrix) xx for which self * xx = vector

        """
        out = ModMatrix([self.solve_vout(vector)])
        return out.get_transpose()

    def solve_vout(self, vector):
        """Input is a ModVector of correct length

        returns False if no solution else returns a ModVector 
        as output

        """
        rref = self.get_rref()
        if self.row_count != vector.length:
            raise TypeError("size mismatch")
        yy = ModMatrix([vector]).get_transpose()
        yy = self.get_basis_change() * yy
        xx = ModVector.null(self.col_count)
        pivot_pos = [(i, rref.rows[i].get_leading_index()) for i in xrange(self.row_count)]
        pivot_cols = [row.get_leading_index() for \
                      row in self.get_rref().rows]
        for pos in pivot_pos:
            xx.vector[pos[1]] = yy.rows[pos[0]][0]
        return xx

    #    def all_solutions(self, vector):

    def get_rank(self):
        counter = 0
        for row in self.get_rref().rows:
            if row.get_leading_index() >= 0:
                counter += 1
        return counter

    def get_kernel(self):
        """This will return a list of column vectors (ModMatrix objects) which
        form a basis for the kernel of the given matrix.

        """
        rref = self.get_rref()
        pivot_pos = [(i, rref.rows[i].get_leading_index()) for i in xrange(self.row_count)]
        pivot_cols = [pair[1] for pair in pivot_pos]
        free_cols = [j for j in xrange(self.col_count) if not (j in pivot_cols)]
        out = []
        for j in free_cols:
            soln = ModVector.null(rref.col_count)
            for basic_pos in pivot_pos:
                if basic_pos[1] < j:
                    soln[basic_pos[1]] = (-rref.rows[basic_pos[0]][j] % opts.prime)
            soln[j] = 1
            out.append(ModMatrix([soln]).get_transpose())
        return out

    def get_kernel_vout(self):
        """This will return a list of vectors (ModVector objects) which
        form a basis for the kernel of the given matrix.

        """
        rref = self.get_rref()
        pivot_pos = [(i, rref.rows[i].get_leading_index()) for i in xrange(self.row_count)]
        pivot_cols = [pair[1] for pair in pivot_pos]
        free_cols = [j for j in xrange(self.col_count) if not (j in pivot_cols)]
        out = []
        for j in free_cols:
            soln = ModVector.null(rref.col_count)
            for basic_pos in pivot_pos:
                if basic_pos[1] < j:
                    soln[basic_pos[1]] = (-rref.rows[basic_pos[0]][j] % opts.prime)
            soln[j] = 1
            out.append(soln)
        return out
        
    def get_inverse(self):
        #check if square matrix
        if self.row_count != self.col_count:
            print "not square"
            return False 
        #check if full rank
        if self.get_rank() != self.row_count:
            print "not invertible"
            return False
        #solve e_i = Ax for all i
        out_columns = []
        for i in range(self.row_count):
            ee = ModVector.null(self.row_count)
            ee[i] = 1
            xx = self.solve_vout(ee)
            out_columns.append(xx)
        #assemble x_i into matrix 
        out_matrix = ModMatrix(out_columns)
        return out_matrix.get_transpose()
