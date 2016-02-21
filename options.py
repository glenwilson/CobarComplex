from pickle_storage import PickleStorage
import cPickle as pickle

class Options(object):
    def __init__(self, prime, degree_bounds, ind, numpcs, coprod_db):
        """prime is an odd prime number

        bounds is a natural number, larger than 5 or so

        ind is a positive integer which divides (prime - 1)

        This program as written only calculates Ext for
        mot. st. alg. over finite fields. We assume we are in the case
        where the Bockstein acts non-trivially on coefficients.

        """
        self.prime = prime
        self.bounds = degree_bounds
        self.numpcs = numpcs
        self.coprod_db = coprod_db
        self.ind = ind
        self.db_flag = False
        if (self.prime - 1) % self.ind or self.ind > (self.prime - 1):
            raise TypeError("Your options are invalid!")
        
    def get_prime(self):
        return self.prime

    def get_bounds(self):
        return self.bounds

    def get_numpcs(self):
        return self.numpcs

    def get_coprod_db(self):
        if not self.db_flag:
            self.coprod_db.load_dictionary()
            self.db_flag = True
        return self.coprod_db

    def close_coprod_db(self):
        self.get_coprod_db().shutdown()
        self.db_flag = False


opts = Options(3, 10, 1, 2, PickleStorage("coprod_db.pickle"))
