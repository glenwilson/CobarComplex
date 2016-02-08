from pickle_storage import PickleStorage
import cPickle as pickle

class Options(object):
    def __init__(self, prime, degree_bounds, numpcs, coprod_db):
        self.prime = prime
        self.bounds = degree_bounds
        self.numpcs = numpcs
        self.coprod_db = coprod_db
        self.db_flag = False
        
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


opts = Options(17, 10, 2, PickleStorage("coprod_db.pickle"))
