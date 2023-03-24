# Official packages
import os 
import unittest

# Third-party packages
import numpy as np

# In-house packages
from BrukerSpSpEpi import BrukerSpSpEpiExp


class TestBrukerSpSpEpi(unittest.TestCase):

    def test_exp_path(self):
        exp_dir = "../20230305_163320_AgroseCylinder2_1_1"
        exp_nbr = 21
        
        exp_21 = BrukerSpSpEpiExp(exp_dir, exp_nbr)

        self.assertEqual(exp_21.exp_path , os.path.join(exp_dir, str(exp_nbr)))

    def test_update_data_paths(self):
        exp_dir = "../20230305_163320_AgroseCylinder2_1_1"
        exp_nbr = 21
        
        exp_21 = BrukerSpSpEpiExp(exp_dir, exp_nbr)

        self.assertEqual(exp_21.data_paths_dict['fid_path'] , os.path.join(exp_dir, str(exp_nbr), "fid"))
        self.assertEqual(exp_21.data_paths_dict['method_path'] , os.path.join(exp_dir, str(exp_nbr), "method"))
        self.assertEqual(exp_21.data_paths_dict['acqp_path'] , os.path.join(exp_dir, str(exp_nbr), "acqp"))
       

if __name__ == '__main__':
    unittest.main()
    
        