# Official packages
import os 
import unittest

# Third-party packages
import numpy as np

# In-house packages
from BrukerSpSpEpi import BrukerSpSpEpiExp


class TestBrukerSpSpEpi(unittest.TestCase):

    def setUp(self) -> None:
        self.exp_dir = "../20230305_163320_AgroseCylinder2_1_1"
        self.exp_nbr = 21
        
        self.spsp_epi = BrukerSpSpEpiExp(self.exp_dir, self.exp_nbr)

    def test_exp_path(self):
        self.assertEqual(self.spsp_epi.exp_path , os.path.join(self.exp_dir, str(self.exp_nbr)))

    def test_update_data_paths(self):
        self.assertEqual(self.spsp_epi.data_paths_dict['fid_path'] , os.path.join(self.spsp_epi.exp_path, "fid"))
        self.assertEqual(self.spsp_epi.data_paths_dict['method_path'] , os.path.join(self.spsp_epi.exp_path, "method"))
        self.assertEqual(self.spsp_epi.data_paths_dict['acqp_path'] , os.path.join(self.spsp_epi.exp_path, "acqp"))
       
    def test_read_param_dicts(self):
        self.assertEqual(self.spsp_epi.method_dict['Method'] , "<User:xji_spsp_epi>")
        self.assertEqual(self.spsp_epi.method_dict['PVM_NAverages'] , 256)
        self.assertTrue(np.array_equal(self.spsp_epi.method_dict['PVM_Matrix'], np.array([120,120]), equal_nan=True))

        PVM_EncSteps1_array = np.linspace(start=-20, stop=59, num=80, endpoint=True, dtype=float)
        self.assertTrue(np.array_equal(self.spsp_epi.method_dict['PVM_EncSteps1'], PVM_EncSteps1_array, equal_nan=True))
        self.assertAlmostEqual(self.spsp_epi.method_dict['PVM_EpiReadOddGrad'] , -0.0805918017183784)

    def test_calc_exp_data_dimensionality(self):
        self.assertAlmostEqual(self.spsp_epi._exp_data_dimensionality_dict['dim_rf_offset'] , (1,))
        self.assertAlmostEqual(self.spsp_epi._exp_data_dimensionality_dict['dim_rf_flip_angle'] , (1,))


if __name__ == '__main__':
    unittest.main()
    
        