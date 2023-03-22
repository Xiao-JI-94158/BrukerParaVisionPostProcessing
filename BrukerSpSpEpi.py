# Official packages
import os

# Third-party packages
import numpy as np

#In-house packages




class BrukerSpSpEpi(object):
    """
        The Class that read, stores, and (post)-processes Bruker SpSpEPI data
    """
    
    def __init__(self, dataset_path:str, exp_nbr:int) -> None:
        """
        """
        
        self.dataset_path = dataset_path
        self.exp_nbr = exp_nbr
        
        self.data_paths_dict = {}
        self._update_data_paths()
        
        self.method_dict = {}
        self.acqp_dict = {}
        
        self._data_dimensionality_dict={}
        
        self.raw_fid = np.array([])
        self.k_space_data = np.array([])
        self.r_space_data = np.array([])

        
        self._update_param_dicts()
        self._read_raw_fid()
        self._reconstruct_k_space_data()
        self._reconstruct_r_space_data()
        


    
    def _update_data_paths(self)->None:
        
        exp_path = os.path.join(self.dataset_path, str(self.exp_nbr)
        if (not (os.isdir(exp_path))):
            raise NotADirectoryError("Given directory of Experiment does not exist")

        fid_path =  os.path.join(exp_path, "fid")
        if (not (os.path.isfile(fid_path))):
            raise FileNotFoundError("Cannot find FID file in the given directory of Experiment")
        self.data_paths_dict['fid_path'] = fid_path

        method_path = os.path.join(exp_path, "method")
        if (not (os.path.isfile(method_path))):
            raise FileNotFoundError("Cannot find METHOD file in the given directory of Experiment")
        self.data_paths_dict['method_path'] = mehtod_path

        acqp_path = os.path.join(exp_path, "acqp")
        if (not (os.path.isfile(acqp_path))):
            raise FileNotFoundError("Cannot find ACQP file in the given directory of Experiment")
        self.data_paths_dict['acqp_path'] = acqp_path

       

    def _update_param_dicts(self)->None:
        """
        
        """
        raise NotImplementedError
        Basic Class that performs post-processing on EPI data acquired from Spectral-Spatial Selective Excitation (SpSp_EPI)
    """

    
    def __init__(self, exp_path) -> None:
        """
            Constructor of the class

    def _read_raw_fid(self)->None:
        """
        """
        raise NotImplementedError
    
    def _reconstruct_k_space_data(self)->None:
        """
        """
        raise NotImplementedError

    def _reconstruct_r_space_data(self)->None:
        """
        """
        raise NotImplementedError
            param : exp_path, the path of dataset of the desired experiment
            return : None
            raise : 
        """
        raise NotImplementedError
        raise FileNotFoundError
        raise IsADirectoryError
        raise 
        pass




#    @property

