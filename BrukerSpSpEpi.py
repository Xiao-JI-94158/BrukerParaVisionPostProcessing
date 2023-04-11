# Official packages
import os
import copy

from typing import List

# Third-party packages
import numpy as np
import pandas as pd

# In-house packages

"""
class SpSpEpiTransientSpace(object):
    def __init__(self, **kwargs):
        
        self._nbr_metabolites       = kwargs.get('nbr_metabolites', None)
        self._metabolite_offset_list= kwargs.get('metabolite_offset_list', None)
        self._nbr_time_points       = kwargs.get('nbr_time_points', None)
        self._dim_k_raw_ro          = kwargs.get('dim_k_raw_ro', None)
        self._dim_k_raw_ph          = kwargs.get('dim_k_raw_ph', None)
        self._dim_r_image_ro        = kwargs.get('dim_r_image_ro', None) 
        self._dim_r_image_ph        = kwargs.get('dim_r_image_ph', None)
        self._fov_ro                = kwargs.get('fov_ro', None)        
        self._fov_ph                = kwargs.get('fov_ph', None)        
        self.transients             = None
"""

TRANSIENT_ENTRIES = {
    "rf_offset"     : "float",
    "time_pts"      : "int",
    "k_ro_pts"      : "int",
    "k_ph_pts"      : "int",
    "k_space_pos"   : "object",
    "k_space_neg"   : "object",
    "r_ro_pts"      : "int",
    "r_ph_pts"      : "int",
    "r_ro_fov_mm"   : "float",
    "r_ph_fov_mm"   : "float",
    "r_space_pos"   : "object",
    "r_space_neg"   : "object",
    "r_space_abs"   : "object"
}

DEFAULT_METABOLITES = ['Urea','Pyruvate','Lactate']

TEST_CHEMICAL = ['Urea']

class BrukerSpSpEpiExp(object):
    """
        Basic Class that that read, stores, and (post-)processes EPI data acquired from Spectral-Spatial Selective Excitation (SpSp_EPI)
    """
    
    def __init__(self, exp_data_path:str, metabolite_list:List[str] = DEFAULT_METABOLITES ) -> None:
        """
        """
        self.metabolite_list = metabolite_list
        self.data_paths_dict = self._update_data_paths(exp_data_path)             
        self.transient_space = self._generate_transient_space()  

        self.param_dict = (self._read_param_dicts(self.data_paths_dict['method']) | self._read_param_dicts(self.data_paths_dict['acqp'])) 

        self._validate() 
        
        self.data = self.reconstruct_transient_space()





    def _update_data_paths(self, exp_data_path)->None:
        data_paths_dict = {}
        
        if (not (os.path.isdir(exp_data_path))):
            raise NotADirectoryError("Given directory of Experiment does not exist")

        fid_path =  os.path.join(exp_data_path, "fid")
        if (not (os.path.isfile(fid_path))):
            raise FileNotFoundError("Cannot find FID file in the given directory of Experiment")
        data_paths_dict['fid'] = fid_path

        method_path = os.path.join(exp_data_path, "method")
        if (not (os.path.isfile(method_path))):
            raise FileNotFoundError("Cannot find METHOD file in the given directory of Experiment")
        data_paths_dict['method'] = method_path

        acqp_path = os.path.join(exp_data_path, "acqp")
        if (not (os.path.isfile(acqp_path))):
            raise FileNotFoundError("Cannot find ACQP file in the given directory of Experiment")
        data_paths_dict['acqp'] = acqp_path
        return data_paths_dict
       

    def _read_param_dicts(self, param_file_path):
        """
        Read a Bruker MRI experiment's method or acqp file to a dictionary.

        Ref: https://github.com/jdoepfert/brukerMRI
        """

        param_dict = {}

        with open(param_file_path, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break

                # when line contains parameter
                if line.startswith('##$'):

                    (param_name, current_line) = line[3:].split('=') # split at "="

                    # if current entry (current_line) is arraysize
                    if current_line[0:2] == "( " and current_line[-3:-1] == " )":
                        value = self._parse_array(f, current_line)

                    # if current entry (current_line) is struct/list
                    elif current_line[0] == "(" and current_line[-3:-1] != " )":

                        # if neccessary read in multiple lines
                        while current_line[-2] != ")":
                            current_line = current_line[0:-1] + f.readline()

                        # parse the values to a list
                        value = [self._parse_single_value(x) for x in current_line[1:-2].split(', ')]

                    # otherwise current entry must be single string or number
                    else:
                        value = self._parse_single_value(current_line)

                    # save parsed value to dict
                    param_dict[param_name] = value

        return param_dict
        

    def _parse_array(self, current_file, line):
        """
        Ref: https://github.com/jdoepfert/brukerMRI
        """
        # extract the arraysize and convert it to numpy
        line = line[1:-2].replace(" ", "").split(",")
        arraysize = np.array([int(x) for x in line])

        # then extract the next line
        vallist = current_file.readline().split()

        # if the line was a string, then return it directly
        try:
            float(vallist[0])
        except ValueError:
            return " ".join(vallist)

        # include potentially multiple lines
        while len(vallist) != np.prod(arraysize):
            vallist = vallist + current_file.readline().split()

        # try converting to int, if error, then to float
        try:
            vallist = [int(x) for x in vallist]
        except ValueError:
            vallist = [float(x) for x in vallist]

        """
        # This block below is the original code from Ref: https://github.com/jdoepfert/brukerMRI
        # For our purpose, we return all numerical types in format of numpy.ndarray, regardless of its length

            # convert to numpy array
            if len(vallist) > 1:
                return np.reshape(np.array(vallist), arraysize)
            # or to plain number
            else:
                return vallist[0]
        """
        return np.reshape(np.array(vallist), arraysize)

    def _parse_single_value(self, val):
        """
        Ref: https://github.com/jdoepfert/brukerMRI
        """
        try: # check if int
            result = int(val)
        except ValueError:
            try: # then check if float
                result = float(val)
            except ValueError:
                # if not, should  be string. Remove  newline character.
                result = val.rstrip('\n')

        return result    

    def _generate_transient_space(self) -> dict[str, pd.DataFrame]:
        transient_space = {}
        for metabolite in self.metabolite_list:
            transient_space[metabolite] = pd.DataFrame( 
                                                {col_name: pd.Series(dtype=col_type) for col_name, col_type in TRANSIENT_ENTRIES.items()}
                                            )
        return transient_space

    def _validate(self):
        _nbr_metabolite = len(self.metabolite_list)
        _nbr_cs_offsite = self.param_dict['NumChemicalShifts']
        if ( _nbr_metabolite != _nbr_cs_offsite ):
            raise ValueError( f'Number of metabolite names ({_nbr_metabolite}) provided by user does not match the number of chemical shift offsets ({_nbr_cs_offsite}) in the raw data.' )
    
    def reconstruct_transient_space(self):
        fid = self._read_raw_fid()
        fid = self._deserialize_raw_fid(fid)
        fid = np.array_split(fid, self.param_dict['PVM_NRepetitions'])
        fid = np.reshape(fid, newshape=(  , -1))

        return fid

    def _read_raw_fid(self) -> np.ndarray:
        _raw_fid_dtype = self.param_dict['GO_raw_data_format']
        if (_raw_fid_dtype == 'GO_32BIT_SGN_INT') :
            fid = np.fromfile(file=self.data_paths_dict['fid'], dtype='int32')
        else:
            raise TypeError( f'Raw FID data in Unknown Datatype ({_raw_fid_dtype})' )
    
    def _deserialize_raw_fid(self, fid) -> np.ndarray:
        return (fid[0::2, ...] + 1j * fid[1::2, ...])

    def _generate_k_space_data(self)->dict:
        """
        When the length of FID doubles the length of raw k-space encoding,
            We assume this is <Double Sampling> type EPI readout.
            Under this assumption:
                No even-line-mirroring is performed, since the two samplings comes from Pos. and Neg lobe of EPI readout respectively.
                Raw k-space encoding is splited into two sub-k-spaces, the Pos and the Neg, for reconstruction, which has exactly opposite phases
                Echo-centers of readout lines in the two sub-k-spaces are aligned
        When the length of FID equals the length of raw k-space encoding,
            It is conventional EPI readout.
            Under this assumption:
                Even-line-mirroring is performed
                Align echo-centers of each line of readout
        """
            
        _k_space_data = {}
        _k_space_encoding_length = (self._exp_data_dim_dict['dim_k_raw_ph'] * self._exp_data_dim_dict['dim_k_raw_ro'])
        if (np.shape(self.fid['deserialized'])[0] == 2 * _k_space_encoding_length):
            _k_space_2d = np.reshape(self.fid['deserialized'], (self._exp_data_dim_dict['dim_k_raw_ph'], -1))
            _k_space_data["Raw"] = _k_space_2d
            _k_space_2d = self._zerofill_fid_2d(_k_space_2d)
            _k_space_data["Pos"], _k_space_data["Neg"] = self._split_fid_2d(_k_space_2d)
            _k_space_data["Pos"] = self._align_echo_center(_k_space_data["Pos"])
            _k_space_data["Neg"] = self._align_echo_center(_k_space_data["Neg"])
            _k_space_data["Neg"] = np.fliplr(_k_space_data["Neg"])
        elif (np.shape(self.fid['deserialized'])[0] == _k_space_encoding_length):
            _k_space_2d = np.reshape(self.fid['deserialized'], (self._exp_data_dim_dict['dim_k_raw_ph'], self._exp_data_dim_dict['dim_k_raw_ro']))
            _k_space_data["Raw"] = _k_space_2d
            _k_space_data["Pos"] = self._align_echo_center(_k_space_2d)
        else:
            raise ValueError(f"Size of Raw FID ({np.shape(self.fid['deserialized'])}) is not equal to the size of partial-phased k-space ({_k_space_encoding_length}), nor the double of it")
        
        return _k_space_data
        

    def _zerofill_fid_2d(self, fid_2d):
        nbr_zerofill_lines = self._exp_data_dim_dict['dim_r_image_ph'] - self._exp_data_dim_dict['dim_k_raw_ph']
        zerofilled_fid_2d = np.pad(fid_2d, ((nbr_zerofill_lines, 0),(0, 0)), 'constant', constant_values=(0))
        return zerofilled_fid_2d

    def _split_fid_2d(self, fid_2d):
        fid_left  = fid_2d[...,:self._exp_data_dim_dict['dim_k_raw_ro']:]
        fid_right = fid_2d[...,self._exp_data_dim_dict['dim_k_raw_ro']::]
        return fid_left, fid_right

    def _align_echo_center(self, fid_2d):
        
        for idx, line in enumerate(fid_2d):
            shift = 60 - np.argmax(np.abs(line))
            fid_2d[idx] = np.roll(line, shift)
        return fid_2d

    def _reconstruct_r_space_data(self)->dict:
        """
        When the k_space_data has both Pos and Neg part
            We assume this is <Double Sampling> type EPI readout.
            Under this assumption:
                FT for two sub-r-space image, the Pos and the Neg.
                Yield the magnitude result from the addition of both magnitude images from the Pos and the Neg.
        When the k_space_data has only Pos part
            It is conventional EPI readout.
            Under this assumption:
                FT for r-space image
                Yield the magnitude result
        """
        _r_space_data = {}
        
        _r_space_data["Pos"] = np.fft.fftshift(np.fft.fft2(self.k_space_data["Pos"]))
        if ('Neg' in self.k_space_data.keys()):
            _r_space_data["Neg"] = np.fft.fftshift(np.fft.fft2(self.k_space_data["Neg"]))

        if ('Neg' in self.k_space_data.keys()):
            _r_space_data["Mag"] = np.abs(_r_space_data["Pos"]) + np.abs(_r_space_data["Neg"])
        else:
            _r_space_data["Mag"] = np.abs(_r_space_data["Pos"])
        
        return _r_space_data

    def _calc_magnitude_proj(self)->dict:
        magnitude_proj = {}
        magnitude_proj['ro'] = np.sum(self.r_space_data['Mag'], axis=0)
        magnitude_proj['ph'] = np.sum(self.r_space_data['Mag'], axis=1)
        return magnitude_proj
    





if __name__ == "__main__":

    pass
