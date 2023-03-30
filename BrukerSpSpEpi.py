# Official packages
import os
import copy

# Third-party packages
import numpy as np

# In-house packages




class BrukerSpSpEpiExp(object):
    """
        Basic Class that that read, stores, and (post)-processes EPI data acquired from Spectral-Spatial Selective Excitation (SpSp_EPI)
    """
    
    def __init__(self, dataset_path:str, exp_nbr:int) -> None:
        """
        """
        
        self.exp_path = os.path.join(dataset_path, str(exp_nbr))
        
        self.data_paths_dict = self._update_data_paths()
        
        
        self.method_dict = self._read_param_dicts(self.data_paths_dict['method_path'])
        self.acqp_dict = self._read_param_dicts(self.data_paths_dict['acqp_path'])
        
        self._exp_data_dim_dict = self._extract_exp_data_dims()
        
        self.fid = {}
        self.fid['raw'] = self._read_raw_fid()
        self.fid['deserialized'] = self._deserialize_raw_fid()

        self.k_space_data = self._reconstruct_k_space_data()
        self.r_space_data = self._reconstruct_r_space_data()

        
    
    def _update_data_paths(self)->None:
        data_paths_dict = {}
        
        if (not (os.path.isdir(self.exp_path))):
            raise NotADirectoryError("Given directory of Experiment does not exist")

        fid_path =  os.path.join(self.exp_path, "fid")
        if (not (os.path.isfile(fid_path))):
            raise FileNotFoundError("Cannot find FID file in the given directory of Experiment")
        data_paths_dict['fid_path'] = fid_path

        method_path = os.path.join(self.exp_path, "method")
        if (not (os.path.isfile(method_path))):
            raise FileNotFoundError("Cannot find METHOD file in the given directory of Experiment")
        data_paths_dict['method_path'] = method_path

        acqp_path = os.path.join(self.exp_path, "acqp")
        if (not (os.path.isfile(acqp_path))):
            raise FileNotFoundError("Cannot find ACQP file in the given directory of Experiment")
        data_paths_dict['acqp_path'] = acqp_path
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

    def _extract_exp_data_dims(self)->dict:
        dim_dict = {}
        
        dim_dict['dim_rf_flip_angle'] = self.method_dict['NumVarPowerExpts']
        dim_dict['dim_rf_offset'] = np.shape(self.method_dict['CSOffsetList'])[0]
        dim_dict['nbr_total_images'] = self.method_dict['PVM_NRepetitions']
        dim_dict['dim_k_raw_ro'] = self.method_dict['PVM_EncMatrix'][0]
        dim_dict['dim_k_raw_ph'] = self.method_dict['PVM_EncMatrix'][1]
        dim_dict['dim_r_image_ro'] = self.method_dict['PVM_Matrix'][0]
        dim_dict['dim_r_image_ph'] = self.method_dict['PVM_Matrix'][1]

        return dim_dict


    def _read_raw_fid(self)->np.ndarray:
        """
        """
        # raise NotImplementedError
        if (self.acqp_dict['GO_raw_data_format'] == 'GO_32BIT_SGN_INT'):
            binary_fid = np.fromfile(file=self.data_paths_dict['fid_path'], dtype='int32')
        else:
            raise TypeError('Only 32bit signed interger is accepted')
        
        return binary_fid
    
    def _deserialize_raw_fid(self)->np.ndarray:
        return (self.fid['raw'][0::2, ...] + 1j * self.fid['raw'][1::2, ...])

    def _reconstruct_k_space_data(self)->dict:
        """
        """
        # raise NotImplementedError
        _k_space_data = {}
        _k_space_encoding_length = (self._exp_data_dim_dict['dim_k_raw_ph'] * self._exp_data_dim_dict['dim_k_raw_ro'])
        if (np.shape(self.fid['deserialized'])[0] == 2 * _k_space_encoding_length):
            _k_space_2d = np.reshape(self.fid['deserialized'], (self._exp_data_dim_dict['dim_k_raw_ph'], -1))
            _k_space_2d = self._zerofill_fid_2d(_k_space_2d)
            _k_space_data["deg0"], _k_space_data["deg90"] = self._split_fid_2d(_k_space_2d)
            _k_space_data["deg0"] = self._align_echo_center(_k_space_data["deg0"])
            _k_space_data["deg90"] = self._align_echo_center(_k_space_data["deg90"])
        elif (np.shape(self.fid['deserialized'])[0] == _k_space_encoding_length):
            _k_space_2d = np.reshape(self.fid['deserialized'], (self._exp_data_dim_dict['dim_k_raw_ph'], self._exp_data_dim_dict['dim_k_raw_ro']))
            _k_space_data["deg0"] = self._align_echo_center(_k_space_2d)
        else:
            raise ValueError("Size of Raw FID is not equal to the size of partial-phased k-space, nor the double of it")
        
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
        aligned_fid_2d = copy.deepcopy(fid_2d)
        for idx, line in enumerate(fid_2d):
            shift = 60 - np.argmax(np.abs(line))
            aligned_fid_2d[idx] = np.roll(line, shift)
        return aligned_fid_2d

    def _reconstruct_r_space_data(self)->dict:
        """
        """
        _r_space_data = {}
        
        _r_space_data["deg0"]  = np.fft.fftshift(np.fft.fft2(self.k_space_data["deg0"]))
        _r_space_data["deg90"] = np.fft.fftshift(np.fft.fft2(self.k_space_data["deg90"]))
        
        return _r_space_data

        




if __name__ == "__main__":

    pass
