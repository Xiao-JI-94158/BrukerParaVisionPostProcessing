# Official packages
import os
import copy

from typing import List

# Third-party packages
import numpy as np
import pandas as pd

# In-house packages

RECONSTRUCTION_PARAMETERS = {
    'does_zerofilling'          = False,
    'does aligning_echo_center' = False
}

TRANSIENT_ENTRIES = {
    "time_pts"       : None,
    "raw_fids"       : None,
    "k_spaces_pos"   : None,
    "k_spaces_neg"   : None,
    "r_images_pos"   : None,
    "r_images_neg"   : None,
    "r_images_abs"   : None
}

DEFAULT_METABOLITES = ['Urea','Pyruvate','Lactate']

TEST_CHEMICAL = ['Urea']

class BrukerSpSpEpiExp(object):
    """
        Basic Class that that read, stores, and (post-)processes EPI data acquired from Spectral-Spatial Selective Excitation (SpSp_EPI)
    """
    
    def __init__(self, exp_data_path:str, metabolite_list:List[str] = DEFAULT_METABOLITES, **kwargs) -> None:
        """
        """
        self.metabolite_list = metabolite_list
        self.data_paths_dict = self._update_data_paths(exp_data_path)             
        
        self.recon_params    = _retrieve_recon_params(kwargs)

        self.param_dict      = (self._read_param_dicts(self.data_paths_dict['method']) | self._read_param_dicts(self.data_paths_dict['acqp'])) 
        self.transient_space = self._generate_transient_space()  
        self._validate() 

        
        self.data = self.reconstruct_transient_space( does_zerofilling=False, does_align_echo_center=False )


    def _retrieve_recon_params(self, kwargs):
        recon_params = copy.deepcopy(RECONSTRUCTION_PARAMETERS)
        recon_params.update((k, kwargs[k]) for k in (recon_params.keys() & kwargs.keys()) )
        return recon_params



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
            transient_space[metabolite] = copy.deepcopy(TRANSIENT_ENTRIES)
        return transient_space

    def _validate(self):
        _nbr_metabolite = len(self.metabolite_list)
        _nbr_cs_offset = self.param_dict['NumChemicalShifts']
        if ( _nbr_metabolite != _nbr_cs_offset ):
            raise ValueError( f'Number of metabolite names ({_nbr_metabolite}) provided by user does not match the number of chemical shift offsets ({_nbr_cs_offset}) in the raw data.' )
    
    def _reconstruct_transient_space(self):
        
        _nbr_metabolite = self.param_dict["NumChemicalShifts"]
        _nbr_repetition = self.param_dict["NR"]
        _nbr_time_pts = int(np.ceil( _nbr_repetition / float(_nbr_metabolite)))
        _TR_sec = self.param_dict['ACQ_repetition_time'] / 1000
        time_pts = np.tile(_TR_sec, reps=_nbr_repetition) + np.tile(self.param_dict['Vd1List'], reps=_nbr_time_pts) # time-interval list
        time_pts = np.concatenate((np.array([0]), time_pts[:-1:])) # ACQ starts at time=0!
        time_pts = np.cumsum(time_pts) # time of transient excution is cumulative sum of time-intervals between consecutive transients
                
        raw_fids = self._read_raw_fid()
        raw_fids = self._deserialize_raw_fid(raw_fids)
        raw_fids = np.array_split(raw_fids, _nbr_repetition)
        #fid = fid.reshape( _nbr_time_pts, _nbr_metabolite, -1)
        #fid = np.swapaxes(fid, axis1=0, axis2=1)
        #fid = np.squeeze(fid)

        pos_k_spaces, neg_k_spaces = zip(*map(self._construct_k_space(), raw_fids))

        def ft2d(x): return np.fft.fftshift(np.fft.fft2(x))

        pos_r_iamges = [ft2d(k_space_pos) for k_space_pos in pos_k_spaces]
        neg_r_iamges = [ft2d(k_space_neg) for k_space_neg in neg_k_spaces]

        
        return NotImplemented


    def _read_raw_fid(self) -> np.ndarray:
        _raw_fid_dtype = self.param_dict['GO_raw_data_format']
        if (_raw_fid_dtype == 'GO_32BIT_SGN_INT') :
            fid = np.fromfile(file=self.data_paths_dict['fid'], dtype='int32')

        else:
            raise TypeError( f'Raw FID data in Unknown Datatype ({_raw_fid_dtype})' )
        
        return fid
    
    def _deserialize_raw_fid(self, fid) -> np.ndarray:
        return (fid[0::2, ...] + 1j * fid[1::2, ...])

    def _construct_k_space(self, fid, does_zerofilling=False, does_align_echo_center=False) -> List[np.ndarray, np.ndarray]:
        dim_k_raw_ro, dim_k_raw_ph = self.param_dict['PVM_EncMatrix']
        dim_r_img_ro, dim_r_img_ph = self.param_dict['PVM_Matrix']
        k_space_encoding_length = (dim_k_raw_ph * dim_k_raw_ro)

        if (len(fid) == 2 * k_space_encoding_length):
            """
            When the length of FID doubles the length of raw k-space encoding,
                We assume this is <Double Sampling> type EPI readout.
                Under this assumption:
                    No even-line-mirroring is performed, since the two samplings comes from Pos. and Neg lobe of EPI readout respectively.
                    Raw k-space encoding is splited into two sub-k-spaces, the Pos and the Neg, for reconstruction, which has exactly opposite phases
                    Echo-centers of readout lines in the two sub-k-spaces are aligned
            """            
            _fid_2d = np.reshape(fid, (dim_k_raw_ph, -1))
            
            _k_space_pos, _k_space_neg = self._split_fid_2d(_fid_2d)

            if does_zerofilling:
                _k_space_pos = self._zerofill_kspace(_k_space_pos)
                _k_space_neg = self._zerofill_kspace(_k_space_neg)
            
            if does_align_echo_center:
                _k_space_pos = self._align_echo_center(_k_space_pos)

                _k_space_neg = self._align_echo_center(_k_space_neg)
                _k_space_neg = np.fliplr(_k_space_neg)

        elif (len(fid) == k_space_encoding_length):
            """
            When the length of FID equals the length of raw k-space encoding,
                It is conventional EPI readout.
                Under this assumption:
                    Even-line-mirroring is performed
                    Align echo-centers of each line of readout
            """
            _k_space_neg = np.array([])
            _k_space_pos = np.reshape(fid, (dim_k_raw_ph, dim_k_raw_ro))
            _k_space_pos[1::2, ::] = _k_space_pos[1::2, ::-1]
            if does_zerofilling:
                _k_space_pos = self._zerofill_kspace(_k_space_pos)
                _k_space_neg = self._zerofill_kspace(_k_space_neg)
            
            if does_align_echo_center:
                _k_space_pos = self._align_echo_center(_k_space_pos)
        else:
            raise ValueError(f"Size of Raw FID ({np.shape(self.fid['deserialized'])}) is not equal to the size of partial-phased k-space ({k_space_encoding_length}), nor the double of it")
        
        return [_k_space_pos, _k_space_neg]
        

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

    def _reconstruct_r_image(self)->dict:
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
