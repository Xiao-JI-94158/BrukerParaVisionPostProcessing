# Official packages
import os
import copy

from typing import List, Dict

# Third-party packages
import numpy as np

# In-house packages

POST_PROCESSING_PARAMETERS = {
    'does_zerofill'             : False,
    'does_align_echo_center'    : False,
    'does_reconstruction'       : False,
    'has_double_sampling'       : False
}

RAW_DATA_SET = {
    'fid'       : 'fid',
    # 'ser'       : 'ser',
    # pdata subdir might not exist due to user config (not performing factory reconstruction)
    '2dseq'     : 'pdata/1/2dseq',
}


RAW_PARAM_SET = {
    'acqp'      : 'acqp',
    'method'    : 'method',
    'visu_pars' : 'visu_pars',
    # pdata subdir might not exist due to user config (not performing factory reconstruction)
    'procs'     : 'pdata/1/procs',
    'reco'      : 'pdata/1/reco'
}

DATA_COLLECTION_TEMPLATE = {
    'time_point_sec': None,
    'fid': None,
    '2dseq': None,
    'k_space': None,
    'r_image': None
}

class BrukerSpSpEpiExp(object):
    """
        Basic Class that that read, stores, and (post-)processes EPI data acquired from Spectral-Spatial Selective Excitation (SpSp_EPI)

        Parameters:
        -----------

        exp_data_path: str
            experimental data path

    """
    
    def __init__(self, exp_dataset_path:str, **kwargs) -> None:
        """
        0. update params for post-processing:
    
        1. validate experiment dataset:
            1.1 data files:
                must        : fid
                optional    : 2dseq
            1.2 param files:
                must        : acqp, method, visu_pars
                optional    : acqu, acqus, procs, reco
        
        2. update dataset_dict['PARAM']:

        3. update data_collection:

        4. perform reconstruction:

        """
        self.post_processing_params = self._update_post_processing_params(kwargs)

        self.dataset = {"DATA": None, "PARAM": None}

        self._validate_dataset_files(exp_dataset_path)

        self._update_dataset_param()       
        
        self._update_dataset_data()
        
        if self.post_processing_params['does_reconstruction']:
            self._update_k_space()
            self._update_r_image()


    def _update_post_processing_params(self, kwargs):
        """
        parse possible tags for post-processing
        """
        recon_params = copy.deepcopy(POST_PROCESSING_PARAMETERS)
        recon_params.update((k, kwargs[k]) for k in (recon_params.keys() & kwargs.keys()) )
        return recon_params  

    def _validate_dataset_files(self, exp_dataset_path):
        """
        Confirm that the given path of experimental dataset is valid
        """
        if (not (os.path.isdir(exp_dataset_path))):
            raise OSError(f"Given directory of Experiment ({exp_dataset_path}) does not exist")
        
        self._validate_data_files(exp_dataset_path)
        self._validate_param_files(exp_dataset_path)

        

    def _validate_data_files(self, exp_dataset_path)->Dict:
        """
        1.1 data files:
            must        : fid
            optional    : 2dseq
        """
        data_dict = RAW_DATA_SET        
        self.dataset['DATA'] = self._complete_abs_path(data_dict, exp_dataset_path)                

    def _validate_param_files(self, exp_dataset_path)->Dict:
        """
        1.2 param files:
            must        : acqp, method, visu_pars
            optional    : acqu, acqus, procs, reco        
        """
        param_dict = RAW_PARAM_SET
        self.dataset['PARAM'] = self._complete_abs_path(param_dict, exp_dataset_path)

    def _complete_abs_path(self, dp_dict, exp_dataset_path):
        """
        """
        ret_dict = copy.deepcopy(dp_dict)
        for key, value in ret_dict.items():
            abs_path = os.path.join(exp_dataset_path, value)

            if ('pdata' in value):
                if (os.path.isfile(abs_path)):
                    ret_dict[key] = abs_path
                else:
                    ret_dict[key] = None
            else:
                if (os.path.isfile(abs_path)):
                    ret_dict[key] = abs_path
                else:
                    raise FileNotFoundError(f"Cannot find {key} file ({abs_path}) in the given directory of Experiment")
        return ret_dict


    def _update_dataset_param(self):
        """
        """
        param_dict = {}
        for key, value in self.dataset['PARAM'].items():
            temp_dict = self._read_param_dicts(value)
            param_dict = (param_dict | temp_dict)
        
        self.dataset['PARAM'] = param_dict
        
    
    def _update_dataset_data(self):
        data = DATA_COLLECTION_TEMPLATE

        data['time_point_sec'] = self._calc_time_pts()

        data['fid'] = self._process_fid()

        if self.dataset['DATA']['2dseq']:
            data['2dseq'] = self._process_2dseq()
        
        self.dataset['DATA'] = data
        

    def _read_param_dicts(self, param_file_path):
        """
        Read a Bruker MRI experiment's parameter files to a dictionary.

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

    def _calc_time_pts(self):
        """
        Calculate the time point of acquisition for each EPI transient.
        The first transient starts at <0.0 second>
        """
        _nbr_metabolite = self.dataset['PARAM']["NumChemicalShifts"]
        _nbr_repetition = self.dataset['PARAM']["NR"]
        _nbr_time_pts = int(np.ceil( _nbr_repetition / float(_nbr_metabolite)))
        _TR_sec = self.dataset['PARAM']['ACQ_repetition_time'] / 1000
        time_pts = np.tile(_TR_sec, reps=_nbr_repetition) + np.tile(self.dataset['PARAM']['Vd1List'], reps=_nbr_time_pts) # time-interval list
        time_pts = np.concatenate((np.array([0]), time_pts[:-1:])) # ACQ starts at time=0!
        time_pts = np.cumsum(time_pts) # time of transient excution is cumulative sum of time-intervals between consecutive transients
        return np.asarray(time_pts)

    def _process_fid(self):
        """
        Read binary fid into cmplx128 format, and partition into transients.
        """
        raw_fids = self._read_raw_fid()
        raw_fids = self._deserialize_raw_fid(raw_fids)
        raw_fids = np.asarray(np.array_split(raw_fids, self.dataset['PARAM']["NR"]))
        return raw_fids

    def _read_raw_fid(self) -> np.ndarray:
        """
        """
        _raw_fid_dtype = self.dataset['PARAM']['GO_raw_data_format']
        if (_raw_fid_dtype == 'GO_32BIT_SGN_INT') :
            fid = np.fromfile(file=self.dataset['DATA']['fid'], dtype='int32')

        else:
            raise TypeError( f'Raw FID data in Unknown Datatype ({_raw_fid_dtype})' )
        return fid
    
    def _deserialize_raw_fid(self, fid) -> np.ndarray:
        fid = np.asarray(fid[0::2, ...] + 1j * fid[1::2, ...])
        fid.astype(np.complex128)
        return fid

    def _process_2dseq(self):
        """
        """
        _raw_2dseq_dtype = self.dataset['PARAM']['VisuCoreWordType']
        _raw_2dseq_b_order = self.dataset['PARAM']['VisuCoreByteOrder']

        if ((_raw_2dseq_dtype == '_16BIT_SGN_INT') and (_raw_2dseq_b_order == 'littleEndian')):
            raw_2dseq = np.fromfile(file=self.dataset['DATA']['2dseq'], dtype='int16')
        
        data_shape = np.append(self.dataset['PARAM']['NR'], self.dataset['PARAM']["PVM_Matrix"])

        raw_2dseq = np.reshape(raw_2dseq, data_shape)
        
        return raw_2dseq
    
    def _update_k_space(self):
        """
        """
        k_space_pos = []
        k_space_neg = []
        
        dim_k_raw_ro, dim_k_raw_ph = self.dataset['PARAM']['PVM_EncMatrix']
        dim_r_img_ro, dim_r_img_ph = self.dataset['PARAM']['PVM_Matrix']
        k_space_encoding_length = (dim_k_raw_ph * dim_k_raw_ro)

        for transient in self.dataset['DATA']['fid']:
            if (len(transient) == 2 * k_space_encoding_length):
                """
                When the length of FID doubles the length of raw k-space encoding,
                    We assume this is <Double Sampling> type EPI readout.
                    Under this assumption:
                        No even-line-mirroring is performed, since the two samplings comes from Pos. and Neg lobe of EPI readout respectively.
                        Raw k-space encoding is splited into two sub-k-spaces, the Pos and the Neg, for reconstruction, which has exactly opposite phases
                        Echo-centers of readout lines in the two sub-k-spaces are aligned
                """   

                self.post_processing_params['has_double_sampling'] = True

                _transient_2d = np.reshape(transient, (dim_k_raw_ph, -1))
                
                _k_transient_pos, _k_transient_neg = self._split_fid_2d(_transient_2d, dim_k_raw_ro)

                if self.post_processing_params['does_zerofill']:
                    _k_transient_pos = self._zerofill_k_transient(_k_transient_pos)
                    _k_transient_neg = self._zerofill_k_transient(_k_transient_neg)
                
                if self.post_processing_params['does_align_echo_center']:
                    _k_transient_pos = self._align_echo_center(_k_transient_pos)
                    _k_transient_neg = self._align_echo_center(_k_transient_neg)
                
                _k_transient_neg = np.fliplr(_k_transient_neg)

                k_space_pos.append( _k_transient_pos)
                k_space_neg.append( _k_transient_neg)

            elif (len(transient) == k_space_encoding_length):
                """
                When the length of FID equals the length of raw k-space encoding,
                    It is conventional EPI readout.
                    Under this assumption:
                        Even-line-mirroring is performed
                        Align echo-centers of each line of readout
                """
                self.post_processing_params['has_double_sampling'] == False

                
                _k_transient_pos = np.reshape(transient, (dim_k_raw_ph, dim_k_raw_ro))
                _k_transient_pos[1::2, ::] = _k_transient_pos[1::2, ::-1] # Even-line-mirroring

                if self.post_processing_params['does_zerofill']:
                    _k_transient_pos = self._zerofill_k_transient(_k_transient_pos)
                    
                if self.post_processing_params['does_align_echo_center']:
                    _k_transient_pos = self._align_echo_center(_k_transient_pos)

                k_space_pos.append( _k_transient_pos)
                

            else:
                self.post_processing_params['has_double_sampling'] = None
                raise ValueError(f"Size of Raw FID ({len(transient)}) is not equal to the size of partial-phased k-space ({k_space_encoding_length}), nor the double of it")


        self.dataset['DATA']['k_space'] = {'Pos': np.asarray(k_space_pos), 'Neg':np.asarray(k_space_neg)}

    def _split_fid_2d(self, fid_2d, dim_k_raw_ro):
        """
        """
        fid_left  = fid_2d[..., :dim_k_raw_ro:]
        fid_right = fid_2d[..., dim_k_raw_ro::]
        return fid_left, fid_right
    
    def _zerofill_k_transient(self, fid_2d):
        """
        """
        dim_k_raw_ro, dim_k_raw_ph = self.dataset['PARAM']['PVM_EncMatrix']
        dim_r_img_ro, dim_r_img_ph = self.dataset['PARAM']['PVM_Matrix']

        nbr_zerofill_ph = dim_r_img_ph - dim_k_raw_ph
        nbr_zerofill_ro = dim_r_img_ro - dim_k_raw_ro

        zerofilled_fid_2d = np.pad(fid_2d, ((nbr_zerofill_ph, nbr_zerofill_ro),(0, 0)), 'constant', constant_values=(0))
        return zerofilled_fid_2d

    def _align_echo_center(self, fid_2d):
        """
        """
        for idx, line in enumerate(fid_2d):
            shift = 60 - np.argmax(np.abs(line))
            fid_2d[idx] = np.roll(line, shift)
        return fid_2d


    def _update_r_image(self):
        """
        """
        
        def ft2d(x): return np.fft.fftshift(np.fft.fft2(x))

        r_image_pos = np.asarray([ft2d(k_transient) for k_transient in self.dataset['DATA']['k_space']['Pos']])

        r_image_abs = np.abs(r_image_pos)

        if self.dataset['DATA']['k_space']['Neg'].any():
            r_image_neg = np.asarray([ft2d(k_transient) for k_transient in self.dataset['DATA']['k_space']['Neg']])
            r_image_abs += np.abs(r_image_neg)

        



        self.dataset['DATA']['r_image'] = {'Pos': r_image_pos, 'Neg': r_image_neg, "Abs": r_image_abs}
    





if __name__ == "__main__":

    pass
