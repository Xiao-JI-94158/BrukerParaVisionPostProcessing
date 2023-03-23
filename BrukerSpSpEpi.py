# Official packages
import os

# Third-party packages
import numpy as np

#In-house packages




class BrukerSpSpEpi(object):
    """
        Basic Class that that read, stores, and (post)-processes EPI data acquired from Spectral-Spatial Selective Excitation (SpSp_EPI)
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
        Read a Bruker MRI experiment's method or acqp file to a
        dictionary.
        """
        param_dict = {}

        with open(filepath, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break

                # when line contains parameter
                if line.startswith('##$'):

                    (param_name, current_line) = line[3:].split('=') # split at "="

                    # if current entry (current_line) is arraysize
                    if current_line[0:2] == "( " and current_line[-3:-1] == " )":
                        value = _parse_array(f, current_line)

                    # if current entry (current_line) is struct/list
                    elif current_line[0] == "(" and current_line[-3:-1] != " )":

                        # if neccessary read in multiple lines
                        while current_line[-2] != ")":
                            current_line = current_line[0:-1] + f.readline()

                        # parse the values to a list
                        value = [_parse_single_value(x) for x in current_line[1:-2].split(', ')]

                    # otherwise current entry must be single string or number
                    else:
                        value = _parse_single_value(current_line)

                    # save parsed value to dict
                    param_dict[param_name] = value

        return param_dict
        

    def _parse_array(current_file, line):

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

        # convert to numpy array
        if len(vallist) > 1:
            return np.reshape(np.array(vallist), arraysize)
        # or to plain number
        else:
            return vallist[0]

    def _parse_single_value(val):

        try: # check if int
            result = int(val)
        except ValueError:
            try: # then check if float
                result = float(val)
            except ValueError:
                # if not, should  be string. Remove  newline character.
                result = val.rstrip('\n')

        return result    

    
    def __init__(self, exp_path) -> None:
        """
            Constructor of the class
        """

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

        



#    @property

if __name__ == "__main__":

    pass
