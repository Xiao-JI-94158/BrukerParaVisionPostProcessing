import os 

import numpy as np


def validate_data_paths(dataset_path, exp_nbr):
    paths_dict = {}
    exp_path = os.path.join(dataset_path, str(exp_nbr))
    if (not (os.path.isdir(exp_path))):
       raise NotADirectoryError
    

    fid_path =  os.path.join(exp_path, "fid")
    if (not (os.path.isfile(fid_path))):
        raise FileNotFoundError
    paths_dict['fid_path'] = fid_path

    method_path = os.path.join(exp_path, "method")
    if (not (os.path.isfile(method_path))):
        raise FileNotFoundError
    paths_dict['method_path'] = method_path

    acqp_path = os.path.join(exp_path, "acqp")
    if (not (os.path.isfile(acqp_path))):
        raise FileNotFoundError
    paths_dict['acqp_path'] = acqp_path

    return paths_dict


def read_param_file(filepath):
    """
    Read a Bruker MRI experiment's method or acqp file to a dictionary.

    Ref: https://github.com/jdoepfert/brukerMRI
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
                    value = parse_array(f, current_line)

                # if current entry (current_line) is struct/list
                elif current_line[0] == "(" and current_line[-3:-1] != " )":

                    # if neccessary read in multiple lines
                    while current_line[-2] != ")":
                        current_line = current_line[0:-1] + f.readline()

                    # parse the values to a list
                    value = [parse_single_value(x)
                             for x in current_line[1:-2].split(', ')]

                # otherwise current entry must be single string or number
                else:
                    value = parse_single_value(current_line)

                # save parsed value to dict
                param_dict[param_name] = value

    return param_dict


def parse_array(current_file, line):
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

    # convert to numpy array
    if len(vallist) > 1:
        return np.reshape(np.array(vallist), arraysize)
    # or to plain number
    else:
        return vallist[0]

def parse_single_value(val):
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

 
    




if __name__ == '__main__':
    dataset_path = "../20230305_163320_AgroseCylinder2_1_1"
    exp_nbr = 21


    paths_dict = validate_data_paths(dataset_path, exp_nbr)
    method_dict = read_param_file(paths_dict['method_path'])
    print(method_dict)
    
        