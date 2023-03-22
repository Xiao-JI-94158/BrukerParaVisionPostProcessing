import os 




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



 
    




if __name__ == '__main__':
    dataset_path = "../20230305_163320_AgroseCylinder2_1_1"
    exp_nbr = 21


    paths_dict = validate_data_paths(dataset_path, exp_nbr)
    
        