def dirs(dirname):
    self = dict()
    # Where the retrieved data will be stored under a runid name
    #self['data_retrieve_dir'] = '/project/MJO_GCSS/hadgem3/data/'
    #self['data_retrieve_dir'] = '/scratch/hadpx/hadgem3/data/'
    # Where will the plots and output data files go?
    self['data_out_dir'] = '/project/MJO_GCSS/hadgem3/data/MJOPy_output'
    return self[dirname]
