import os
import glob
import uuid
import iris

def _pp2nc_regrid(combined_ppfile, ncfile):
    cube = iris.load_cube(combined_ppfile)

    #regrid to 2.5 deg
    #cube = _regrid2obs(cube)

    # tropics alone
    cube = cube.intersection(latitude=(-30, 30), longitude=(0, 360))
    print(cube)
    iris.save(cube, ncfile, netcdf_format="NETCDF3_CLASSIC")

def model_data_retrieve(varnames, control=None, expt=None, netcdf=True):

    runs = [control, expt]

    for run in runs:
        runid = run['runid']
        start_date = run['start_date']
        end_date = run['end_date']

        # data directory to retrieve data
        data_root = os.path.join(run['data_retrieve_dir'], runid)
        if not os.path.exists(data_root):
            os.makedirs(data_root)

        # Temporary filter file to hold the MOOSE filters
        filter_tmp_file = os.path.join('/var/tmp/diags_filter_temp_level_%s.txt' % uuid.uuid4())

        for varname in varnames:
            filtList = glob.glob(os.path.join('data', varname + '*.filter'))
            for f in filtList:
                # Get template Filter file name
                # Get variable name from file name
                # Get MOOSE data stream from file name
                ffile = os.path.basename(f)
                varname = ffile.split(".")[0]
                stream = ffile.split(".")[1]

                # Final file name of retrieved data as a single file
                combined_ppfile = os.path.join(data_root, runid + '_' + varname + '.pp')
                ncfile = combined_ppfile + '.nc'

                if not os.path.exists(ncfile):
                    # open the filter query file
                    fx = open(f, 'r')

                    # open a fresh temporary filter query file with the time info
                    if os.path.exists(filter_tmp_file):
                        os.remove(filter_tmp_file)

                    w = open(filter_tmp_file, 'w')
                    lines = fx.readlines()
                    for line in lines:
                        if line.strip() == "begin":
                            w.write(line)
                            w.write('T1>={%s}\n' % start_date)
                            w.write('T2<={%s}\n' % end_date)
                        else:
                            w.write(line)
                    fx.close()
                    w.close()

                    print('*' * 40)
                    print(varname)
                    print(filter_tmp_file)
                    os.system('cat ' + filter_tmp_file)
                    print('*' * 40)

                    if not os.path.exists(combined_ppfile):
                        print('%s not found. Retrieving from MASS...' % combined_ppfile)
                        moodir = 'moose:/crum/%s/ap%s.pp' % (runid, stream)
                        command = 'moo select -C %s %s %s' % (filter_tmp_file, moodir, combined_ppfile)
                        print(command)
                        stat = os.system(command)
                        print('Sort pp headers...')
                        os.system('./utils/sortpp %s' % combined_ppfile)

                    if netcdf:
                        if not os.path.exists(ncfile):
                            print('Converting to netcdf files...')
                            if os.path.exists(combined_ppfile):
                                _pp2nc_regrid(combined_ppfile, ncfile)
                        if os.path.exists(ncfile):
                            os.remove(combined_ppfile)
                        else:
                            print('PP to netcdf conversion failed.')
                else:
                    print('%s file exist. Skipping retrieval...' %ncfile)
                    break

if __name__ == '__main__':
    control = {}
    control['runid'] = 'u-ab680'
    control['start_date'] = '1989/01/01'
    control['end_date'] = '2009/12/29'
    control['data_retrieve_dir'] = '/scratch/hadpx/hadgem3/data/SEAPy'

    expt = {}
    expt['runid'] = 'u-ab680'
    expt['start_date'] = '1989/01/01'
    expt['end_date'] = '2009/12/29'
    expt['data_retrieve_dir'] = '/scratch/hadpx/hadgem3/data/SEAPy'

    obs = {}
    obs['runid'] = 'obs'
    obs['start_date'] = '1989/01/01'
    obs['end_date'] = '2009/12/29'
    obs['data_retrieve_dir'] = '/project/MJO_GCSS/hadgem3/data/'

    print(control, expt)
    varnames = ['U850', 'V850', 'PRECIP', 'SST']
    model_data_retrieve(varnames, control=control, expt=expt)
