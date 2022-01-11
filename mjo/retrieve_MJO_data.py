#!/usr/local/sci/bin/python2.7
import sys, os
import matplotlib
import glob
import uuid
import iris
import data_paths

# use the Agg environment to generate an image rather than outputting to screen
matplotlib.use('Agg')

def _regrid2obs(my_cubes):
    base_cube = iris.load_cube('obs_144x73_grid.nc')[0]
    
    #print base_cube.regrid
    base_cube.coord('latitude').coord_system = None
    base_cube.coord('longitude').coord_system = None
    my_cubes.coord('latitude').coord_system = None
    my_cubes.coord('longitude').coord_system = None
    
    # For lat/lon regriding, make sure coordinates have bounds
    if my_cubes.coord('longitude').bounds is None:
        my_cubes.coord('longitude').guess_bounds()
    if my_cubes.coord('latitude').bounds is None:
        my_cubes.coord('latitude').guess_bounds()
    if base_cube.coord('longitude').bounds is None:
        base_cube.coord('longitude').guess_bounds()
    if base_cube.coord('latitude').bounds is None:
        base_cube.coord('latitude').guess_bounds()
    reg_cube = my_cubes.regrid(base_cube, iris.analysis.Linear())
    return reg_cube

def _pp2nc_regrid(combined_ppfile, ncfile):
    cube = iris.load_cube(combined_ppfile)
    
    #regrid to 2.5 deg
    cube = _regrid2obs(cube)
    
    # tropics alone
    #cube = cube.extract(region([0, -30, 360, 30]))
    print(cube)
    iris.save(cube, ncfile)

def model_data_retrieve(control=None, expt=None):

    runs = [control, expt]

    netcdf = True
    for run in runs:
        runid = run['runid']
        start_date = run['start_date']
        end_date = run['end_date']

        # data
        data_root = os.path.join(run['data_retrieve_dir'], runid)

        # Temporary filter file to hold the filters
        filter_tmp_file = os.path.join('/var/tmp/diags_filter_temp_level_%s.txt' % uuid.uuid4())
        
        if not os.path.exists(data_root):
            os.makedirs(data_root)
        
        varnames = ['OLR', 'U850', 'U200', 'PRECIP', 'SST']

        for varname in varnames:
            filtList = glob.glob(os.path.join('data', varname + '*.filter'))
            for f in filtList:
                fFile = os.path.basename(f)
                varname = fFile.split(".")[0]
                stream = fFile.split(".")[1]
                combined_ppfile = os.path.join(data_root, runid + '_' + varname + '.pp')
                combined_ncfile = os.path.join(data_root, runid + '_' + varname + '.pp.nc')
                
                fx = open(f, 'r')             # open the filter query file
                if os.path.exists(filter_tmp_file):
                    os.remove(filter_tmp_file)
                w = open(filter_tmp_file, 'w')    # open a temporary filter query file with the time info
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

                if not os.path.exists(combined_ncfile):
                    if not (os.path.exists(combined_ppfile)):
                        print('%s not found. Retrieving from MASS...' % combined_ppfile)
                        moodir = 'moose:/crum/%s/ap%s.pp' % (runid, stream)
                        command = 'moo select -C %s %s %s' % (filter_tmp_file, moodir, combined_ppfile)
                        print(command)
                        stat = os.system(command)
                        print('Sort pp headers...')
                        os.system('./utils/sortpp %s' % combined_ppfile)

                    if netcdf:
                        if not os.path.exists(combined_ncfile):
                            print('Converting to netcdf files...')
                            if os.path.exists(combined_ppfile):
                                _pp2nc_regrid(combined_ppfile, combined_ncfile)
                        if os.path.exists(combined_ncfile):
                            os.remove(combined_ppfile)
                        else:
                            print('PP to netcdf conversion failed.')

if __name__ == '__main__':
    control = {}
    control['runid'] = 'u-ab680'
    control['start_year'] = 1989
    control['end_year'] = 2009

    expt = {}
    expt['runid'] = 'u-ab680'
    expt['start_year'] = 1989
    expt['end_year'] = 2009

    print(control, expt)
    model_data_retrieve(control=control, expt=expt)
