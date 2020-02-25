#!/usr/bin/env python
#########################################
#
#  Plotting Hart's phase space with output in
#	Hodges's TRACK program's format
#
#  JHPS
#
# C:  06/07/2013
# LU: 07/01/2016
#########################################
# 
# Usage:
#
#   - To run with default settings, run
#	$ python plot_ps.py -cl DATAFILE
# 
#   - Otherwise to select optional
#      settings on the fly, run:
#	 $ python plot_ps.py
#
#  Default columns set here:
#
#    *** HAVE MADE COUNTING START FROM 1 ****
#
vtl_col = 5  # VTL column
vtu_col = 6  # VTU column
b_col = 7    # B column
lat_col = 3   # Latitude column
inten_col = 4 # Intensity column
#
#
#########################################

##########
# Imports
##########

import pylab as plt
import sys

#########################
# Initial Configuration
#########################

cml_ticker = 0

try:
  stp = sys.argv[1]
  if stp == '-cl':
    f = open(sys.argv[2],'r')
    plot_start = 0
    plot_end = -1
  else:
    print('***ERROR***:COMMAND LINE INPUT NOT UNDERSTOOD.\n\n IF YOU WOULD LIKE COMMAND LINE SETUP USE -cl +DATAFILE (DEFAULT SETUP USED)')
    quit()

except:
  cml_ticker = True
  plot_infil = raw_input('Enter datafile to plot:\n\n')
  try:
    f = open(plot_infil,'r') 
  except:
    print('Could not find requested file')
    quit()
        
############ Get Data Columns ################

header = True

for line in f:
  e = line.replace('\n','').replace(' ','S').replace('&','').split('S')
  e = [elem for elem in e if elem != ''] 
  if header == True:
    if e[0] == 'TRACK_NUM':
      header = False
              
  elif header == False:
    if e[0] == 'TRACK_ID':
      pass      
    elif not e[0] == 'POINT_NUM':
      num_cols = len(e)
      break
f.close()

if cml_ticker == True:
  f = open(plot_infil,'r') 
else:
  f = open(sys.argv[2],'r')

##### REMEMBER have made the counting starts from one, and negative numbers count back from the end of the list.

if cml_ticker == True:
  vtl_col = raw_input('\nWhich column (of the %d available) is VTL? (Default: %d):\n\n' %(num_cols,vtl_col))
  try:
    vtl_col = int(vtl_col) - 1
  except:
    print('***ERROR***: PLEASE ENTER INPUT IN REQUESTED FORMAT')
    quit() 

if cml_ticker == True:
  vtu_col = raw_input('\nWhich column (of the %d available) is VTU? (Default: %d):\n\n' %(num_cols,vtu_col))
  try:
    vtu_col = int(vtu_col) - 1
  except:
    print('***ERROR***: PLEASE ENTER INPUT IN REQUESTED FORMAT')
    quit() 

if cml_ticker == True:
  b_col = raw_input('\nWhich column (of the %d available) is b? (Default: %d):\n\n' %(num_cols,b_col))
  try:
    b_col = int(b_col) - 1
  except:
    print('***ERROR***: PLEASE ENTER INPUT IN REQUESTED FORMAT')
    quit() 
    
if cml_ticker == True:
  lat_col = raw_input('\nWhich column (of the %d available) is latitude? (Default: %d):\n\n' %(num_cols,lat_col))
  try:
    lat_col = int(lat_col) - 1
  except:
    print('***ERROR***: PLEASE ENTER INPUT IN REQUESTED FORMAT')
    quit()
    
if cml_ticker == True:
  inten_col = raw_input('\nWhich column (of the %d available) is intensity? (Default: %d):\n\n' %(num_cols,inten_col))
  try:
    inten_col = int(inten_col) - 1
  except:
    print('***ERROR***: PLEASE ENTER INPUT IN REQUESTED FORMAT')
    quit() 

#####################
# Reading in Data
#####################

header = True

all_tracks_list = []
track_data = []
track_ids = []

for line in f:
  e = line.replace('\n','').replace(' ','S').replace('&','').split('S')
  e = [elem for elem in e if elem != ''] 

  if header == True:                # Skipping header by waiting for 'TRACK_NUM' to appear
    if e[0] == 'TRACK_NUM':
      header = False
        
        
  elif header == False:
    if e[0] == 'TRACK_ID':
      track_ids.append(e[-1])
      track_name = 'track_' + str(int(track_ids[-1]))
 
      ### Pushing points out for sorting
      all_tracks_list.append(track_data)
      track_data = []     
 
    elif not e[0] == 'POINT_NUM':
      track_data.append(e)

all_tracks_list.append(track_data)   # Need to push the last one too, which isn't signaled by a header
all_tracks_list = all_tracks_list[1:] # remove empty place holder

# Sorting data into dict
all_tracks_data = {}

for tck in range(len(track_ids)):
  track_name = 'track_' + str(track_ids[tck])
  track_data = all_tracks_list[tck]
  ### sorting for hart parameters
  lats = []
  vtl = []
  vtu = []
  b = []
  intensities = []
  for pnt in track_data:
    lats.append(float(pnt[lat_col]))
    vtl.append(float(pnt[vtl_col]))
    vtu.append(float(pnt[vtu_col]))
    b.append(float(pnt[b_col]))
    intensities.append(float(pnt[inten_col]))
      
    new_track = {
    		    'lats' : lats,
                    'vtl' : vtl,
                    'vtu' : vtu,
                    'b'   : b,
                    'intensities' : intensities
                        }

    all_tracks_data[track_name] = new_track

######################
# Final configuration
######################


###### Get which track to plot

track_chosen = 1
if cml_ticker == True:
  track_chosen = raw_input('\nWhich track (out of the %d) would you like to plot? (Default: 1):\n\n' %len(all_tracks_data))
  try:
    track_chosen = int(track_chosen)
  except:
    print('***ERROR***: PLEASE ENTER INPUT IN REQUESTED FORMAT')
    quit()

track_chosen = 'track_' + str(track_chosen)


## Get num points
plot_start = 0
plot_end = -1
if cml_ticker == True:
  plot_check = raw_input('\nWould you like to plot all points in the track? (Y/N):\n\n')
  plot_start = 0
  plot_end = -1
  if not plot_check == 'Y':
    if plot_check == 'N':
      plot_start = int(raw_input('\nWhat is the first trackpoint to plot? (1-%d):\n\n'  %len(all_tracks_data[track_chosen]['b']))) - 1
      plot_end   = int(raw_input('\nWhat is the last trackpoint to plot? (1-%d):\n\n'   %len(all_tracks_data[track_chosen]['b']))) - 1
    else:
      print('***ERROR***: PLEASE ENTER INPUT IN REQUESTED FORMAT')
      quit()
      

## Saving
plots_save = 'Y'
if cml_ticker == True:
  plots_save = raw_input('\nWould you like to save the output as postscript files? (Y/N):\n\n')
  if not plots_save in ['Y','N']:
    print('***ERROR***: PLEASE ENTER INPUT IN REQUESTED FORMAT')
    quit()      
      
      

# Setting up the plotting data

intensities = all_tracks_data[track_chosen]['intensities'][plot_start:plot_end]
B = all_tracks_data[track_chosen]['b'][plot_start:plot_end]
vtl = all_tracks_data[track_chosen]['vtl'][plot_start:plot_end]
vtu = all_tracks_data[track_chosen]['vtu'][plot_start:plot_end]
lats = all_tracks_data[track_chosen]['lats'][plot_start:plot_end]


############
# Plotting
############

# VTL vs B

f1 = plt.figure(1)

xs = []					# x = 0 line
xs = xs + [0]*(250-len(xs))

ys = []
ys = ys + [0]*(900-len(ys))		# y = 0 line


plt.plot(vtl,B,'k')
plt.scatter(vtl,B,c=intensities,s=60)

plt.plot(xs,range(-50,200,1),'k-')
plt.plot(range(-600,300,1),ys,'k-')
plt.ylabel('B [925-500hPa Storm-Relative Thickness Symmetry]')
plt.ylim(-25,125)
plt.xlabel('-VTL [925-700hPa Thermal Wind]')
plt.xlim(-600,300)

plt.annotate('A',xy=(vtl[0],B[0]),xytext=(-5,5),ha='right',textcoords='offset points')                   # labeling start and end points
plt.annotate('Z',xy=(vtl[-1],B[-1]),xytext=(-5,5),ha='right',textcoords='offset points')

cb =plt.colorbar()
cb.set_label('Intensity')

plt.text(-500,110,'ASYMMETRIC COLD-CORE')
plt.text(100,110,'ASYMMETRIC WARM-CORE')
plt.text(-500,-10,'SYMMETRIC COLD-CORE')
plt.text(100,-10,'SYMMETRIC WARM-CORE')


if plots_save == 'Y':
  saveas = track_chosen + '.VTL_B.ps'
  plt.savefig(saveas, format='ps')


# VTL vs VTU

f2 = plt.figure(2)

xs = []				          # x = 0 line
xs = xs + [0]*(900-len(xs))

ys = []					  # y = 0 line
ys = ys + [0]*(900-len(ys))


plt.plot(vtl,vtu,'k')
plt.scatter(vtl,vtu,c=intensities,s=60)

plt.plot(xs,range(-600,300,1),'k-')
plt.plot(range(-600,300,1),ys,'k-')
plt.ylabel('-VTU [700-400hPa Thermal Wind]')
plt.ylim(-600,300)
plt.xlabel('-VTL [925-700hPa Thermal Wind]')
plt.xlim(-600,300)

cb = plt.colorbar()
cb.set_label('Intensity')


plt.annotate('A',xy=(vtl[0],vtu[0]),xytext=(-5,5),ha='right',textcoords='offset points')		   # labeling start and end points
plt.annotate('Z',xy=(vtl[-1],vtu[-1]),xytext=(-5,5),ha='right',textcoords='offset points')

plt.text(-500,220,'SHALLOW COLD-CORE')
plt.text(100,220,'DEEP WARM-CORE')
plt.text(-500,-400,'DEEP COLD-CORE')
plt.text(100,-400,'SHALLOW WARM-CORE')


if plots_save == 'Y':
  saveas = track_chosen + '.VTL_VTU.ps'
  plt.savefig(saveas, format='ps')


# B vs latitude

f3 = plt.figure(3)
				
ys = [10]*(len(lats))			# Transition threshold line

plt.scatter(lats,B,c=intensities,s=20)
plt.plot(lats,ys,'k')
plt.ylabel('B [925-500hPa Storm-Relative Thickness Symmetry]')
plt.xlabel('Latitude')
plt.xlim(lats[0],lats[-1])

cb = plt.colorbar()
cb.set_label('Intensity')


if plots_save == 'Y':
  saveas = track_chosen + '.B_LAT.ps'
  plt.savefig(saveas, format='ps')

plt.show()

#######
# END
#######
