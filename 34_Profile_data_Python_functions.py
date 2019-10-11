#
# Reads data from thredds
# This function is 'operated' by 'read_profile_list()' (see script '34_Profile_data_R_functions.R')
#
# Input:
#   input_lon, input_lat: - desired position for the data. The function finds the closest point in the model
#   variables: variables that we want, as a vector of strings. E.g. ('temp', 'salt')
#   filepath: url for server
#   
# Returns data as a Python dict, which becomes a named list in R  
# List components:
#   time_unix: numeric - time given as the number of seconds since 01-01-1970 
#   variables: a 1-dim array of variable names (same as those give as input). To find names, use variable_info(), see below  
#   values: a list with one item per variable. Each item is a matrix (depth x time) of values
#   i,j: indices for x,y-position in the model 
#   lon,lat: actual position for the x,y-position in the model  
#   Z: bottom and top of each depth layer, starting from the bottom. Length = number of layers plus 1
# 
def read_profile_python(input_lon, input_lat, variables, filepath):
  import numpy as np     # Package for scientific computing 
  import pandas as pd    # Package for data frames 
  import matplotlib.pyplot as plt   
  from datetime import datetime,timedelta
  from netCDF4 import Dataset #  This is handy for working with netCDF files
  import roppy

  # Connect to server
   # The OPENDAP URL
  filehandle = Dataset(filepath) # open for reading 
  grid = roppy.SGrid(filehandle) # Create a grid object for our file - used to get depth layers

  # Create handle for variables
  lon = filehandle.variables['lon_rho']
  lat = filehandle.variables['lat_rho']
  ocean_time = filehandle.variables['ocean_time']

  # We need to find the grid location closest to the selected coordinates.
  # For the sake of simplicity, let's calculate the sum of the absolute 
  # differences between all grid points latitude and longitude and the 
  # selected coordinates. 
  position_diff = np.abs( lat[:] - input_lat ) + np.abs( lon[:] - input_lon )

  # This line will find the indices of the minimum value in 
  i, j = np.unravel_index( position_diff.argmin(), position_diff.shape )

  print('Grid indices of grid point closest to selected location: {}, {}\n'.format(i, j))
  print('Grid point longitude: {}'.format(lon[i,j]))
  print('Grid point latitude: {}'.format(lat[i,j]))

  # Depth layer (uses the grid object created by roppy.SGrid)
  Z = grid.z_w[: ,i,j]
  
  # max of the time dimension
  max_time = len(filehandle.dimensions["ocean_time"])

  # Pick times (for now, we choose all of them)
  i_time = np.arange(0,max_time,1) 

  # temp = potential temperaturen
  # salt = salinity
  # N1_p = phosphate/phosphorus
  # N3_n = nitrate/nitrogen
  # chla = chlorophyll A

  loc_time = filehandle.variables['ocean_time'][i_time]
  
  # loc_variable = filehandle.variables[variable_name][i_time,:,i,j]
  
  value_list = []
  for var in variables:
    values = [filehandle.variables[var][i_time,:,i,j]]
    value_list.append(values)

  # Turn times into datetimes
  # expected_format = 'seconds since %Y-%m-%d %H:%M:%S'
  # timeref = datetime.strptime(ocean_time.units, expected_format) # Creating a datetime object for the reference time
  # datetime_list = np.array([ timeref + timedelta(seconds=t) for t in loc_time  ])
  # loc_datetime = np.array(datetime_list)

  # result = {'time_unix':loc_time, 'values':loc_variable, 'i':i, 'j':j, 'lon':lon[i,j], 'lat':lat[i,j], 'Z':Z}
  result = {'time_unix':loc_time, 'variables':variables, 'values':value_list, 'i':i, 'j':j, 'lon':lon[i,j], 'lat':lat[i,j], 'Z':Z}

  return result


#
# Gets unit for a given variable
# Used by variable_info()
# Needs its own function, since not all variables have 'unit' (in contrast to 'name' and 'long_name')
#
def get_unit(variable):
  try:
    unit = filehandle.variables[variable].units
  except:
    unit = "NA"
  return unit

#
# Returns a pandas DataFrame (in R: data.frame) with 3 variables: name, long name and unit  
# Is used to fine 'name'
#
def variable_info():
  import numpy as np
  import pandas as pd
  from netCDF4 import Dataset #  This is handy for working with netCDF files

  # Connect to server
  filepath = 'http://thredds.met.no/thredds/dodsC/metusers/arildb/MARTINI800_prov_v2.ncml' # The OPENDAP URL
  filehandle = Dataset(filepath) # open for reading 
  keys = filehandle.variables.keys()
  print('Number of variables:', len(keys))
  print()

  names = [(filehandle.variables[key].name) for key in keys]
  #names = np.array(names)

  long_names = [(filehandle.variables[key].long_name) for key in keys]
  #long_names = np.array(long_names)
  
  units = [get_unit(key) for key in keys]
  
  # result = pd.DataFrame({'keys':keys, 
  #                        'long_names':long_names})
  result = pd.DataFrame({'name':np.array(names), 'long_name':np.array(long_names), 'unit':np.array(units)})
  
  #
  # UNCOMMENT TO SEE THE WHOLE LIST 
  #
  # for key in keys:
  #   print(key, '=', filehandle.variables[key].long_name)
  
  return result


