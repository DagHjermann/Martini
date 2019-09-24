
# Getting data 
def get_surface_data(station_code, long_sel, lat_sel, day_interval):

  # We need to find the grid location closest to Torbjornskjaer.
  # For the sake of simplicity, let's calculate the sum of the absolute 
  # differences between all grid points latitude and longitude and the 
  # Torbjornskjaer coordinates. 
  
  position_diff = np.abs( lat[:] - lat_sel) + np.abs( lon[:] - long_sel)
  
  # This line will find the indices of the minimum value in 
  i, j = np.unravel_index( position_diff.argmin(), position_diff.shape )

  # Defining the times for which we want data
  # E.g. if day_interval = 30, we get 0, 30, 60, ... , 360
  times = range(0, 364, int(day_interval))   # need to 'make sure' day_interval is an integer value
  
  # return times
  #return filehandle.variables['ocean_time'][times]

  # temp = potential temperaturen
  # salt = salinity
  # N1_p = phosphate/phosphorus
  # N3_n = nitrate/nitrogen
  # chla = chlorophyll A
  
  # The data are picked from the top layer 
  #   (layer number 0 in the model)
  # We make the function return a pandas DataFrame,
  #   which in R will automatically be converted to a R data.frame
  data_loc = pd.DataFrame({
    'STATION_CODE': station_code,
    'i': i,
    'j': j,
    'Time_num': filehandle.variables['ocean_time'][times],
    'Temp': filehandle.variables['temp'][times,41,i,j],
    'Salinity': filehandle.variables['salt'][times,41,i,j],
    'PO4': filehandle.variables['N1_p'][times,41,i,j],
    'NO3': filehandle.variables['N3_n'][times,41,i,j]
    })

  return data_loc
  
