#
# Read profile data from OpenDAP server
# Returns a data frame in 'long' format (suitable for ggplot)
#
# This function uses the functions 'read_profile_list' (which in turn uses the Python function 'read_profile_python')
#   and 'profile_list2dataframe'
# This function is used by the function 'read_profiles' (which calls 'read_profile' repeatedly)
#
read_profile <- function(input_lon, input_lat, 
                         variable_names,
                         server_url = 'https://thredds.met.no/thredds/dodsC/metusers/arildb/MARTINI800_prov_v2.ncml'){
  X <- read_profile_list(input_lon, input_lat, variable_names, server_url)
  profile_list2dataframe(X)
}


#
# Read data using Python function
# Returns list (including a time variable)
#
read_profile_list <- function(input_lon, input_lat, 
                              variable_names, 
                              server_url = 'https://thredds.met.no/thredds/dodsC/metusers/arildb/MARTINI800_prov_v2.ncml'){
  result <- read_profile_python(input_lon, input_lat, variable_names, server_url)
  # Remove one level from value list
  for (i in seq_along(result$values))
    result$values[[i]] <- result$values[[i]][[1]]
  # Set names of value list
  names(result$values) <- result$variables
  # Make time variable ready to use in plots
  result$time <- as.POSIXct(result$time_unix, origin = "1970-01-01", tz = "UTC")
  result
}


#
# Takes list returned by 'read_profile_list()' and returns a data frame
#   in 'long' format (suitable for ggplot)
#
profile_list2dataframe <- function(list){
  depth_lo <- head(list$Z,-1)              # all values except the last
  depth_hi <- tail(list$Z,-1)              # all values except the first
  depth_mid <- (depth_lo + depth_hi)/2     # mid depth
  
  n_var <- length(list$values)
  df <- matrix(NA, nrow = length(list$values[[1]]), ncol = n_var)
  for (i in 1:n_var)
    df[,i] <- as.numeric(list$values[[i]])
  df <- as.data.frame(df)
  colnames(df) <- list$variables
  
  # Make data frame
  data.frame(
    Time = list$time,                                       # time has length 364, but will just repeat itself
    Depth_lo = rep(depth_lo, each = length(list$time)),     # value no. 1 is repeteated 364 times, etc.
    Depth_hi = rep(depth_hi, each = length(list$time)),     #    "
    Depth_mid = rep(depth_mid, each = length(list$time)),   #    "
    df
  )
  
}

#
# Read profile data for several stations  
# 
# Compulsory input:
#   input_data: data frame with at least 3 variables: one 'ID' column (e.g., station name), longitude and latitude
#   id_column: the name of the variable to be used for ID (e.g. id_column = "Station")
#   variable_names: text vector of variable names tha twe want to extract
# Optional input:
#   lon_column: if not given, the function tries to search for a variable name including 'lon'; if not found, it returns an error message
#   lat_column: if not given, the function tries to search for a variable name including 'lat'; if not found, it returns an error message
#   server_url: URL to the thredds server
# Returns:
#   A data frame similar to the one returned by 'read_profile()', but with an extra colymn named 'ID'
#
read_profiles <- function(input_data,  
                          id_column,
                          variable_names,
                          lon_column = NULL,
                          lat_column = NULL,
                          server_url = 'https://thredds.met.no/thredds/dodsC/metusers/arildb/MARTINI800_prov_v2.ncml'){
  
    k <- grep("lon", colnames(input_data), ignore.case = TRUE)[1]
    if (length(k) == 0)
      stop("You need to specify the name of the variable with longitude using lon_column")
    lon_column <- colnames(input_data)[k]
    
    k <- grep("lat", colnames(input_data), ignore.case = TRUE)[1]
    if (length(k) == 0)
      stop("You need to specify the name of the variable with latitude using lon_column")
    lat_column <- colnames(input_data)[k]
    
    result <- 1:nrow(input_data) %>% map_df(~profile_one_station(., input_data, id_column, lon_column, lat_column))
    
    result
}

  

# Function for getting data from a single station - for use in 'read_profiles()'   
# - the main reason that we make a special function for this is to add a variable for Station
# The function's input is the row number 'i' of 'df_stations'
# The function returns a data frame with the data (using read_profile) and an extra variable for Station  
profile_one_station <- function(i, input_data, id_column, lon_column, lat_column){
  
    data <- 
    # Get data 
    read_profile(input_lon = df_stations[i, lon_column], 
                 input_lat = df_stations[i, lat_column], 
                 c('N1_p','N3_n')
    ) %>%
    # Add 'ID' variable 
    mutate(ID = df_stations[i, id_column]) %>%
    # Puts ID variable as variable number 1 (not really needed)
    select(ID, everything())
  
  # The function returns this as result
  data
  
}


