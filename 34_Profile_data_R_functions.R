#
# Read profile data from OpenDAP server
# Returns a data frame in 'long' format (suitable for ggplot)
#
read_profile <- function(input_lon, input_lat, 
                         variable_names,
                         server_url = 'http://thredds.met.no/thredds/dodsC/metusers/arildb/MARTINI800_prov_v2.ncml'){
  X <- read_profile_list(10.5268, 59.0267, c('temp','salt'), server_url)
  profile_list2dataframe(X)
}


#
# Read data using Python function
# Returns list (including a time variable)
#
read_profile_list <- function(input_lon, input_lat, 
                              variable_names, 
                              server_url = 'http://thredds.met.no/thredds/dodsC/metusers/arildb/MARTINI800_prov_v2.ncml'){
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


