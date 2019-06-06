#
# Getting hard- and soft-bottom Økokyst data
#

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 1 Libraries ++ ----
# 
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

library(niRvana)
library(dplyr)
library(ggplot2)
library(lubridate)

# Summarizes for instance years into a string - see example
summarize_sequence <- function(x){
  x <- sort(unique(x))
  dx <- diff(x)
  df <- tibble(
    x = x,
    index = cumsum(c(1, dx) > 1) + 1)
  df %>% 
    group_by(index) %>%
    summarize(Min = min(x),Max = max(x)) %>%
    mutate(Summ = ifelse(Min < Max, paste0(Min,"-",Max), Min)) %>%
    ungroup() %>%
    summarize(Summ = paste0(Summ, collapse = ",")) %>%
    pull(Summ)
}


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
# Set username + password
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

set_credentials()


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 2 Get projects ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

df_projects <- get_projects()   # we call it 'df_projects' (the default name used by 'get_stations_from_project')

# Search for Økokyst projects
proj_names <- grep("økokyst", df_projects$PROJECT_NAME, value = TRUE, ignore.case = TRUE)

# Filter for Skagerrak
proj_names <- grep("skagerrak",proj_names, value = TRUE, ignore.case = TRUE)

# Add extras (KYO + indre and ytre Oslofjord)
proj_names <- c(
  proj_names, 
  "Nordsjøen Nord",
  grep("indre oslo", df_projects$PROJECT_NAME, value = TRUE, ignore.case = TRUE), 
  grep("ytre oslo", df_projects$PROJECT_NAME, value = TRUE, ignore.case = TRUE), 
  grep("kyo", df_projects$PROJECT_NAME, value = TRUE, ignore.case = TRUE)
)


# get PROJECT_ID
proj_id <- df_projects %>% 
  filter(PROJECT_NAME %in% proj_names) %>%
  pull(PROJECT_ID)

# Get extra PROJECT_IDs from Steilene (to find Indre Oslofjord projects)
df <- get_nivabase_data("select * from NIVADATABASE.PROJECTS_STATIONS where STATION_NAME like 'Steilene%'") 
proj_id_extra <- df_projects %>% filter(PROJECT_ID %in% df$PROJECT_ID) %>% pull(PROJECT_ID)

proj_id <- c(proj_id, proj_id_extra) %>% unique()

#
# Other projects to check....?
# grep("kyst", df_projects$PROJECT_NAME, value = TRUE, ignore.case = TRUE)
#

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 3. Get stations
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

df_stations <- get_nivabase_selection(
  "PROJECT_ID, STATION_ID, STATION_CODE, STATION_NAME",
  "PROJECTS_STATIONS",
  "PROJECT_ID",
  proj_id
)   

df_stations_geomid <- get_nivabase_selection(
  "STATION_ID, GEOM_REF_ID",
  "STATIONS",
  "STATION_ID",
  df_stations$STATION_ID
)   

df_stations_pos <- get_nivabase_selection(
  "SAMPLE_POINT_ID, LATITUDE, LONGITUDE",
  "SAMPLE_POINTS",
  "SAMPLE_POINT_ID",
  df_stations_geomid$GEOM_REF_ID,
  owner = "NIVA_GEOMETRY"
)   

df_stations <- df_stations %>% 
  left_join(df_stations_geomid) %>%
  left_join(df_stations_pos, by = c("GEOM_REF_ID" = "SAMPLE_POINT_ID"))

# df_stations %>% filter(STATION_CODE == "514")
df_stations <- df_stations %>%
  filter(LONGITUDE > 7.34 & LATITUDE < 60)

nrow(df_stations)  # 1247

# Make interactive map if you wish 
make_map <- FALSE
if (make_map){
  library(leaflet)
  leaflet() %>%
    addTiles() %>%  # Default OpenStreetMap map tiles
    addMarkers(lng = df_stations$LONGITUDE, lat = df_stations$LATITUDE,
               popup = paste(df_stations$STATION_CODE, df_stations$STATION_NAME))
}

# df_stations %>% filter(STATION_NAME == "Nakkholmen")
  
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 4. Get hard-bottom data (HB_PARAMETER_VALUES only) ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Parameter definitions
df_hb_pardef <- get_nivabase_data("select * from NIVADATABASE.HB_PARAMETER_DEFS")


df_hb_parvalues <- get_nivabase_selection(
  "STATION_ID, SAMPLE_DATE, DEPTH1, DEPTH2, PARAMETER_ID, NIVA_TAXON_ID, VALUE",
  "HB_PARAMETER_VALUES",
  "STATION_ID",
  df_stations$STATION_ID
)   
nrow(df_hb_parvalues)  # 15054

df_hb_parvalues %>% filter(STATION_ID == 66198)

test <- get_nivabase_selection(
  "STATION_ID, SAMPLE_DATE, DEPTH1, DEPTH2, PARAMETER_ID, NIVA_TAXON_ID, VALUE",
  "HB_PARAMETER_VALUES",
  "STATION_ID",
  66198
)   


# Add defs
df_hb_parvalues <- df_hb_parvalues %>%
  left_join(df_hb_pardef %>% select(PARAMETER_ID, NAME, UNIT, DESCRIPTION))

# Add species
species <- unique(df_hb_parvalues$NIVA_TAXON_ID)
species <- species[!is.na(species)] 

df_species <- get_nivabase_selection(
  "NIVA_TAXON_ID, LATIN_NAME",
  "TAXONOMY",
  "NIVA_TAXON_ID",
  species
)   

df_hb_parvalues <- df_hb_parvalues %>%
  left_join(df_species)

# Tables
# xtabs(~year(SAMPLE_DATE), df_hb_parvalues)  # 1990-2018, few data the last years
# xtabs(~year(SAMPLE_DATE) + STATION_ID, df_hb_parvalues)


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 5.Get soft-bottom data (fauna indices) ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Parameter definitions
df_bb_indexdef <- get_nivabase_data("select * from NIVADATABASE.BB_INDEXES_DESCRIPTION")

df_bb_grabcollection <- get_nivabase_selection(
  "GRAB_COLLECTION_ID, STATION_ID, SAMPLE_POINT_ID, SAMPLE_DATE, DEPTH",
  "BB_GRAB_COLLECTION",
  "STATION_ID",
  df_stations$STATION_ID
)   
nrow(df_bb_grabcollection)

df_bb_grabsample <- get_nivabase_selection(
  "GRAB_COLLECTION_ID, GRAB_ID, GRAB_NR, SAMPLE_POINT_ID, SAMPLE_AREA",
  "BB_GRAB_SAMPLE",
  "GRAB_COLLECTION_ID",
  df_bb_grabcollection$GRAB_COLLECTION_ID
)   
nrow(df_bb_grabsample)

df_bb_indexvalues <- get_nivabase_selection(
  "GRAB_COLLECTION_ID, GRAB_ID, INDEX_ID, NIVA_TAXON_ID, VALUE, VALUE_ID",
  "BB_INDEX_VALUES",
  "GRAB_ID",
  df_bb_grabsample$GRAB_ID
)   
nrow(df_bb_indexvalues)  # 34128

# Add defs
df_bb_indexvalues <- df_bb_indexvalues %>%
  left_join(df_bb_indexdef %>% select(INDEX_ID, INDEX_NAME))

# Add sample date and depth
df_bb_indexvalues <- df_bb_indexvalues %>%
  left_join(df_bb_grabcollection %>% select(GRAB_COLLECTION_ID, STATION_ID, SAMPLE_DATE, DEPTH))

# head(df_bb_indexvalues)

# Tables
# xtabs(~year(SAMPLE_DATE), df_bb_indexvalues)  # 1990-2018, few data 2017-2018
# xtabs(~year(SAMPLE_DATE) + STATION_ID, df_bb_indexvalues) # fewer data per station 2013-2018, change of stations 2017-2018


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 6. Get sediment chemistry ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

df_sed_samplingmethods <- get_nivabase_data("select * from NIVADATABASE.SED_SAMPLING_METHODS")

df_sedsamples <- get_nivabase_selection(
  "STATION_ID, SAMPLE_ID, SAMPLE_DATE, WATER_DEPTH, SAMPLE_METHOD_ID",
  "SEDIMENT_SAMPLES",
  "STATION_ID",
  df_stations$STATION_ID
)   
nrow(df_sedsamples) # 66

df_sedslices <- get_nivabase_selection(
  "SAMPLE_ID, SLICE_ID, DEPTH1, DEPTH2",
  "SEDIMENT_SLICES",
  "SAMPLE_ID",
  df_sedsamples$SAMPLE_ID
)   
nrow(df_sedslices) # 96


df_sed_chemvalues <- get_nivabase_selection(
  "SLICE_ID, METHOD_ID, VALUE, FLAG1, DETECTION_LIMIT, QUANTIFICATION_LIMIT",
  "SEDIMENT_CHEMISTRY_VALUES",
  "SLICE_ID",
  df_sedslices$SLICE_ID
)   
nrow(df_sed_chemvalues) # 8408

df_sed_chemdef <- get_nivabase_selection(
  "METHOD_ID, NAME, UNIT",
  "METHOD_DEFINITIONS",
  "METHOD_ID",
  unique(df_sed_chemvalues$METHOD_ID)
)   
nrow(df_sed_chemdef) # 8

# Add chemical names, depths, STATION_ID etc.
df_sed_chemvalues <- df_sed_chemvalues %>%
  left_join(df_sed_chemdef) %>%
  left_join(df_sedslices) %>%
  left_join(df_sedsamples)

# df_chemistrydef <- get_nivabase_data("select * from NIVADATABASE.METHOD_DEFINITIONS")
# df_sed_samplingmethods <- get_nivabase_data("select * from NIVADATABASE.SED_SAMPLING_METHODS")

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 7. Make results ready for export ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

df_hb_parvalues <- df_hb_parvalues %>%
  left_join(df_stations %>% select(STATION_ID, STATION_CODE, STATION_NAME, LONGITUDE, LATITUDE, PROJECT_ID)) %>%
  left_join(df_projects %>% select(PROJECT_ID, PROJECT_NAME)) %>%
  mutate(YEAR = year(SAMPLE_DATE), MONTH = month(SAMPLE_DATE)) %>%
  select(STATION_CODE, STATION_NAME, PROJECT_ID, PROJECT_NAME, YEAR, MONTH, SAMPLE_DATE, DEPTH1, DEPTH2, NAME, LATIN_NAME, VALUE, everything())

df_hb_parvalues_since2014 <- df_hb_parvalues %>%
  filter(YEAR >= 2014)

df_bb_indexvalues <- df_bb_indexvalues %>%
  left_join(df_stations %>% select(STATION_ID, STATION_CODE, STATION_NAME, LONGITUDE, LATITUDE, PROJECT_ID)) %>%
  left_join(df_projects %>% select(PROJECT_ID, PROJECT_NAME)) %>%
  mutate(YEAR = year(SAMPLE_DATE), MONTH = month(SAMPLE_DATE)) %>%
  select(STATION_CODE, STATION_NAME, PROJECT_ID, PROJECT_NAME, GRAB_COLLECTION_ID, GRAB_ID, YEAR, MONTH, SAMPLE_DATE, DEPTH, INDEX_NAME, VALUE, everything())

df_sed_chemvalues <- df_sed_chemvalues %>%
  left_join(df_stations %>% select(STATION_ID, STATION_CODE, STATION_NAME, LONGITUDE, LATITUDE, PROJECT_ID)) %>%
  left_join(df_projects %>% select(PROJECT_ID, PROJECT_NAME)) %>%
  mutate(YEAR = year(SAMPLE_DATE), MONTH = month(SAMPLE_DATE)) %>%
  select(STATION_CODE, STATION_NAME, PROJECT_ID, PROJECT_NAME, 
         SAMPLE_DATE, SAMPLE_ID, WATER_DEPTH, 
         SLICE_ID, DEPTH1, DEPTH2, 
         NAME, UNIT, VALUE, FLAG1, DETECTION_LIMIT, QUANTIFICATION_LIMIT, everything())

# number of rows
for (df in list(df_hb_parvalues, df_hb_parvalues_since2014, df_bb_indexvalues, df_sed_chemvalues)) 
  print(nrow(df))

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 8. Check data ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# df <- df_hb_parvalues
df <- df_hb_parvalues_since2014
xtabs(~STATION_NAME + PROJECT_NAME, df)

df %>%
  filter(STATION_NAME %in% "Torgersøy") 

df_map <- df_hb_parvalues_since2014 %>%
  mutate(Year = year(SAMPLE_DATE)) %>%
  group_by(STATION_CODE, STATION_NAME) %>%
  summarize(LONGITUDE = mean(LONGITUDE), LATITUDE = mean(LATITUDE),
            Years = summarize_sequence(Year),
            No_species = length(unique(NIVA_TAXON_ID))) %>%
  mutate(Popup_text = paste0(STATION_CODE, ": ", STATION_NAME, "<br>", Years, "<br>", No_species, " species"))

library(leaflet)
leaflet() %>%
  addTiles() %>%  # Default OpenStreetMap map tiles
  addMarkers(lng = df_map$LONGITUDE, lat = df_map$LATITUDE,
             popup = df_map$Popup_text)


df_hb_parvalues %>%
  filter(STATION_NAME %in% "Vindvik" & year(SAMPLE_DATE) == 2018)   # nedre voksegrense (overall)



xtabs(~INDEX_NAME, df_bb_indexvalues)

df_map <- df_hb_parvalues_since2014 %>%
  mutate(Year = year(SAMPLE_DATE)) %>%
  group_by(STATION_CODE, STATION_NAME) %>%
  summarize(LONGITUDE = mean(LONGITUDE), LATITUDE = mean(LATITUDE),
            Years = summarize_sequence(Year),
            No_indices = length(unique(NAME))) %>%
  mutate(Popup_text = paste0(STATION_CODE, ": ", STATION_NAME, "<br>", Years, "<br>", No_indices, " indices"))

library(leaflet)
leaflet() %>%
  addTiles() %>%  # Default OpenStreetMap map tiles
  addMarkers(lng = df_map$LONGITUDE, lat = df_map$LATITUDE,
             popup = df_map$Popup_text)



#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 9. Save/export ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

openxlsx::write.xlsx(df_hb_parvalues, "Test_data/data_hardbottom.xlsx")
openxlsx::write.xlsx(df_hb_parvalues_since2014, "Test_data/data_hardbottom_since2017.xlsx")
openxlsx::write.xlsx(df_bb_indexvalues, "Test_data/data_softbottom.xlsx")
openxlsx::write.xlsx(df_sed_chemvalues, "Test_data/data_sedimentchemistry.xlsx")

# Copied to 
# \\niva-of5\OSL-Data-NIVA\Avdeling\Vass\316_Miljøinformatikk\Prosjekter\180116_Martini\Datasets


