#
# Checking NIVAbasen data for Økokyst data
#

#
# 1 Libraries ++ ----
# 

library(niRvana)
library(dplyr)
library(ggplot2)
library(lubridate)

# Set username + password
set_credentials()


# 
# 2 Find hard bottom and soft bottom tables ---- 
#
get_nivabase_data("select TABLE_NAME from ALL_TAB_COLUMNS where OWNER = 'NIVADATABASE' and column_name = 'STATION_ID'")  
get_nivabase_data("select TABLE_NAME from ALL_TAB_COLUMNS where OWNER = 'NIVADATABASE' and column_name = 'GRAB_COLLECTION_ID'")  
get_nivabase_data("select TABLE_NAME from ALL_TAB_COLUMNS where OWNER = 'NIVADATABASE' and column_name = 'GRAB_ID'")  
get_nivabase_data("select TABLE_NAME from ALL_TAB_COLUMNS where OWNER = 'NIVADATABASE' and column_name = 'PARAMETER_ID'")  
get_nivabase_data("select TABLE_NAME from ALL_TAB_COLUMNS where OWNER = 'NIVADATABASE' and column_name = 'INDEX_ID'")  
get_nivabase_data("select TABLE_NAME from ALL_TAB_COLUMNS where OWNER = 'NIVADATABASE' and column_name = 'NIVA_TAXON_ID'")  
# 3         BB_GRAB_COLLECTION
# 4             BB_GRAB_SAMPLE
# 5              BB_PARAM_STAT
# 18      HB_OBSERVATION_SITES
# 19       HB_PARAMETER_VALUES
# 20             HB_PARAM_STAT

get_nivabase_data("select * from NIVADATABASE.STATIONS where rownum < 4")  # incl. SITE_ID STATION_ID SAMPLE_POINT_ID

# Hard-bottom
get_nivabase_data("select * from NIVADATABASE.HB_OBSERVATION_SITES where rownum < 4")  # incl. SITE_ID STATION_ID SAMPLE_POINT_ID
get_nivabase_data("select * from NIVADATABASE.HB_PARAMETER_VALUES where rownum < 4")   # incl. STATION_ID SAMPLE_DATE DEPTH1 DEPTH2 NIVA_TAXON_ID 
get_nivabase_data("select * from NIVADATABASE.HB_PARAM_STAT where rownum < 4")         # as above but FIRST and LAST, plus XMIN and XMAX
get_nivabase_data("select * from NIVADATABASE.HB_PARAMETER_DEFS")         # as above but FIRST and LAST, plus XMIN and XMAX
get_nivabase_data("select * from NIVADATABASE.HB_PARAMETERS_TAXA where rownum < 4")         # as above but FIRST and LAST, plus XMIN and XMAX
get_nivabase_data("select * from NIVADATABASE.HB_PARAMS_TAXA_INSERTS where rownum < 4")         # as above but FIRST and LAST, plus XMIN and XMAX
get_nivabase_data("select * from NIVADATABASE.TAXONOMY where rownum < 4")         # as above but FIRST and LAST, plus XMIN and XMAX

get_nivabase_data("select * from NIVADATABASE.HB_PARAMETER_VALUES where STATION_ID = 66198")   # incl. STATION_ID SAMPLE_DATE DEPTH1 DEPTH2 NIVA_TAXON_ID 

# Soft-bottom
get_nivabase_data("select * from NIVADATABASE.BB_GRAB_COLLECTION where rownum < 4") # GRAB_COLLECTION_ID, STATION_ID, SAMPLE_POINT_ID, SAMPLE_DATE 
get_nivabase_data("select * from NIVADATABASE.BB_GRAB_SAMPLE where rownum < 4")     # GRAB_COLLECTION_ID, GRAB_ID
get_nivabase_data("select * from NIVADATABASE.BB_PARAM_STAT where rownum < 4")      # as above but FIRST and LAST, plus XMIN and XMAX
get_nivabase_data("select * from NIVADATABASE.BB_INDEX_VALUES where rownum < 4")    # GRAB_COLLECTION_ID, GRAB_ID
get_nivabase_data("select * from NIVADATABASE.BB_SPECIE_COUNTS where rownum < 4")   # GRAB_ID BB_SPECIE_ID COUNT
get_nivabase_data("select * from NIVADATABASE.BB_INDEXES_DESCRIPTION")        
get_nivabase_data("select * from NIVADATABASE.V_BB_GRAB_COLL_SAMPLE where rownum < 4")        

#
# 3 Get data set with all projects ----
#

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


df <- get_nivabase_data("select * from NIVADATABASE.PROJECTS_STATIONS where STATION_NAME like 'Steilene%'") 
df_projects %>% filter(PROJECT_ID %in% df$PROJECT_ID) %>% View()
df_projects %>% filter(PROJECT_ID %in% df$PROJECT_ID) %>% pull(PROJECT_NAME)
df_projects %>% filter(PROJECT_ID %in% df$PROJECT_ID) %>% pull(PROJECT_ID)


#
# Other projects to check....?
#
grep("kyst", df_projects$PROJECT_NAME, value = TRUE, ignore.case = TRUE)
# "KystØstfold"                                                                    
# "Miljøtilstanden i Lillesands kystområder"                                       
# "Miljøtilstanden i Risørs kystområder før igangsetting av nytt renseanlegg¿"     
# "Miljøtilstanden i Tvedestrands kystområder før igangsetting av nytt renseanlegg"
# "Bløtbunnfaunaundersøkelser langs Vestfoldkysten"    


# grep("okokyst", df_projects$PROJECT_NAME, value = TRUE, ignore.case = TRUE)

# Get the name of the Økokyst project we are looking for

#
# 4 Check one project ----
# 

i <- 2
# Get a list of the stations
df_stations <- get_stations_from_project(proj_names[i], exact = TRUE)

#
# . a Hard-bottom ----
#

# Parameter definitions
df_hb_pardef <- get_nivabase_data("select * from NIVADATABASE.HB_PARAMETER_DEFS")

#
# . . a1 HB_PARAMETER_VALUES ----
#
df_hb_parvalues <- get_nivabase_selection(
  "STATION_ID, SAMPLE_DATE, DEPTH1, DEPTH2, PARAMETER_ID, NIVA_TAXON_ID, VALUE",
  "HB_PARAMETER_VALUES",
  "STATION_ID",
  df_stations$STATION_ID
)   
nrow(df_hb_parvalues)  # 5107

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
xtabs(~year(SAMPLE_DATE), df_hb_parvalues)  # 1990-2018, few data the last years
xtabs(~year(SAMPLE_DATE) + STATION_ID, df_hb_parvalues)

# Examples
df_hb_parvalues %>%
  filter(STATION_ID == 49907 & year(SAMPLE_DATE) == 2009)   # relativ forekomst per meter
df_hb_parvalues %>%
  filter(STATION_ID == 57623 & year(SAMPLE_DATE) == 2009)   # nedre voksegrense (per taxon) + antall arter, shannon index
df_hb_parvalues %>%
  filter(STATION_ID == 49907 & year(SAMPLE_DATE) == 2018)   # nedre voksegrense (overall)



#
# . . a2 HB_PARAM_STAT ----
#
df_hb_statvalues <- get_nivabase_selection(
  "*",
  "HB_PARAM_STAT",
  "STATION_ID",
  df_stations$STATION_ID
)  
nrow(df_hb_statvalues)  # 203

# Add defs
df_hb_statvalues <- df_hb_statvalues %>%
  left_join(df_hb_pardef %>% select(PARAMETER_ID, NAME, UNIT, DESCRIPTION))

xtabs(~year(FIRST), df_hb_statvalues)  
xtabs(~year(FIRST) + STATION_ID, df_hb_statvalues)

# Example
df_hb_statvalues %>%
  filter(STATION_ID == 49907 & year(FIRST) == 2016)   # nedre voksegrense (overall) - nEQR

#
# . b Soft-bottom (fauna indices) ----
#

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
nrow(df_bb_indexvalues)  # 9766

# Add defs
df_bb_indexvalues <- df_bb_indexvalues %>%
  left_join(df_bb_indexdef %>% select(INDEX_ID, INDEX_NAME))

# Add sample date and depth
df_bb_indexvalues <- df_bb_indexvalues %>%
  left_join(df_bb_grabcollection %>% select(GRAB_COLLECTION_ID, STATION_ID, SAMPLE_DATE, DEPTH))

# head(df_bb_indexvalues)

# Tables
xtabs(~year(SAMPLE_DATE), df_bb_indexvalues)  # 1990-2018, few data 2017-2018
xtabs(~year(SAMPLE_DATE) + STATION_ID, df_bb_indexvalues) # fewer data per station 2013-2018, change of stations 2017-2018

# Examples
# 15-16 indeces per grab sample
df_bb_indexvalues %>%
  filter(STATION_ID == 48480 & year(SAMPLE_DATE) == 2012) %>% xtabs(~GRAB_ID + GRAB_COLLECTION_ID, .)  # 2 grab collection, each with 8 grabs
df_bb_indexvalues %>%
  filter(STATION_ID == 48480 & year(SAMPLE_DATE) == 2016) %>% xtabs(~GRAB_ID + GRAB_COLLECTION_ID, .)  # 1 grab collection, with 4 grabs


#
# . c Sediment chemistry ----
#
df_chemistrydef <- get_nivabase_data("select * from NIVADATABASE.METHOD_DEFINITIONS")
df_sed_samplingmethods <- get_nivabase_data("select * from NIVADATABASE.SED_SAMPLING_METHODS")

df_sedsamples <- get_nivabase_selection(
  "STATION_ID, SAMPLE_ID, SAMPLE_DATE, WATER_DEPTH, SAMPLE_METHOD_ID",
  "SEDIMENT_SAMPLES",
  "STATION_ID",
  df_stations$STATION_ID
)   
nrow(df_sedsamples) # 22

df_sedslices <- get_nivabase_selection(
  "SAMPLE_ID, SLICE_ID, DEPTH1, DEPTH2",
  "SEDIMENT_SLICES",
  "SAMPLE_ID",
  df_sedsamples$SAMPLE_ID
)   
nrow(df_sedslices) # 40


df_sed_chemvalues <- get_nivabase_selection(
  "SLICE_ID, METHOD_ID, VALUE, FLAG1, DETECTION_LIMIT, QUANTIFICATION_LIMIT",
  "SEDIMENT_CHEMISTRY_VALUES",
  "SLICE_ID",
  df_sedslices$SLICE_ID
)   
nrow(df_sed_chemvalues) # 70

df_sed_chemdef <- get_nivabase_selection(
  "METHOD_ID, NAME, UNIT",
  "METHOD_DEFINITIONS",
  "METHOD_ID",
  unique(df_sed_chemvalues$METHOD_ID)
)   
nrow(df_sed_chemdef) # 7

# Add chemical names, depths, STATION_ID etc.
df_sed_chemvalues <- df_sed_chemvalues %>%
  left_join(df_sed_chemdef) %>%
  left_join(df_sedslices) %>%
  left_join(df_sedsamples)

# df_chemistrydef <- get_nivabase_data("select * from NIVADATABASE.METHOD_DEFINITIONS")
# df_sed_samplingmethods <- get_nivabase_data("select * from NIVADATABASE.SED_SAMPLING_METHODS")

#
# . d Combine results (except df_hb_statvalues) ----
#
colnames(df_hb_parvalues)
colnames(df_bb_indexvalues)
colnames(df_sed_chemvalues)

df_hb_parvalues <- df_hb_parvalues %>%
  left_join(df_stations %>% select(STATION_ID, STATION_CODE, STATION_NAME)) %>%
  select(STATION_CODE, STATION_NAME, SAMPLE_DATE, DEPTH1, DEPTH2, NAME, LATIN_NAME, VALUE, everything())

df_bb_indexvalues <- df_bb_indexvalues %>%
  left_join(df_stations %>% select(STATION_ID, STATION_CODE, STATION_NAME)) %>%
  select(STATION_CODE, STATION_NAME, GRAB_COLLECTION_ID, GRAB_ID, SAMPLE_DATE, DEPTH, INDEX_NAME, VALUE, everything())

df_sed_chemvalues <- df_sed_chemvalues %>%
  left_join(df_stations %>% select(STATION_ID, STATION_CODE, STATION_NAME)) %>%
  select(STATION_CODE, STATION_NAME, 
         SAMPLE_DATE, SAMPLE_ID, WATER_DEPTH, 
         SLICE_ID, DEPTH1, DEPTH2, 
         NAME, UNIT, VALUE, FLAG1, DETECTION_LIMIT, QUANTIFICATION_LIMIT, everything())





