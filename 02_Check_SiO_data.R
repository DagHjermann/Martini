library(niRvana)
library(dplyr)
library(lubridate)
library(ggplot2)

set_credentials()
get_nivabase_data("select * from NIVADATABASE.PROJECTS_STATIONS where rownum < 4")   

# Ger station ID
get_nivabase_data("select * from NIVADATABASE.PROJECTS_STATIONS where STATION_NAME like 'Arendal St%'")   
st_id <- 42403

# Get samples
get_nivabase_data("select * from NIVADATABASE.WATER_SAMPLES where rownum < 4")   
samples <- get_nivabase_selection(
  "WATER_SAMPLE_ID, STATION_ID, SAMPLE_DATE, DEPTH1, DEPTH2",
  "WATER_SAMPLES",
  "STATION_ID",
  st_id
)
nrow(samples)

# Get silicate methods from WC_PARAMETER_DEFINITIONS (not used)
df <- get_nivabase_data("select * from NIVADATABASE.WC_PARAMETER_DEFINITIONS where NAME like 'Si%'")   
df

# Get silicate methods from METHOD_DEFINITIONS
df2 <- get_nivabase_data("select * from NIVADATABASE.METHOD_DEFINITIONS where NAME like 'SiO%'")   
df2
method_ids <- df2$METHOD_ID

# Get all silicate data (Arendal or not)
df_data <- get_nivabase_selection(
  "WATER_SAMPLE_ID, METHOD_ID,VALUE",
  "WATER_CHEMISTRY_VALUES",
  "METHOD_ID",
  method_ids
)
nrow(df_data)

# Pick silicate data from Arendal
df_data2 <- df_data %>%
  filter(WATER_SAMPLE_ID %in% samples$WATER_SAMPLE_ID) %>%
  left_join(samples)
  
nrow(df_data2)

# Explore
df_data2 %>% xtabs(~year(SAMPLE_DATE), .)
df_data2 %>% xtabs(~year(SAMPLE_DATE) + DEPTH1, .)
df_data2 %>% xtabs(~year(SAMPLE_DATE) + METHOD_ID, .)

# Check methods used at Arendal
df2 %>% filter(METHOD_ID %in% c(7973,11668))
# METHOD_ID NAME UNIT                           LABORATORY METHOD_REF
# 1      7973 SiO2   µM Havforskningsinstituttet, Flødevigen       <NA>
# 2     11668 SiO3   µM Havforskningsinstituttet, Flødevigen       <NA>
#   DESCR
# 1 Opprettet i forbindelse med import av data fra Kystovervåkingen, 1990-2007. Merk: Parameternavnet i den opprinnelige basen er SiO3, men dataene importeres med navnet SiO2 etter beskjed fra Birger Bjerkeng og Jan Magnusson.
# 2                                                                                                                                                                                                                           <NA>
#   ENTERED_BY        ENTERED_DATE MATRIX CAS IUPAC BZNO BASIS_ID SUBSTANCE_ID MATRIX_ID DETECTION_LIMIT UNCERTAINTY QUANTIFICATION_LIMIT ANALYSIS STANDARD
# 1        ABM 2008-10-09 14:41:36   <NA>  NA    NA   NA       NA          631        NA              NA          NA                   NA     <NA>     <NA>
# 2        JAN 2011-04-28 13:44:26   <NA>  NA    NA   NA       NA           NA        NA              NA          NA                   NA     <NA>     <NA>
  
# Test plot in year where both methods are used
# Spoiler: the data are identical for SiO2 and SiO3 (not wrong, as micromol is used)
df_data2 %>%
  filter(year(SAMPLE_DATE) %in% 2005 & DEPTH1 %in% 0) %>%
  ggplot(aes(SAMPLE_DATE, VALUE, color = factor(METHOD_ID))) + 
  geom_line() +
  facet_grid(. ~ METHOD_ID)

# Check some data
df_data2 %>%
  filter(year(SAMPLE_DATE) %in% 2005) %>% arrange(METHOD_ID, SAMPLE_DATE, DEPTH1) %>% 
  head(20)
df_data2 %>%
  filter(year(SAMPLE_DATE) %in% 2011) %>% arrange(METHOD_ID, SAMPLE_DATE, DEPTH1) %>% 
  head(20)
df_data2 %>%
  filter(year(SAMPLE_DATE) %in% 1991) %>% arrange(METHOD_ID, SAMPLE_DATE, DEPTH1) %>% 
  head(20)


  