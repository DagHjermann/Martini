#
# Checking data exported from Vannmiljø
#
# 1. VannmiljoEksport_2019-05-27_okokyst.xlsx              - all Økokyst data (NIVA and IMR)
# 2. VannmiljoEksport_2019-05-27_havforsuring.xlsx         - all Havforsuring data
# 3. VannmiljoEksport_2019-05-27_indre_ytre_oslofjord.xlsx - data from "Oslofjorden" (selected i station part) 
#                                                            and "Tiltaksrettet overvåkning" 

library(readxl)
library(dplyr)
library(ggplot2)
library(lubridate)
library(sp)        # SpatialPoints(), CRS(), spTransform()

crs_longlat <- "+proj=longlat"
crs_utm <- "+proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m"

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 1. Read and fix data ----
# Okokyst data (including Økokyst klima from IMR) + Oslofjorden (Indre + Ytre)
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

dat_okokyst <- read_excel("Input_data/VannmiljoEksport_2019-05-27_okokyst.xlsx", col_types = "text")
dat_of <- read_excel("Input_data/VannmiljoEksport_2019-05-27_indre_ytre_oslofjord.xlsx", col_types = "text")

# Combine
dat <- rbind(dat_okokyst, dat_of)

# Seems that dyp has been set to NA where Dyp is really zero
dat$Ovre_dyp[is.na(dat$Ovre_dyp)] <- 0
dat$Nedre_dyp[is.na(dat$Nedre_dyp)] <- 0

# Convert some character variables to numbers (first have to change decimal comma to 'full stop')

vars <- c("Verdi", "Ovre_dyp", "Nedre_dyp", "UTM33 Ost (X)", "UTM33 Nord (Y)")

for (var in vars){
  x <- as.numeric(sub(",", ".", dat[[var]], fixed = TRUE))
  cat(var, "- check that number is zero:", sum(is.na(x)), "\n")  # 0, so all ar ok
  dat[[var]] <- x
  }

# Make time variable
dat$Time <- ymd_hms(dat$Tid_provetak)
dat$Year <- year(dat$Time)
dat$Month <- month(dat$Time)

# Fix Siktedyp depth to zero
sel <- dat$Parameter_navn == "Siktedyp"; sum(sel); # table(dat$Nedre_dyp[sel])
dat$Nedre_dyp[sel] <- 0

#
# Get ordinary long/lat (transform from UTM zone 33 to Long,Lat)
#

SP.utm <- SpatialPoints(dat[,c("UTM33 Ost (X)", "UTM33 Nord (Y)")], 
                        proj4string=CRS(crs_utm)
)
SP.longlat <- spTransform(SP.utm, CRS(crs_longlat))
dat$Long <- SP.longlat@coords[,1]
dat$Lat <- SP.longlat@coords[,2]

#
# Spatial filter
#

dat <- dat %>%
  filter(Long > 7.5 & Lat < 60)

#
# Test plot
#

map_df <- dat %>% 
  count(Vannlokalitet, Long, Lat) %>% 
  mutate(Popup = paste0(Vannlokalitet, "<br>Long,lat = ", round(Long, 4), ", ", round(Lat, 4), ""))

library(leaflet)
leaflet() %>%
  addTiles() %>%  # Default OpenStreetMap map tiles
  addMarkers(lng = map_df$Long, lat = map_df$Lat,
             popup = map_df$Popup)

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 2. Hard bottom data ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Get station list
stations_hardbottom <- dat %>%
  filter(Parameter_navn %in% c("Relativ forekomst makroalger (skala)") |
           grepl("Maksdypindeks makroalger", Parameter_navn)) %>%
  group_by(Vannlokalitet_kode, Vannlokalitet, Long, Lat)  %>%
  summarize(Year_min = min(Year), Year_max = max(Year), 
            Oppdragstaker = paste(unique(Oppdragstaker), collapse = ";")
            )
stations_hardbottom

# Pick all data that are in station list
data_hardbottom <- dat %>%
  filter(Vannlokalitet_kode %in% stations_hardbottom$Vannlokalitet_kode,
         Vannlokalitet %in% stations_hardbottom$Vannlokalitet)

# Positions
library(leaflet)
leaflet() %>%
  addTiles() %>%  # Default OpenStreetMap map tiles
  addMarkers(lng = stations_hardbottom$Long, 
             lat = stations_hardbottom$Lat,
             popup = stations_hardbottom$Vannlokalitet)

table(data_hardbottom$Oppdragstaker) # NIVA + NIVA/HI

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 3. Soft bottom data ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Get station list
stations_softbottom <- dat %>%
  filter(grepl("diversitetsindeks", Parameter_navn)) %>%
  group_by(Vannlokalitet_kode, Vannlokalitet, Long, Lat)  %>%
  summarize(Year_min = min(Year), Year_max = max(Year), 
            Oppdragstaker = paste(unique(Oppdragstaker), collapse = ";")
  )
stations_softbottom

# Pick all data that are in station list
data_softbottom <- dat %>%
  filter(Vannlokalitet_kode %in% stations_softbottom$Vannlokalitet_kode,
         Vannlokalitet %in% stations_softbottom$Vannlokalitet)

# Positions
library(leaflet)
leaflet() %>%
  addTiles() %>%  # Default OpenStreetMap map tiles
  addMarkers(lng = stations_softbottom$Long, 
             lat = stations_softbottom$Lat,
             popup = stations_softbottom$Vannlokalitet)

table(data_softbottom$Oppdragstaker) # also DnV (and a little bit Norconsult)

xtabs(~Parameter_navn + Oppdragstaker, data_softbottom)

# Check data
xtabs(~Parameter_navn, data_softbottom %>% filter(Oppdragstaker %in% "DnV"))
xtabs(~Parameter_navn + Vannlokalitet, data_softbottom %>% filter(Oppdragstaker %in% "DnV"))
xtabs(~Parameter_navn, data_softbottom %>% filter(Oppdragstaker %in% "Norconsult AS"))
xtabs(~Parameter_navn, data_softbottom %>% filter(grepl("NIVA", Oppdragstaker)))

df <- data_softbottom %>%
  group_by(Vannlokalitet) %>%
  summarize(Oppdragstakere = paste(unique(Oppdragstaker), collapse = ", "),
            n_Oppdragstaker = length(unique(Oppdragstaker)))
nrow(df)
df

df <- data_softbottom %>%
  filter(Vannlokalitet %in% "Arendal, BT44") %>%
  group_by(Parameter_navn) %>%
  summarize(Oppdragstakere = paste(unique(Oppdragstaker), collapse = ", "),
            n_Oppdragstaker = length(unique(Oppdragstaker)))
nrow(df)
df %>% filter(n_Oppdragstaker > 1)

df <- data_softbottom %>%
  filter(Vannlokalitet %in% "Arendal, BT44" & Parameter_navn == "Individantall marin bløtbunnsfauna (takson) per arealenhet")
nrow(df)
df %>% filter(Year == 2015) %>% as.data.frame() %>% head()
xtabs(~VitenskapligNavn, df) %>% sort(decreasing = TRUE) %>% head()

# Species data are from NIVA, NIVA/HI plus Norconsult
data_softbottom %>%
  filter(Vannlokalitet %in% "Arendal, BT44" & VitenskapligNavn %in% "Abra nitida") %>%
  ggplot(aes(Year, Verdi, color = Oppdragstaker)) + geom_point()

# Diversity data are from NIVA and NIVA/HI only
data_softbottom %>%
  filter(Vannlokalitet %in% "Arendal, BT44" & Parameter_navn %in% "Hurlberts diversitetsindeks (ES100) marin bløtbunnsfauna") %>%
  ggplot(aes(Year, Verdi, color = Oppdragstaker)) + geom_point()



#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 4. Save data ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

saveRDS(data_hardbottom, file = "Data/04_vannmiljo_hardbottom.rds")
saveRDS(data_softbottom, file = "Data/04_vannmiljo_softbottom.rds")


