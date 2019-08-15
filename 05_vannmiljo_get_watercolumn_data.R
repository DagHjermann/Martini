#
# Collect data exported from Vannmiljø and export it
# 
# Actually we export all the data, not only water column data but also the soft + hard bottom data
#
# 1. VannmiljoEksport_2019-05-27_okokyst.xlsx              - all Økokyst data (NIVA and IMR)
# 2. VannmiljoEksport_2019-05-27_havforsuring.xlsx         - all Havforsuring data
# 3. VannmiljoEksport_2019-05-27_indre_ytre_oslofjord.xlsx - data from "Oslofjorden" (selected i station part) 
#                                                            and "Tiltaksrettet overvåkning" 
# Exported to K:
# "\\niva-of5\OSL-Data-NIVA\Avdeling\Vass\316_Miljøinformatikk\Prosjekter\180116_Martini\Datasets"

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
dat_acid <- read_excel("Input_data/VannmiljoEksport_2019-05-27_havforsuring.xlsx")
dat_of <- read_excel("Input_data/VannmiljoEksport_2019-05-27_indre_ytre_oslofjord.xlsx", col_types = "text")

# Combine
dat <- rbind(dat_okokyst, dat_acid, dat_of)

# Seems that dyp has been set to NA where Dyp is really zero
dat$Ovre_dyp[is.na(dat$Ovre_dyp)] <- 0
dat$Nedre_dyp[is.na(dat$Nedre_dyp)] <- 0

# Convert some character variables to numbers (first have to change decimal comma to 'full stop')

vars <- c("Verdi", "Ovre_dyp", "Nedre_dyp", "UTM33 Ost (X)", "UTM33 Nord (Y)")

for (var in vars){
  x <- as.numeric(sub(",", ".", dat[[var]], fixed = TRUE))
  cat(var, "- check that number is zero:", sum(is.na(x)), "\\n")  # 0, so all ar ok
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


write.csv(dat, "\\\\niva-of5\\OSL-Data-NIVA\\Avdeling\\Vass\\316_Miljøinformatikk\\Prosjekter\\180116_Martini\\Datasets\\Vannmiljo_data.csv")

# dat <- read.csv("\\\\niva-of5\\OSL-Data-NIVA\\Avdeling\\Vass\\316_Miljøinformatikk\\Prosjekter\\180116_Martini\\Datasets\\Vannmiljo_data.csv")

tab <- xtabs(~VitenskapligNavn, dat)
length(tab)
tab %>% sort(decreasing = TRUE) %>% .[1:30]
