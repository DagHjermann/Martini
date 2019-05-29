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
# 1. Okokyst data ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

dat <- read_excel("Input_data/VannmiljoEksport_2019-05-27_okokyst.xlsx")

# Seems that dyp has been set to NA where Dyp is really zero
dat$Ovre_dyp[is.na(dat$Ovre_dyp)] <- 0
dat$Nedre_dyp[is.na(dat$Nedre_dyp)] <- 0

# Convert Verdi (value) to numbers (first have to change decimal comma to 'full stop')
x <- as.numeric(sub(",", ".", dat$Verdi, fixed = TRUE))
sum(is.na(x))  # 0, so all ar ok
dat$Verdi <- x

x <- as.numeric(sub(",", ".", dat$Ovre_dyp, fixed = TRUE))
sum(is.na(value))  # 0, so all ar ok
dat$Ovre_dyp <- x

x <- as.numeric(sub(",", ".", dat$Nedre_dyp, fixed = TRUE))
sum(is.na(value))  # 0, so all ar ok
dat$Nedre_dyp <- x

# Make time variable
dat$Time <- ymd_hms(dat$Tid_provetak)

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

map_df <- dat %>% 
  count(Vannlokalitet, Long, Lat) %>% 
  mutate(Popup = paste0(Vannlokalitet, "<br>Long,lat = ", round(Long, 4), ", ", round(Lat, 4), ""))

#
# Test plot
#
library(leaflet)
leaflet() %>%
  addTiles() %>%  # Default OpenStreetMap map tiles
  addMarkers(lng = map_df$Long, lat = map_df$Lat,
             popup = map_df$Popup)

map_df_sel <- map_df %>%
  filter(Long > 7.5 & Lat < 60)

leaflet() %>%
  addTiles() %>%  # Default OpenStreetMap map tiles
  addMarkers(lng = map_df_sel$Long, lat = map_df_sel$Lat,
             popup = map_df_sel$Popup)

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# .a A couple of tables ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

xtabs(~Vannlokalitet, dat) %>% length()  # 203
xtabs(~Vannlokalitet, dat) %>% sort() %>% rev() %>% names() %>% dput()

xtabs(~year(Time), dat)

# 1990-2012: only 'Maksdypindeks makroalger'
xtabs(~Parameter_navn, dat %>% filter(year(Time) <= 2008))   # Maksdypindeks makroalger, beskyttet kyst/fjord (MSMDI3) 
xtabs(~Vannlokalitet, dat %>% filter(year(Time) <= 2008))    # Tromøy nord, HT113 
xtabs(~Parameter_navn, dat %>% filter(year(Time) %in% 2009:2012))   # Maksdypindeks makroalger
xtabs(~Vannlokalitet, dat %>% filter(year(Time) %in% 2009:2012))    # 10 locations
# dat %>% filter(year(Time) == 1990) %>% View()

xtabs(~Parameter_navn, dat) %>% length()  # 88
xtabs(~Parameter_navn, dat) %>% sort() %>% rev() %>% .[1:20]
# 1. Relativ forekomst makroalger (skala)
# 2. Individantall planteplankton (takson) per volumenhet + Forekomst av arter/taksa planteplankton
# 3. Temperatur etc. 

xtabs(~Medium_navn, dat)
# Saltvann Sediment saltvann 
#   318159               278

xtabs(~Ovre_dyp, dat) %>% length()  # 1326

xtabs(~Oppdragstaker, dat)

# Number of stations per parameter
dat %>%
  count(Vannlokalitet, Parameter_navn) %>%
  count(Parameter_navn) %>%
  arrange(desc(n)) %>% View()

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# .b Add Stasjonstype ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Make file for manual edit
# Manually added "Stasjonstype" to this file and saved as "Parameter_types.csv"

if (FALSE){
dat %>%
  count(Vannlokalitet, Parameter_navn, Enhet) %>%
  count(Parameter_navn, Enhet) %>%
  arrange(desc(n)) %>% 
  write.csv2("Input_data/Parameter_type01.csv")  
}

#
# Add Parametertype
#
# dat <- dat %>% select(-Stasjonstype)
df_paramtypes <- read.csv2("Input_data/Parameter_types.csv")
dat <- dat %>%
  left_join(df_paramtypes %>% select(-n))

# Check "deviant" stations
tab <- xtabs(~Vannlokalitet + Parametertype, dat)
tab[apply(tab>0,1,sum) > 1, ]

#
# Add Stasjonstype
#
# Make file for manual edit
# Manually added "Stasjonstype" to this file and saved as "Parameter_types.csv"
if (FALSE){
  dat %>%
  ungroup() %>%
  count(Vannlokalitet, Parametertype, .drop = FALSE) %>%
  tidyr::spread(Parametertype, n) %>%
  write.csv2("Input_data/Station_type_01.csv")  
}

#
# Add Stasjonstype
#
df_stationtype <- read.csv2("Input_data/Station_types.csv")
dat <- dat %>%
  left_join(df_stationtype %>% select(Vannlokalitet, Stasjonstype))

xtabs(~Stasjonstype + Oppdragstaker, dat)

#
# a. Hardbunn ----
#
dat %>%
  filter(Stasjonstype %in% "Hardbunn" & year(Time) >= 2013) %>%
  count(Vannlokalitet, Parameter_navn) %>%
  count(Vannlokalitet) %>% View()

dat %>%
  filter(Stasjonstype %in% "Hardbunn" & year(Time) >= 2013) %>%
  count(Vannlokalitet, Parameter_navn) %>%
  count(Parameter_navn) %>% View()

dat %>%
  filter(Vannlokalitet %in% "Homborøya, HR105") %>%
  count(Time, Parameter_navn) %>%
  count(Time)
dat %>%
  filter(Vannlokalitet %in% "Homborøya, HR105") %>%
  filter(year(Time) == 2016) %>% # View()
  count(Parameter_navn)

dat %>%
  filter(Parameter_navn %in% c("Dekningsgrad makroalger (skala)", "Relativ forekomst makroalger (skala)") |
           grepl("Maksdypindeks makroalger", Parameter_navn)) %>%
  mutate(Param_short = substr(Parameter_navn, 1, 7)) %>%
  xtabs(~Vannlokalitet + Param_short, .)

# Only data from 0 m depth, example
dat %>%
  filter(Vannlokalitet %in% "Alstein") %>% View()

# Real hard-bottom data, example
dat %>%
  filter(Vannlokalitet %in% "Veslekalven, HT3") %>% View()

# Real hard-bottom data
dat %>%
  filter(Parameter_navn %in% c("Relativ forekomst makroalger (skala)") |
           grepl("Maksdypindeks makroalger", Parameter_navn)) %>%
  filter(year(Time) >= 2013)  %>%
  mutate(Param_short = substr(Parameter_navn, 1, 7)) %>%
  xtabs(~Vannlokalitet + year(Time), .)

#
# ..b Bløtbunn ----
#
dat %>%
  filter(Stasjonstype %in% c("Bløtbunn", "Bløtbunn, Hydrologi") & year(Time) >= 2013) %>%
  count(Vannlokalitet, Parameter_navn) %>%
  count(Vannlokalitet) %>% View("No of parameters")

dat %>%
  filter(Stasjonstype %in% c("Bløtbunn", "Bløtbunn, Hydrologi") & year(Time) >= 2013) %>%
  count(Vannlokalitet, Parameter_navn) %>%
  count(Parameter_navn) %>% View("No of stations")

dat %>%
  filter(grepl("diversitetsindeks", Parameter_navn) & year(Time) >= 2013) %>%
  xtabs(~Vannlokalitet + year(Time), .)



dat %>%
  filter(Stasjonstype %in% "Hardbunn" & 
           year(Time) >= 2012) %>%
  xtabs(~Vannlokalitet + year(Time), .)

dat %>%
  filter(Stasjonstype %in% "Hardbunn" & 
           year(Time) >= 2017) %>%
  xtabs(~Vannlokalitet + Oppdragstaker, .)

dat %>%
  filter(Stasjonstype %in% c("Bløtbunn", "Bløtbunn, Hydrologi") & 
           year(Time) >= 2017) %>%
  xtabs(~Vannlokalitet + year(Time), .)

dat %>%
  filter(Stasjonstype %in% c("Bløtbunn", "Bløtbunn, Hydrologi") & 
           year(Time) >= 2017) %>%
  xtabs(~Vannlokalitet + Oppdragstaker, .)


dat %>%
  filter(Stasjonstype %in% "Hardbunn" & 
           year(Time) >= 2013 & 
           Oppdragstaker %in% c("Havforskningssinstituttet","HI/NIVA")) %>%
  xtabs(~Vannlokalitet + year(Time), .)


dat %>%
  filter(Stasjonstype %in% "Hardbunn" & 
           year(Time) >= 2013 & 
           Oppdragstaker %in% c("Havforskningssinstituttet","HI/NIVA")) %>%
  xtabs(~Vannlokalitet + year(Time), .)

dat %>%
  filter(Stasjonstype %in% c("Bløtbunn", "Bløtbunn, Hydrologi") & 
           year(Time) >= 2013 & 
           Oppdragstaker %in% c("Havforskningssinstituttet","HI/NIVA")) %>%
  xtabs(~Vannlokalitet + year(Time), .)


dat %>%
  filter(Stasjonstype %in% "Hardbunn" & year(Time) >= 2013) %>%
  xtabs(~Vannlokalitet + year(Time), .)

dat %>%
  filter(Stasjonstype %in% "Hardbunn" & year(Time) >= 2013) %>%
  xtabs(~Vannlokalitet + year(Time), .)

  

# Error in VT10 (hydrography station)?
# Seems to be one single Totalnitrogen measured in g/kg (as bløtbunn) 2018-06-14
dat %>% filter(grepl("VT10", Vannlokalitet)) %>%
  xtabs(~Parameter_navn + Parametertype, .)
# dat %>% filter(grepl("VT10", Vannlokalitet)) %>%
#   filter(Parameter_navn == "Totalnitrogen") %>%
#   select(Time, Ovre_dyp, Nedre_dyp, Parameter_navn, Enhet, Verdi, Parametertype)  %>%
#   arrange(Time, Ovre_dyp, Nedre_dyp)  %>%
#   View()
tab1 <- dat %>% filter(grepl("VT10", Vannlokalitet)) %>%
  filter(Parameter_navn == "Totalnitrogen" & Enhet == "g/kg N t.v.") %>%
  xtabs(~Time + paste(Ovre_dyp, Nedre_dyp, sep = "-"), .)
tab2 <- dat %>% filter(grepl("VT10", Vannlokalitet)) %>%
  filter(Parameter_navn == "Totalnitrogen" & Enhet == "µg/l N") %>%
  xtabs(~Time + paste(Ovre_dyp, Nedre_dyp, sep = "-"), .)
tab1
# tab2

# Korsen - bløtbunn, but contains "Suspendert tørrstoff"?
dat %>% filter(grepl("BR113", Vannlokalitet)) %>%
  xtabs(~Parameter_navn + Parametertype + Vannlokalitet, .)

# Stations that are in fact both soft bottom and hydrology, which is OK (examples)
dat %>% filter(grepl("BT71", Vannlokalitet)) %>%
  xtabs(~Parameter_navn + Parametertype, .)
dat %>% filter(grepl("Broemsneset, BR114/VR52", Vannlokalitet)) %>%
  xtabs(~Parameter_navn + Parametertype, .)


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# .c A hydrography station ----
# From Økokyst klima 
# Both planteplankton and hydrography
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

dat %>%
  filter(grepl("VT52", Vannlokalitet)) %>% #  View()
  xtabs(~Parameter_navn, .)

dat %>%
  filter(grepl("VT52", Vannlokalitet) & Time %in% ymd_hms("2017-08-22 00:00:00")) %>% #  View()
  xtabs(~Parameter_navn, .)

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# .d Hydrography, plankton data ----
# Two parameters
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# 0-30 m depth
# Value is always 1
dat_occurence <- dat %>%
  filter(grepl("VT52", Vannlokalitet) & Time %in% ymd_hms("2017-08-22 00:00:00")) %>% #  View()
  filter(Parameter_navn %in% "Forekomst av arter/taksa planteplankton")
dat_occurence %>%
  select(VitenskapligNavn, Ovre_dyp, Nedre_dyp, Verdi, Enhet) %>% View()

# 5 m depth
# Value is some number
dat_density <- dat %>%
  filter(grepl("VT52", Vannlokalitet) & Time %in% ymd_hms("2017-08-22 00:00:00")) %>% #  View()
  filter(Parameter_navn %in% "Individantall planteplankton (takson) per volumenhet")
dat_density %>%
  select(VitenskapligNavn, Ovre_dyp, Nedre_dyp, Verdi, Enhet) %>% View()

# Overlap between the names almost (!) 100%
names1 <- dat_occurence %>% pull(VitenskapligNavn)
names2 <- dat_density %>%   pull(VitenskapligNavn)
mean(names1 %in% names2)  # 85 %
mean(names2 %in% names1)  # 100 %
names1[!names1 %in% names2]

# Plot the 30 most common algae
dat_density %>%
  top_n(30, Verdi) %>%
  ggplot(aes(VitenskapligNavn, Verdi)) + 
  geom_col() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0), plot.margin = margin(0.5, 4, 0.5, 0.5, "cm"))


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# .e Hydrography, physics/chemistry ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

dat_sel <- dat %>%
  filter(grepl("VT52", Vannlokalitet) & Time %in% ymd_hms("2017-08-22 00:00:00"))

dat_sel %>%
  group_by(Parameter_navn) %>%
  summarize(max(Nedre_dyp))

# dat_sel %>% filter(Parameter_navn %in% "Totalfosfor") %>% View()
# dat_sel %>% filter(Parameter_navn %in% "Oksygen") %>% View()

gg <- dat_sel %>% 
  filter(VitenskapligNavn %in% "Artsuavhengig") %>% # View()
  mutate(Dyp = (Ovre_dyp + Nedre_dyp)/2) %>%
  arrange(Parameter_navn, Dyp) %>%
  ggplot(aes(Verdi, -Dyp)) +
  geom_path() + geom_point() +
  facet_wrap(~Parameter_navn, scales = "free_x")
gg
gg + coord_cartesian(ylim = c(-31,0))



#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 2. Havforsuring data ----
# VannmiljoEksport_2019-05-27_havforsuring.xlsx
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

dat <- read_excel("Input_data/VannmiljoEksport_2019-05-27_havforsuring.xlsx")

# Seems that dyp has been set to NA where Dyp is really zero
dat$Ovre_dyp[is.na(dat$Ovre_dyp)] <- 0
dat$Nedre_dyp[is.na(dat$Nedre_dyp)] <- 0

# Convert Verdi (value) to numbers (first have to change decimal comma to 'full stop')
x <- as.numeric(sub(",", ".", dat$Verdi, fixed = TRUE))
sum(is.na(x))  # 0, so all ar ok
dat$Verdi <- x

x <- as.numeric(sub(",", ".", dat$Ovre_dyp, fixed = TRUE))
sum(is.na(value))  # 0, so all ar ok
dat$Ovre_dyp <- x

x <- as.numeric(sub(",", ".", dat$Nedre_dyp, fixed = TRUE))
sum(is.na(value))  # 0, so all ar ok
dat$Nedre_dyp <- x

# Make time variable
dat$Time <- ymd_hms(dat$Tid_provetak)

# Fix Siktedyp depth to zero
sel <- dat$Parameter_navn == "Siktedyp"; sum(sel); # table(dat$Nedre_dyp[sel])
dat$Nedre_dyp[sel] <- 0

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# .a A couple of tables ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

xtabs(~Vannlokalitet, dat) %>% length()  # 323

# All names
# xtabs(~Vannlokalitet, dat) %>% sort() %>% rev() %>% names() %>% dput()

xtabs(~year(Time), dat) # 2010 - 2016

xtabs(~Parameter_navn, dat) %>% length()  # 13

xtabs(~Parameter_navn, dat)

# Parameter_navn
# Fosfat (ufiltrert)             Løst reaktivt silikat       Løst uorganisk karbon (DIC) 
#               1738                              1736                              2523 
# Metningsgrad aragonitt (Omega Ar)   Metningsgrad kalsitt (Omega Ca)                            Nitrat 
#                              2480                              2481                              1562 
# Nitrat + nitritt  Partialtrykk karbondioksid (CO2)                                pH 
#             137                               289                               255 
# pH (totalskala)                         Salinitet                        Temperatur 
#            2474                              2558                              2564 
# Total alkalitet 
#            2502


xtabs(~Medium_navn, dat)  # Saltvann

xtabs(~Ovre_dyp, dat) %>% length()  # 446

xtabs(~Oppdragstaker, dat)  # note: both Havforskningsinstituttet and HI given

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# 3. Indre og Ytre Oslofjord ----
#    VannmiljoEksport_2019-05-27_indre_ytre_oslofjord.xlsx - data from "Oslofjorden" (selected i station part) 
#                                                            and "Tiltaksrettet overvåkning" 
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

dat_of <- read_excel("Input_data/VannmiljoEksport_2019-05-27_indre_ytre_oslofjord.xlsx", col_types = "text")

# Convert Verdi (value) to numbers (first have to change decimal comma to 'full stop')
x <- as.numeric(sub(",", ".", dat_of$Verdi, fixed = TRUE))
sum(is.na(x))  # 0, so all ar ok
dat_of$Verdi <- x

x <- as.numeric(sub(",", ".", dat_of$Ovre_dyp, fixed = TRUE))
sum(is.na(value))  # 0, so all ar ok
dat_of$Ovre_dyp <- x

x <- as.numeric(sub(",", ".", dat_of$Nedre_dyp, fixed = TRUE))
sum(is.na(value))  # 0, so all ar ok
dat_of$Nedre_dyp <- x

# Seems that dyp has been set to NA where Dyp is really zero
sel <- is.na(dat_of$Ovre_dyp); sum(sel)
dat_of$Ovre_dyp[sel] <- 0
sel <- is.na(dat_of$Nedre_dyp); sum(sel)
dat_of$Nedre_dyp[sel] <- 0

# Make time variable
dat_of$Time <- ymd_hms(dat_of$Tid_provetak)

# Fix Siktedyp depth to zero
sel <- dat_of$Parameter_navn == "Siktedyp"; sum(sel); # table(dat_of$Nedre_dyp[sel])
dat_of$Nedre_dyp[sel] <- 0

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# .a A couple of tables ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

xtabs(~Vannlokalitet, dat_of) %>% length()  # 411
xtabs(~Vannlokalitet, dat_of) %>% sort() %>% rev() %>% names() %>% dput()

# Lots of Bjørvika, Pipervika, Lohavn stations, but with little data (see below)
tab <- xtabs(~Vannlokalitet, dat_of)

x <- "Bjørvika"
sel <- grepl(x, names(tab)); sum(sel); mean(sel)
sel <- grepl(x, dat_of$Vannlokalitet); sum(sel); mean(sel)   # 31 % of stations, 1.4% of data

x <- "Pipervika"
sel <- grepl(x, names(tab)); sum(sel); mean(sel)
sel <- grepl(x, dat_of$Vannlokalitet); sum(sel); mean(sel)   # 10 % of stations, 0.4% of data

x <- "Lohavn"
sel <- grepl(x, names(tab)); sum(sel); mean(sel)
sel <- grepl(x, dat_of$Vannlokalitet); sum(sel); mean(sel)   # 6 % of stations, 0.2% of data

x <- "Bispevika"
sel <- grepl(x, names(tab)); sum(sel); mean(sel)
sel <- grepl(x, dat_of$Vannlokalitet); sum(sel); mean(sel)   # 2.6 % of stations, 0.1% of data

# Years
xtabs(~year(Time), dat_of)  # 1995 - 2018, but only N=2 in 2012...

xtabs(~Parameter_navn, dat_of) %>% length()  # 284!
xtabs(~Parameter_navn, dat_of) %>% sort() %>% rev() %>% .[1:40]  # includes contaminants
# 1. Relativ forekomst makroalger (skala)
# 2. Individantall planteplankton (takson) per volumenhet + Forekomst av arter/taksa planteplankton
# 3. Temperatur etc. 

xtabs(~Medium_navn, dat_of)    # inlcudes biota tissue, particles, sediment
# Saltvann Sediment saltvann 
#   318159               278

xtabs(~Ovre_dyp, dat_of) %>% length()  # 360

xtabs(~Oppdragstaker, dat_of)   # dnV, NGI, NIVA, Norconsult dominating


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# .b Add Stasjonstype ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

df_paramtypes <- read.csv2("Input_data/Parameter_types.csv")
dat_of <- dat_of %>%
  left_join(df_paramtypes %>% select(-n))

# Check "deviant" stations
tab <- xtabs(~Vannlokalitet + addNA(Parametertype), dat_of)
tab[apply(tab>0,1,sum) > 1, ]


# Check one  station with quite a lot of both Hydrography and Parametertype = NA (usually contaminants)
dat_of %>% filter(Vannlokalitet %in% "Indre Drammensfjorden (DH2)") %>%
  count(Parameter_navn, Parametertype) %>% View()

dat_of %>% filter(Vannlokalitet %in% "Indre Drammensfjorden (DH2)") %>%
  select(Time, Ovre_dyp, Nedre_dyp, Parameter_navn, Parametertype, Verdi, Enhet) %>% 
  arrange(Parameter_navn, Parametertype, Time, Ovre_dyp, Nedre_dyp, Verdi, Enhet) %>% 
  View()

# Tabulate hydrology stations (in order)
tab <- xtabs(~Vannlokalitet, dat_of %>% filter(Parametertype %in% "Hydrologi" & Medium_navn %in% "Saltvann"))
tab <- sort(tab, decreasing = TRUE)
tabdf <- data.frame(tab)
tabdf$Cumul_percent <- cumsum(tabdf$Freq)/sum(tabdf$Freq)*100
tabdf

tab2 <- dat_of %>% 
  filter(Vannlokalitet %in% tabdf$Vannlokalitet & year(Time) >= 2007) %>%
  xtabs(~Vannlokalitet + year(Time), .)

# Tabulate soft bottom stations (in order)
tab <- xtabs(~Vannlokalitet, dat_of %>% filter(Parametertype %in% "Bløtbunn" & Medium_navn %in% "Saltvann"))
tab <- sort(tab, decreasing = TRUE)
tabdf <- data.frame(tab)
tabdf$Cumul_percent <- cumsum(tabdf$Freq)/sum(tabdf$Freq)*100
tabdf

tab2 <- dat_of %>% 
  filter(Vannlokalitet %in% tabdf$Vannlokalitet) %>%
  xtabs(~Vannlokalitet + year(Time) + Oppdragstaker, .)
tab2



# ..a hard-bottom data ----
dat_of %>%
  filter(Parameter_navn %in% c("Relativ forekomst makroalger (skala)") |
           grepl("Maksdypindeks makroalger", Parameter_navn)) %>%
  filter(year(Time) >= 2013)  %>%
  mutate(Param_short = substr(Parameter_navn, 1, 7)) %>%
  xtabs(~Vannlokalitet + year(Time), .)

#
# ..b Bløtbunn ----
#


dat_of %>%
  filter(grepl("diversitetsindeks", Parameter_navn) & year(Time) >= 2013) %>%
  xtabs(~Vannlokalitet + year(Time), .)
