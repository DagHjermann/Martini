#
# Adding 2017 data (from Vannmiljø, long format) to Phil's data (wide format)
# Also: the Vannmiljø data are in weight per L, Phil's data are in micromol per L
#
# dat = Original data set from Phil
# dat_add = data from Vannmiljø. 
#   NOTE 1: also contains 2016, which is used to check 
#   NOTE 2: "Nedre_dyp" = øvre dyp, which equals Depth_M in Phil's data
#

library(readxl)
library(dplyr)
library(ggplot2)
library(lubridate)

#
# Conversion to micromol, constants needed
#
# From
#   http://www.ices.dk/marine-data/tools/Pages/Unit-conversions.aspx
# 1 µg P/l = 1/MW P = 0.032285 µmol/l
# 1 µg N/l = 1/MW N = 0.071394 µmol/l
# 1 µg Si/l = 1/MW Si = 0.035606 µmol/l
# MW SiO3 = 76.083820 µg/l
# MW Si = 28.085530 µg/l
# MW C = 12.011 µg/l

convert_P <- 0.032285 
convert_N <- 0.071394  
convert_C <-  0.083257      # = 1/12.011
convert_SiO2 <- 0.0166444   # = 1/60.08   # from marelac::molweight("SiO2")/1000
convert_Si <- 0.035606      # = 1/28.08

# A sidetrack: Also check marelac package
marelac::molweight("P")/1000     # ok
marelac::molweight("N")/1000     # OK
marelac::molweight("SiO2")/1000  # 

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Data ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

dat <- read_excel("Input_data/Arendal_allvars_1990_2016c.xlsx", na = "NaN")
dat

dat_add <- read_excel("Input_data/VannmiljoEksport_vannreg.xlsx")
# Two of these warnings, they are just the two last lines, so nevermind:
# Expecting logical in AA4269 / R4269C27: got 'AUTO'
# tail(dat_add %>% select(Parameter_navn, Medium_navn))
# xtabs(~Medium_navn, dat_add)

# Seems that dyp has been set to NA where Dyp is really zero
dat_add$Ovre_dyp[is.na(dat_add$Ovre_dyp)] <- 0
dat_add$Nedre_dyp[is.na(dat_add$Nedre_dyp)] <- 0

# Convert Verdi (value) to numbers (first have to change decimal comma to 'full stop')
value <- as.numeric(sub(",", ".", dat_add$Verdi, fixed = TRUE))
sum(is.na(value))  # 0, so all ar ok
dat_add$Verdi <- value

# Make time variable
dat_add$Time <- ymd_hms(dat_add$Tid_provetak)

# Fix Siktedyp depth (was set to 23!!)
sel <- dat_add$Parameter_navn == "Siktedyp"; sum(sel); # table(dat_add$Nedre_dyp[sel])
dat_add$Nedre_dyp[sel] <- 0

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# A couple of tables ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

dat %>%
  select(Year, Temp_C:Secchi_m) %>%
  tidyr::gather("Parameter", "Value", -Year) %>%
  filter(Year >= 2015 & !is.na(Value)) %>%
  xtabs(~Parameter + Year, .)

xtabs(~Parameter_navn + year(Time), dat_add)
# Nitrat + nitritt, the three "Partikulært" variables: 2017 only
# Suspendert tørrstoff: 2016 only
# "Individantall planteplankton (takson) per volumenhet" has more records, but will not be included

# Units
dat_add %>%
  count(Parameter_navn, Enhet) %>% View("Unit")

# Depths, for Nitrat and fosfat
xtabs(~year(Time) + Nedre_dyp, dat_add %>% filter(Parameter_navn == "Nitrat"))
xtabs(~year(Time) + Nedre_dyp, dat_add %>% filter(Parameter_navn == "Fosfat (ufiltrert)"))
xtabs(~year(Time) + Nedre_dyp, dat_add %>% filter(Parameter_navn == "Løst reaktivt silikat"))
# No 2017 data for depth 50 and 75?

# Depths, for O2
xtabs(~year(Time) + Nedre_dyp, dat_add %>% filter(Parameter_navn == "Oksygen"))
# No 2017 data for depth 0, 5, 10, 20?


# Siktedyp (Sechhi depth) given "nedre dyp" = 23 (!)
# FIxed above
# xtabs(~year(Time) + Nedre_dyp, dat_add %>% filter(Parameter_navn == "Siktedyp"))
# dat_add %>% filter(Parameter_navn == "Siktedyp") %>% select(Time, Parameter_navn, Verdi, Enhet)

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Define parameters ----
#
# NOTE: very uncertain how or whether to convert the POP/PON/POC data.
# Now they are NOT transformed at all (original numbers used)
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# xtabs(~Parameter_navn, dat_add)
# xtabs(~Parameter_navn, dat_add) %>% names() %>% dput()
# colnames(dat) %>% dput()

parnames <- c("Temp_C" = "Temperatur",
              "Salinity_PSU" = "Salinitet",
              "Density_SigmaT" = "Tetthet",
              "O2_mlL" = "Oksygen",
              "O2sat_percent" = "Oksygenmetning",
              "PO4_umolL" = "Fosfat (ufiltrert)",
              "NO2_umolL" = "Nitritt", 
              "NO3_umolL" = "Nitrat",
              "NO2_NO3_umolL" = "Nitrat + nitritt",
              "NH4_umolL" = "Ammonium",
              "SiO4_umolL" = "Løst reaktivt silikat",
              # "N_P" = "",            # N_P = (NO2 + NO3)/PO4 - will be added later
              "Chla_ugL" = "Klorofyll a",
              "TotP_umolL" = "Totalfosfor", 
              "TotN_umolL" = "Totalnitrogen", 
              "POP_umolL" = "Partikulært organisk fosfor (POP)", 
              "PON_umolL" = "Partikulært organisk nitrogen (PON)", 
              "POC_umolL" = "Partikulært organisk karbon (POC)",  
              "TSM_mgL" = "Suspendert tørrstoff",
              "Secchi_m" = "Siktedyp")

conversion <- c("Temp_C" = 1,
                "Salinity_PSU" = 1,
                "Density_SigmaT" = 1,
                "O2_mlL" = 1,
                "O2sat_percent" = 1,
                "PO4_umolL" = convert_P,
                "NO2_umolL" = convert_N, 
                "NO3_umolL" = convert_N,
                "NO2_NO3_umolL" = convert_N,
                "NH4_umolL" = convert_N,
                "SiO4_umolL" = convert_Si,
                # "N_P" = "",
                "Chla_ugL" = 1,
                "TotP_umolL" = convert_P, 
                "TotN_umolL" = convert_N, 
                "POP_umolL" = 1, # convert_P,   # see above!
                "PON_umolL" = 1, # convert_N,   # see above!  
                "POC_umolL" = 1, # convert_C,   # see above!
                "TSM_mgL" = 1,
                "Secchi_m" = 1)


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Conversion to micromol, test ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

df1 <- dat %>% 
  filter(Year == 2016 & Month == 12 & Day == 6 & Depth_m == 0)
df2 <- dat_add %>% 
  filter(Time == ymd_hms("2016-12-06 00:00:00") & Ovre_dyp == 0) %>%
  select(Vannlokalitet, Time, Parameter_id, Parameter_navn, Ovre_dyp, Nedre_dyp, Verdi, Enhet)
# df2 %>% View()

# Show
df2 %>% filter(Parameter_navn %in% c("Fosfat (ufiltrert)", "Nitrat","Løst reaktivt silikat"))

# P substance example (given as weight of P) - confirms that conversion works
x <- df2 %>% filter(Parameter_navn == "Fosfat (ufiltrert)") %>% pull(Verdi)
x*convert_P
df1$PO4_umolL
# The other way round
df1$PO4_umolL*30.973762  # 30.97 = MW (molecular weight) of P
x

# N substance example (given as weight of N) - confirms that conversion works
x <- df2 %>% filter(Parameter_navn == "Nitrat") %>% pull(Verdi)
x*convert_N
df1$NO3_umolL
# The other way round
df1$NO3_umolL*14.006720   # MW of N = 1/convert_N
x

# SiO4 (given as weight of SiO2) - can make conversion work, see below
x <- df2 %>% filter(Parameter_navn == "Løst reaktivt silikat") %>% pull(Verdi)
x*convert_SiO2    # 1.45 - should fit, but doesn't
x*convert_Si      # 3.09 - fits
df1$SiO4_umolL    # 3.10
# The other way round
df1$SiO4_umolL*28.085530    # 87.21 - using the molecular weight of Si (not SiO2) - fits what's inside Vannmiljø
x
# So the conversion to Vannmiljø is wrong (either using the wrong constant or setting the wrong unit)

# TSM - OK
df2 %>% filter(Parameter_navn == "Suspendert tørrstoff")   # unit = mg/l
x <- df2 %>% filter(Parameter_navn == "Suspendert tørrstoff") %>% pull(Verdi)
x
df1$TSM_mgL

# Chl A - OK 
df2 %>% filter(Parameter_navn == "Klorofyll a")   # unit = ug/l
x <- df2 %>% filter(Parameter_navn == "Klorofyll a") %>% pull(Verdi)
x
df1$Chla_ugL




#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Make df2 (= dat_add with converted names/units)  ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Function for makiung df2
#
# Input: named string, e.g. c("Temp_C" = "Temperatur")
# Output: Data frame suitable for adding to data using left_join
#   Data frame has columns Time, Nedre_dyp and converted value (for chemistry: mol) of variable (with new name)
#
make_newpar <- function(par){
  # Parameter to add from Data 2 (Vannmiljø)
  df2_newpar <- dat_add %>%
    filter(Parameter_navn %in% par) %>%
    select(Time, Nedre_dyp, Verdi, Enhet) %>%
    mutate(Value = Verdi*conversion[names(par)]) %>%
    select(Time, Nedre_dyp, Value)
  colnames(df2_newpar)[3] <- names(par)
  df2_newpar
}

make_newpar(c("Temp_C" = "Temperatur")) %>% head()
make_newpar(parnames[9]) %>% head()
make_newpar(parnames[9]) %>% filter(Nedre_dyp == 0) %>% ggplot(aes(Time, NO2_NO3_umolL)) + geom_line()


# Data 1 (from Phil) - just for plotting (df_plot below)
df1 <- dat
df1$Time <- with(df1, ymd_hms(paste(Year, Month, Day, "00:00:00")))

# Data 2 (Vannmiljø) - "start data" without data columns
# Contains all comnbinations of Time and Nedre_dyp for the 'parnames' variables
df2 <- dat_add %>%
  filter(Parameter_navn %in% parnames) %>%
  count(Time, Nedre_dyp) %>%
  mutate(Year = year(Time),
         Month = month(Time),
         Day = mday(Time)) %>%
  select(-n)                   # We don't need this
nrow(df2)
xtabs(~year(Time), df2)

# Add parameter columns to "start data" - loops through all variables
for (i in 1:length(parnames)){
  df2_newpar <- make_newpar(parnames[i])
  df2 <- df2 %>% left_join(df2_newpar, by = c("Time", "Nedre_dyp"))
  }

head(df2)

#
# Fill out NO2_NO3_umolL
#
sel <- is.na(df2$NO2_NO3_umolL); sum(sel)
df2$NO2_NO3_umolL[sel] <- df2$NO2_umolL[sel] + df2$NO3_umolL[sel]


#
# Make N_P
#
df2$N_P <- df2$NO2_NO3_umolL/df2$PO4_umolL

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Compare using plots (2016 only) ----
#
# Set1 = Excel sheet from Phil
# Set2 = Vannmiljø data
#
# Cannot compare POP, POC, PON (not given for 2016 in Vannmiljø)
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

#
# Check one by one
#
df_plot <- bind_rows(
    df1 %>% mutate(Set = "set1"),
    df2 %>% mutate(Set = "set2") %>% rename(Depth_m = Nedre_dyp)
  ) %>%
  filter(year(Time) == 2016) # & Depth_m == 50)
df_plot

ggplot(df_plot, aes(Time, Temp_C, color = Set)) +
  geom_line() + facet_grid(Depth_m~Set)
ggplot(df_plot, aes(Time, NO2_NO3_umolL, color = Set)) +
  geom_line() + facet_grid(Depth_m~Set)
ggplot(df2, aes(Time, NO2_NO3_umolL, color)) +
  geom_line() + facet_grid(Nedre_dyp~.)

#
# Check 5 by 5
#
df_plot2 <- df_plot %>%
  tidyr::gather("Param", "Value", Temp_C:Secchi_m)

variables <- colnames(dat %>% select(Temp_C:Secchi_m))
length(variables)  # 20

# (Compare plots below with this table:)
xtabs(~Parameter_navn + year(Time), dat_add)

df_plot2 %>% 
  filter(Depth_m == 0 & Param %in% variables[1:5]) %>%
  ggplot(aes(Time, Value, color = Set)) +
  geom_line() + facet_grid(Param~Set, scales = "free")
# OK

df_plot2 %>% 
  filter(Depth_m == 0 & Param %in% variables[6:10]) %>%
  ggplot(aes(Time, Value, color = Set)) +
  geom_line() + facet_grid(Param~Set, scales = "free")
# OK, including the calculated NO2_NO3 for set2

df_plot2 %>% 
  filter(Depth_m == 0 & Param %in% variables[11:15]) %>%
  ggplot(aes(Time, Value, color = Set)) +
  geom_line() + facet_grid(Param~Set, scales = "free")
# OK, including the calculated N_P for Set2

df_plot2 %>% 
  filter(Depth_m == 0 & Param %in% variables[16:20]) %>%
  ggplot(aes(Time, Value, color = Set)) +
  geom_line() + facet_grid(Param~Set, scales = "free")
# POP, POC, PON not given for 2016


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Make dat_updated
#  I.e. add df2 (2017 data only) to dat ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

dat_updated <- bind_rows(
   dat,
   df2 %>% 
     mutate(Year = year(Time), Month = month(Time), Day = mday(Time)) %>%
     filter(Year == 2017 & Nedre_dyp %in% c(0,5,10,20,30,50,75)) %>% 
     rename(Depth_m = Nedre_dyp)
  ) %>%
  select(-Time)


dat_updated
dat_updated %>% tail(4)

# dat_updated %>% as.data.frame() %>% tail(1)

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Check time series visually ----
#
# We plot 2014-2017 onwards only
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

df_plot3 <- dat_updated %>%
  mutate(Time = ymd_hms(paste(Year, Month, Day, "00:00:00"))) %>%
  tidyr::gather("Param", "Value", Temp_C:Secchi_m)

# Variables to loop over
variables <- colnames(dat_updated %>% select(Temp_C:Secchi_m))
length(variables)  # 20

# Save a bunch of plots
for (i in 1:length(variables)){
# for (i in 16:19){
  if (i %in% 16:19){
    plottitle <- sprintf("%02i: %s (NOTE: 'mol/L' numbers until 2016, 'mg/L' numbers in 2017!!)", i,  variables[i])
  } else {
    plottitle <- sprintf("%02i: %s", i,  variables[i])
  }
  gg <- df_plot3 %>% 
    filter(Year >= 2014 & Param %in% variables[i]) %>%
    ggplot(aes(Time, Value, color = Depth_m)) +
    geom_line() + 
    facet_wrap(~Depth_m) +
    ggtitle(plottitle)
  ggsave(sprintf("Figures/01/%02i_%s.png", i,  variables[i]), gg)
  }

# i <- 18
# variables[i]
# df_plot3 %>% filter(Year >= 2014 & Param %in% variables[i]) %>%
#   group_by(Depth_m, Year) %>%
#   summarize(mean(Value, na.rm = TRUE)) %>% View()
# 
# dat_add %>% filter(year(Time) >= 2014 & Parameter_navn %in% "Partikulært organisk karbon (POC)") %>%
#   group_by(Nedre_dyp, year(Time)) %>%
#   summarize(mean(Verdi, na.rm = TRUE)) %>% View()



#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Save ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

openxlsx::write.xlsx(dat_updated, "Data/01_Arendal_allvars_1990_2017.xlsx")

