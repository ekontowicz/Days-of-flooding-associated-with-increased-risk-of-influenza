######################################################
# Code used for data wrangling and analysis
# The goal of this research was to estimate the 
# association between flooding and influenza
# when accounting for other known population and
# environmental determinants of influenza transmission
######################################################

#Libraries
library(tidyverse)
library(R2OpenBUGS)
library(maps)
library(maptools)
library(sf)
library(spdep)
library(rgdal)
library(coda)

#### FUNCTIONS ####
#Creates the Study year index value tha truns from May of the first year to April of the following year
# Used to ensure proper temporal alignment of flooding and influenza seasons.
study_year_cat <- function(dateVal, minYr = 2005){
  dateVal <- as.Date(dateVal)
  if ( length(dateVal) != 1 || class(dateVal) != "Date"){
    stop("That's not what this function is for")
  }
  for (yr in minYr:(minYr + 25)){
    if (dateVal >= as.Date(paste0(yr, "-05-01")) && 
        dateVal <= as.Date(paste0(yr+1, "-04-30"))){
      return(yr - minYr + 1)
    }
  }
  stop("No matching year found for date:", dateVal)
}

#### DATA WRANLING INFLUENZA DIAGNOSIS ####
#Private insurance database
diagcodes <- read.csv("FluDiagnosisCodes_allinfo.csv")

#### CLEANING DIAGCODES ####
range(diagcodes$date)

#Some missing and mislabeled 3 digit ZCTAs in the raw data sheet. Correcting those errors
diagcodes <- diagcodes %>% 
  mutate(zip3 = zip, 
         date = lubridate::ymd(paste(year, month, "01"))) %>% 
  filter("2005-06-01" <= date & date < "2017-01-01") %>% 
  arrange(zip)

diagcodes$zip3[diagcodes$zip3 %in% 23] <- 527
diagcodes$zip3[diagcodes$zip3 %in% 52] <- 522
diagcodes$zip3[diagcodes$zip3 %in% 81] <- 525
diagcodes$zip3[diagcodes$zip3 %in% 128] <- 582
diagcodes$zip3[diagcodes$zip3 %in% 201] <- 501
diagcodes$zip3[diagcodes$zip3 %in% 206] <- 506
diagcodes$zip3[diagcodes$zip3 %in% 224] <- 524
diagcodes$zip3[diagcodes$zip3 %in% 226] <- 526
diagcodes$zip3[diagcodes$zip3 %in% 252] <- 525
diagcodes$zip3[diagcodes$zip3 %in% 257] <- 525
diagcodes$zip3[diagcodes$zip3 %in% 424] <- 524
diagcodes$zip3[diagcodes$zip3 %in% 529] <- 527
diagcodes$zip3[diagcodes$zip3 %in% 530] <- 503
diagcodes$zip3[diagcodes$zip3 %in% 532] <- 523
diagcodes$zip3[diagcodes$zip3 %in% 541] <- 514
diagcodes$zip3[diagcodes$zip3 %in% 550] <- 503
diagcodes$zip3[diagcodes$zip3 %in% 554] <- 500
diagcodes$zip3[diagcodes$zip3 %in% 560] <- 506
diagcodes$zip3[diagcodes$zip3 %in% 570] <- 523
diagcodes$zip3[diagcodes$zip3 %in% 571] <- 511
diagcodes$zip3[diagcodes$zip3 %in% 572] <- 504
diagcodes$zip3[diagcodes$zip3 %in% 575] <- 515
diagcodes$zip3[diagcodes$zip3 %in% 582] <- 502
diagcodes$zip3[diagcodes$zip3 %in% 601] <- 502
diagcodes$zip3[diagcodes$zip3 %in% 614] <- 514
diagcodes$zip3[diagcodes$zip3 %in% 618] <- 526
diagcodes$zip3[diagcodes$zip3 %in% 622] <- 522
diagcodes$zip3[diagcodes$zip3 %in% 681] <- 515
diagcodes$zip3[diagcodes$zip3 %in% 712] <- 500
diagcodes$zip3[diagcodes$zip3 %in% 824] <- 524
diagcodes$zip3[diagcodes$zip3 %in% 982] <- 501

#These are not ZCTAs in Iowa so set to NA for later removal 
diagcodes$zip3[diagcodes$zip3 %in% c(24, 170, 333, 386, 495, 577, 605, 612, 630, 641, 744, 833)] <- NA


range(diagcodes$date)

#Adding study year index
diagcodes$study_year_idx <- sapply(diagcodes$date, study_year_cat)

#Reaggregrating to 3 digit ZCTA and study year index
diagcodes_3agg <- diagcodes %>% 
  group_by(study_year_idx, zip3) %>% 
  summarise(total_flu_diag = sum(flu_total, na.rm = T),
            total_ILI = sum(ILIness_total, na.rm = T),
            total_vacc = sum(vacc_total, na.rm = T),
            total_AsmAtck = sum(AsmAtck_total, na.rm = T),
            total_Asma = sum(Asma_total, na.rm = T),
            total_pat = sum(pat_vis_total, na.rm = T),
            total_live = sum(mem_live_total, na.rm = T))

#Removing thos spatials units not in Iowa 
diagcodes_3agg2 <- subset(diagcodes_3agg, zip3 != "NA")
diagcodes_3agg2_509 <- diagcodes_3agg2[diagcodes_3agg2$zip3 == "509",]
diagcodes_3agg3 <- anti_join(diagcodes_3agg2, diagcodes_3agg2_509)
sum(is.na(diagcodes_3agg3)) # 0 missing

#Hygenic Lab aggregrated data
flu_lab <- read_csv(file = "SHLData_3ZipAgg.csv")%>% 
  mutate(study_year_idx = stdy_yr_idx) %>% select(-stdy_yr_idx)
sum(is.na(flu_lab)) # 0 missing

flu_all <- diagcodes_3agg3 %>% group_by(zip3, study_year_idx) %>% 
  left_join(., flu_lab, by = c("zip3", "study_year_idx")) %>% ungroup()
sum(is.na(flu_all)) # 16 missing data points
#This comes from times were SHL did not have a test to report. setting these NAs to 0
sum(is.na(flu_all$total_labpos))
sum(is.na(flu_all$total_labtest))

flu_all[is.na(flu_all)] <- 0
sum(is.na(flu_all))  #0 missing data now

flu_all <- flu_all %>% mutate(flu_diag_lab = total_flu_diag + total_labpos,
                              ILI_diag_lab = total_ILI + total_labpos) %>% 
  select(-total_pat, -total_live)
write_csv(flu_all, path = "FluDiag&SHL_3ZCTA_YrIdx_Aggregated.csv")

#### DATA WRANGLING WEATHER ####
SynopMetar_MonthAgg <- read_csv(file = "SynopMetar_Combo_5ZCTAMonth_Agg_long.csv") %>% 
  select(Zip_Clean, stdy_year_idx, season, date,
         avg_AirTemp, avg_DewPointTemp, avg_RelHumid,
         avg_AblHumid_Kgm3, avg_AblHumid_gm3) %>% 
  mutate(avg_AirTemp.F = ((avg_AirTemp * 9/5) + 32),   # convert C to F
         avg_DewPointTemp.F = ((avg_DewPointTemp * 9/5) + 32)) %>% 
  select(-avg_AirTemp, -avg_DewPointTemp) %>% #dump the old measures so I can rename the new ones to the same name. will making the code below easier so i don't have to change everything
  mutate(avg_AirTemp = avg_AirTemp.F,
         avg_DewPointTemp = avg_DewPointTemp.F) %>% 
  select(-avg_AirTemp.F, -avg_DewPointTemp.F)

SynopMetar_MonthAgg <- SynopMetar_MonthAgg %>% select(Zip_Clean, stdy_year_idx, season, date, avg_AirTemp, avg_RelHumid, avg_AblHumid_gm3)
SynopMetar_MonthAgg$season <- NA
SynopMetar_MonthAgg$season[lubridate::month(SynopMetar_MonthAgg$date) %in% c(5, 6, 7)] <- "flood"
SynopMetar_MonthAgg$season[lubridate::month(SynopMetar_MonthAgg$date) %in% c(10, 11, 12, 1, 2, 3)] <- "flu"
SynopMetar_MonthAgg$season[lubridate::month(SynopMetar_MonthAgg$date) %in% c(8, 9)] <- "Aug&Sept"
SynopMetar_MonthAgg$season[lubridate::month(SynopMetar_MonthAgg$date) %in% c(4)] <- "April"

SynopMetar_zip3IdxAgg <- SynopMetar_MonthAgg %>%
  mutate(zip3 = as.numeric(substr(Zip_Clean, 1, 3)),
         study_year_idx = stdy_year_idx) %>% 
  select(-stdy_year_idx) %>% 
  group_by(zip3, study_year_idx, season) %>% 
  summarize(avg_TAVG = mean(avg_AirTemp),
            sd_TAVG = sd(avg_AirTemp),
            avg_RHum = mean(avg_RelHumid),
            sd_RHum = sd(avg_RelHumid),
            avg_AHum = mean(avg_AblHumid_gm3),
            sd_AHum = sd(avg_AblHumid_gm3))
#there is missing data. this is coming from sd columns. Happening becuase only one vaule for that zip3 and season, thus can't have an sd. going to set to 0
SynopMetar_zip3IdxAgg[is.na(SynopMetar_zip3IdxAgg)] <- 0


#Other station measures
Weather_MonthAgg <- read_csv(file = "Weather_5ZCTA_MonthAverages.csv") %>% 
  select(-season) #Want to get rid of old season data
Weather_MonthAgg$season <- NA
Weather_MonthAgg$season[lubridate::month(Weather_MonthAgg$date) %in% c(5, 6, 7)] <- "flood"
Weather_MonthAgg$season[lubridate::month(Weather_MonthAgg$date) %in% c(10, 11, 12, 1, 2, 3)] <- "flu"
Weather_MonthAgg$season[lubridate::month(Weather_MonthAgg$date) %in% c(8, 9)] <- "Aug&Sept"
Weather_MonthAgg$season[lubridate::month(Weather_MonthAgg$date) %in% c(4)] <- "April"

#Aggregrate up to the 3zip, study_year_idx and season level
Weather_MonthAgg2 <- Weather_MonthAgg %>% 
  group_by(zip3, study_year_idx, season) %>% 
  summarise(avg_TMAX = mean(TMAX.F, na.rm = T),
            sd_TMAX = sd(TMAX.F, na.rm = T),
            avg_TMIN = mean(TMIN.F, na.rm = T),
            sd_TMIN = sd(TMIN.F, na.rm = T),
            avg_PRCP = mean(PRCP, na.rm = T),
            sd_PRCP = sd(PRCP, na.rm = T),
            total_PRCP = sum(PRCP, na.rm = T),
            avg_SNOW = mean(SNOW, na.rm = T),
            sd_SNOW = sd(SNOW, na.rm = T),
            total_SNOW = sum(SNOW, na.rm = T))
#there is missing data. this is coming from sd columns. Happening becuase only one vaule for that zip3 and season, thus can't have an sd. going to set to 0
Weather_MonthAgg2[is.na(Weather_MonthAgg2)] <- 0

#Both weather datasets now are aggregrated to the 3zip study year season level. will now merge and then make wide
AllWeather <- SynopMetar_zip3IdxAgg %>% 
  left_join(., Weather_MonthAgg2, by = c("zip3", "study_year_idx", "season"))

#Saving the long data set above
# write_csv(AllWeather, path = "AllWeather_3zipYearFloodFluSeasonAgg_long.csv")

#Now making the long data set into wide
AllWeather <- read_csv("AllWeather_3zipYearFloodFluSeasonAgg_long.csv")


wide_avgPRCP <- AllWeather %>% 
  select(zip3, study_year_idx, avg_PRCP, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = avg_PRCP)
colnames(wide_avgPRCP)[c(-1, -2)] <- paste0(colnames(wide_avgPRCP)[c(-1, -2)], ".avgPRCP")

wide_sdPRCP <- AllWeather %>% 
  select(zip3, study_year_idx, sd_PRCP, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = sd_PRCP)
colnames(wide_sdPRCP)[c(-1, -2)] <- paste0(colnames(wide_sdPRCP)[c(-1, -2)], ".sdPRCP")

wide_totalPRCP <- AllWeather %>% 
  select(zip3, study_year_idx, total_PRCP, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = total_PRCP)
colnames(wide_totalPRCP)[c(-1, -2)] <- paste0(colnames(wide_totalPRCP)[c(-1, -2)], ".totalPRCP")

wide_avgSNOW <- AllWeather %>% 
  select(zip3, study_year_idx, avg_SNOW, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = avg_SNOW)
colnames(wide_avgSNOW)[c(-1, -2)] <- paste0(colnames(wide_avgSNOW)[c(-1, -2)], ".avgSNOW")

wide_sdSNOW <- AllWeather %>% 
  select(zip3, study_year_idx, sd_SNOW, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = sd_SNOW)
colnames(wide_sdSNOW)[c(-1, -2)] <- paste0(colnames(wide_sdSNOW)[c(-1, -2)], ".sdSNOW")

wide_totalSNOW <- AllWeather %>% 
  select(zip3, study_year_idx, total_SNOW, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = total_SNOW)
colnames(wide_totalSNOW)[c(-1, -2)] <- paste0(colnames(wide_totalSNOW)[c(-1, -2)], ".totalSNOW")

wide_avgTMAX <- AllWeather %>% 
  select(zip3, study_year_idx, avg_TMAX, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = avg_TMAX)
colnames(wide_avgTMAX)[c(-1, -2)] <- paste0(colnames(wide_avgTMAX)[c(-1, -2)], ".avgTMAX")

wide_sdTMAX <- AllWeather %>% 
  select(zip3, study_year_idx, sd_TMAX, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = sd_TMAX)
colnames(wide_sdTMAX)[c(-1, -2)] <- paste0(colnames(wide_sdTMAX)[c(-1, -2)], ".sdTMAX")

wide_avgTMIN <- AllWeather %>% 
  select(zip3, study_year_idx, avg_TMIN, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = avg_TMIN)
colnames(wide_avgTMIN)[c(-1, -2)] <- paste0(colnames(wide_avgTMIN)[c(-1, -2)], ".avgTMIN")

wide_sdTMIN <- AllWeather %>% 
  select(zip3, study_year_idx, sd_TMIN, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = sd_TMIN)
colnames(wide_sdTMIN)[c(-1, -2)] <- paste0(colnames(wide_sdTMIN)[c(-1, -2)], ".sdTMIN")

wide_avgTAVG <- AllWeather %>% 
  select(zip3, study_year_idx, avg_TAVG, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = avg_TAVG)
colnames(wide_avgTAVG)[c(-1, -2)] <- paste0(colnames(wide_avgTAVG)[c(-1, -2)], ".avgTAVG")

wide_sdTAVG <- AllWeather %>% 
  select(zip3, study_year_idx, sd_TAVG, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = sd_TAVG)
colnames(wide_sdTAVG)[c(-1, -2)] <- paste0(colnames(wide_sdTAVG)[c(-1, -2)], ".sdTAVG")

wide_avgRhum <- AllWeather %>% 
  select(zip3, study_year_idx, avg_RHum, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = avg_RHum)
colnames(wide_avgRhum)[c(-1, -2)] <- paste0(colnames(wide_avgRhum)[c(-1, -2)], ".avgRhum")

wide_sdRhum <- AllWeather %>% 
  select(zip3, study_year_idx, sd_RHum, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = sd_RHum)
colnames(wide_sdRhum)[c(-1, -2)] <- paste0(colnames(wide_sdRhum)[c(-1, -2)], ".sdRhum")

wide_avgAhum <- AllWeather %>% 
  select(zip3, study_year_idx, avg_AHum, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = avg_AHum)
colnames(wide_avgAhum)[c(-1, -2)] <- paste0(colnames(wide_avgAhum)[c(-1, -2)], ".avgAhum")

wide_sdAhum <- AllWeather %>% 
  select(zip3, study_year_idx, sd_AHum, season) %>% 
  group_by(zip3, study_year_idx) %>% 
  spread(key = season, value = sd_AHum)
colnames(wide_sdAhum)[c(-1, -2)] <- paste0(colnames(wide_sdAhum)[c(-1, -2)], ".sdAhum")

weather_wide <- wide_avgPRCP %>% 
  left_join(., wide_sdPRCP, by = c("zip3", "study_year_idx")) %>% 
  left_join(., wide_totalPRCP, by = c("zip3", "study_year_idx")) %>% 
  left_join(., wide_avgSNOW, by = c("zip3", "study_year_idx")) %>% 
  left_join(., wide_sdSNOW, by = c("zip3", "study_year_idx")) %>% 
  left_join(., wide_totalSNOW, by = c("zip3", "study_year_idx")) %>% 
  left_join(., wide_avgTMAX, by = c("zip3", "study_year_idx")) %>% 
  left_join(., wide_sdTMAX, by = c("zip3", "study_year_idx")) %>% 
  left_join(., wide_avgTMIN, by = c("zip3", "study_year_idx")) %>% 
  left_join(., wide_sdTMIN, by = c("zip3", "study_year_idx")) %>% 
  left_join(., wide_avgTAVG, by = c("zip3", "study_year_idx")) %>% 
  left_join(., wide_sdTAVG, by = c("zip3", "study_year_idx")) %>% 
  left_join(., wide_avgRhum, by = c("zip3", "study_year_idx")) %>% 
  left_join(., wide_sdRhum, by = c("zip3", "study_year_idx")) %>% 
  left_join(., wide_avgAhum, by = c("zip3", "study_year_idx")) %>% 
  left_join(., wide_sdAhum, by = c("zip3", "study_year_idx"))
# write_csv(weather_wide, path = "AllWeather_3zipStudySeasonIdxAgg.csv")

#### DATA WRANGLING WEATHER ####
StreamData <- readxl::read_xlsx("AllStream_Short.xlsx")
#Some gauges that formly did not have flood markers. updating those below
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5388310"] <- 10
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5464315"] <- 15
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5463050"] <- 89
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5463500"] <- 14
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5464780"] <- 16
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "6605750"] <- 90
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5420500"] <- 17
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5422000"] <- 11
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5484500"] <- 17
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5416900"] <- 14
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5453100"] <- 15
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5418400"] <- 12
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5418500"] <- 24
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5454220"] <- 10
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5464420"] <- 12.5
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5488500"] <- 14
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5489000"] <- 25
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "6602400"] <- 23
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "6809500"] <- 18
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "6811800"] <- 11.5
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "6600030"] <- 87
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "6811875"] <- 92.5
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "6817300"] <- 26.5
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5422600"] <- 12
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "6483495"] <- 16
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5471000"] <- 21.5
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5451770"] <- 12.5
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "5481300"] <- 19
StreamData$Flood_Stage_Height_Decimal[StreamData$Station_ID == "6602190"] <- 80

StreamZip <- read_csv(file = "StreamGauge_ZCTALxns.csv")[-2:-3]

flood <- StreamData %>% group_by(Station_ID) %>%  left_join(., StreamZip) %>%
  mutate(AvgAbove = ifelse(dayavg > Flood_Stage_Height_Decimal, yes = 1, no = 0),
         MaxAbove = ifelse(maxhght > Flood_Stage_Height_Decimal, yes = 1, no = 0),
         MinAbove = ifelse(minhght > Flood_Stage_Height_Decimal, yes = 1, no = 0)) %>% 
  select(date, Station_ID, Zip_Clean, AvgAbove, MaxAbove, MinAbove) %>% 
  mutate(MonthYear = lubridate::ymd(paste(lubridate::year(date), lubridate::month(date), "O1"))) %>% 
  arrange(Station_ID, date)


max_overlap_length <- function(sequence){
  is_flooded = 0
  overlap_ct = 0
  max_overlap_ct = 0
  for (i in 1:length(sequence)){
    # Case 1: Already flooding
    if (sequence[i] == 1 && is_flooded == 1){
      overlap_ct <- overlap_ct + 1
    }
    # Case 2: New flood event
    if (sequence[i] == 1 && is_flooded == 0){
      is_flooded <- 1
    }
    # Case 3: Newly ended flood
    if (sequence[i] == 0 && is_flooded == 1){
      is_flooded = 0
      max_overlap_ct = max(max_overlap_ct, overlap_ct)
      overlap_ct <- 0
    }
    # Case 4: No flood anywhere
    if (sequence[i] == 0 && is_flooded == 0){
      # Don't have to do anything
    }
  }
  max_overlap_ct = max(max_overlap_ct, overlap_ct)
  return(max_overlap_ct + 1)
}


flood_consec <- flood %>% 
  mutate(avg_overlap1 = 1*(AvgAbove & c(AvgAbove[2:length(AvgAbove)],0)),  #current day and next day flooded?
         avg_overlap2 = 1*(AvgAbove & c(0, AvgAbove[1:(length(AvgAbove)-1)])),  #Day before and current day flooded?
         avg_overlap = as.numeric(1*(avg_overlap1 | avg_overlap2)),
         zip3 = substr(Zip_Clean, 1, 3)) 

flood_consec$study_year_idx <- sapply(flood_consec$date, study_year_cat)
#Had to convert NAs to 0's in order to make function run. will reset these back to NAs shortly.
flood_consec$avg_overlap[is.na(flood_consec$avg_overlap)] <- 0

#Directly computes for each 3zip study_year_idx
allstream_summary <-flood_consec %>% 
  group_by(zip3, study_year_idx) %>% 
  summarize(avg_max_consec = max_overlap_length(avg_overlap),
            avg_days_consec = sum(avg_overlap)) %>% 
  ungroup() %>% 
  mutate(zip3 = as.numeric(zip3)) %>% 
  arrange(zip3, study_year_idx)


#Alternatively, can try to aggregrate in the same way i did flood data initially
#Calculate consecutive days per month, per station
consec_summary <-flood_consec %>% 
  group_by(MonthYear, Station_ID, Zip_Clean) %>% 
  summarize(avg_max_consec = max_overlap_length(avg_overlap),
            avg_days_consec = sum(avg_overlap),
            totalAvgAbove = sum(AvgAbove, na.rm = T),
            totalMaxAbove = sum(MaxAbove, na.rm = T),
            totalMinAbove = sum(MinAbove, na.rm = T)) %>% 
  arrange(MonthYear, Station_ID, Zip_Clean)

flood_consec_month_agg2 <- consec_summary %>% 
  group_by(MonthYear, Zip_Clean) %>% 
  summarise(avg_max_consec2 = mean(avg_max_consec, na.rm = T),
            avg_days_consec2 = mean(avg_days_consec, na.rm = T),
            AvgAbove = mean(totalAvgAbove, na.rm = T),
            avgMax = mean(totalMaxAbove, na.rm = T),
            avgMin = mean(totalMinAbove, na.rm = T),
            n = n()) %>% 
  mutate(zip3 = substr(Zip_Clean, 1, 3),
         study_year_idx = sapply(MonthYear, study_year_cat))


sum(is.na(flood_consec_month_agg2)) #only 10 missing data now from 5 observatoins. all from the stations that did not have Zipcodes. 
flood_consec_month_agg2[is.na(flood_consec_month_agg2$Zip_Clean),] #this is what accounts for the missing flooding data I am going to exclude this data and save files
range(flood_consec_month_agg2$n) #Range is 1 to 6. Makes sense as this will be the count of how many statoins per ZCTA

flood_consec_month_agg2 <- na.exclude(flood_consec_month_agg2) #this ultimately gets rids of the stream gauges that have any missing data. For this particular case is getting rid fo the guages that don't have Zip_CLean
sum(is.na(flood_consec_month_agg2))

flood_consec_3ZCTA_Month_Agg <- flood_consec_month_agg2 %>% 
  group_by(zip3, MonthYear, study_year_idx) %>% 
  summarize(avgAvgAbove = mean(AvgAbove, na.rm = T),
            totalAvgAbove = sum(AvgAbove, na.rm = T),
            avgAvgMax = mean(avgMax, na.rm = T),
            totalMaxAbove = sum(avgMax, na.rm = T),
            avgAvgMin = mean(avgMin, na.rm = T),
            totalMinAbove = sum(avgMin, na.rm = T),
            AvgMaxConsec = mean(avg_max_consec2, na.rm = T),
            TotalMaxConsec = sum(avg_max_consec2, na.rm = T),
            AvgDaysConsec = mean(avg_days_consec2, na.rm = T),
            TotalDaysConsec = sum(avg_days_consec2, na.rm = T),
            n = n())

flood_3ZCTA_YrIDX_Agg <- flood_consec_3ZCTA_Month_Agg %>% 
  group_by(zip3, study_year_idx) %>% 
  summarize(totalAvgAbove2 = sum(avgAvgAbove, na.rm = T),
            sd_avgAvgAbove = sd(avgAvgAbove),
            countAvgAbove = sum(totalAvgAbove, na.rm = T),
            totalMaxAbove2 = sum(avgAvgMax, na.rm = T),
            sd_avgMaxAbove = sd(avgAvgMax),
            countMaxAbove = sum(totalMaxAbove, na.rm = T),
            totalMinAbove2 = sum(avgAvgMin, na.rm = T),
            sd_avgMinAbove = sd(avgAvgMin),
            countMinAbove = sum(totalMinAbove, na.rm = T),
            total_AvgMaxConsec = sum(AvgMaxConsec, na.rm = T),
            sd_AvgMaxConsec = sd(AvgMaxConsec, na.rm = T),
            count_MaxConsec = sum(TotalMaxConsec, na.rm = T),
            total_AvgDaysConsec = sum(AvgDaysConsec, na.rm = T),
            sd_AvgDaysConsec = sd(AvgDaysConsec, na.rm = T),
            count_total_DaysConsec = sum(TotalDaysConsec, na.rm = T),
            n = n())

# write_csv(flood_3ZCTA_YrIDX_Agg, path = "flood_3ZCTA_YearIdx_Agg.csv")

#### DATA WRANGLING CENSUS ####
afo <- readxl::read_excel(path = "IAOnly_AFO_ZCTA.xlsx")
producers <- read_csv(file = "2017 Producer.csv")

# Aggregrating the total number of pigs and birds per 5ZCTA
AFO_5ZCTA_Agg <- afo %>%
  mutate(Id2 = ZCTA5CE10,
         zip3 = substr(ZCTA5CE10, 1, 3)) %>% 
  select(uniqueID, MAPLABELNA, Id2, zip3, opStatus, Swine, CattleDair, CattleBeef, Chickens, Turkeys, Horses, SheepLGoat, Other) %>% 
  group_by(Id2) %>% 
  summarize(total_pig = sum(Swine),
            total_bird = sum(Chickens) + sum(Turkeys))

# Aggregating the total number of pigs and birds per 3ZCTA
AFO_3ZCTA_Agg <- afo %>%
  mutate(Id2 = ZCTA5CE10,
         zip3 = substr(ZCTA5CE10, 1, 3)) %>% 
  select(uniqueID, MAPLABELNA, Id2, zip3, opStatus, Swine, CattleDair, CattleBeef, Chickens, Turkeys, Horses, SheepLGoat, Other) %>% 
  group_by(zip3) %>% 
  summarize(total_pig = sum(Swine),
            total_bird = sum(Chickens) + sum(Turkeys))


#Aggregating Producer information
prodicers_5ZCTA_Agg <- producers %>% 
  mutate(Id2 = `Zip Code`,
         zip3 = substr(`Zip Code`, 1, 3)) %>% 
  group_by(Id2) %>% 
  summarize(Total_producers = sum(Value[`Data Item` == "PRODUCERS, (ALL) - NUMBER OF PRODUCERS"]),
            Total_Male_prod = sum(Value[`Data Item` == "PRODUCERS, MALE - NUMBER OF OPERATIONS"]),
            Total_Female_prod = sum(Value[`Data Item` == "PRODUCERS, FEMALE - NUMBER OF OPERATIONS"]),
            Total_PrimaryOcc = sum(Value[`Data Item` == "PRODUCERS, PRIMARY OCCUPATION, FARMING - NUMBER OF PRODUCERS"]),
            Total_LiveOn = sum(Value[`Data Item` == "PRODUCERS, RESIDENCE, ON OPERATION - NUMBER OF PRODUCERS"]))


#this is to get the population level data from the 2010 Census
pop <- readxl::read_xlsx(path =  "2010ZCTACensusData.xlsx")

agegrpandgender <- read_csv(file = "2010_AgeGroups&Sex.csv")
pop <- left_join(pop, agegrpandgender, by = "Id2") %>% 
  left_join(., prodicers_5ZCTA_Agg, by = "Id2") %>% 
  left_join(., AFO_5ZCTA_Agg, by = "Id2")

colnames(pop)
#Getting the Square mile data for Pop density calcs
Zip3_Asqmi <- readxl::read_xlsx(path = "3ZCTA_Areasqmi.xlsx") %>% arrange(zip3)

pop2 <- pop %>% 
  mutate(zip3 = as.numeric(substr(Id2, 1, 3))) %>% 
  select(zip3, Total_Population, Total_households_occupied, `Number - Both sexes; Total population - Under 5 years.x`, #Under 5 at risk pop
         `Number - Both sexes; Total population - 75 to 79 years.x`, `Number - Both sexes; Total population - 80 to 84 years.x`, 
         `Number - Both sexes; Total population - 85 to 89 years.x`, `Number - Both sexes; Total population - 90 years and over.x`, #at risk pops
         
         Total_producers, Total_PrimaryOcc, Total_LiveOn, #All the producer info (so who is a farmer, or works with animals)
         total_pig, total_bird #This is the count for the total number of animals that are flu hosts
  ) %>% 
  mutate(Total_older75 = (`Number - Both sexes; Total population - 75 to 79 years.x` + 
                            `Number - Both sexes; Total population - 80 to 84 years.x` +
                            `Number - Both sexes; Total population - 85 to 89 years.x` + 
                            `Number - Both sexes; Total population - 90 years and over.x`)) %>% 
  left_join(., Zip3_Asqmi, by = "zip3")

pop3 <- pop2 %>% 
  group_by(zip3, Area_sqmi) %>%
  summarise(total_pop = sum(Total_Population),
            Under5 = sum(`Number - Both sexes; Total population - Under 5 years.x`),
            Older75 = sum(Total_older75),
            total_prod = sum(Total_producers, na.rm = T),
            Total_pig = sum(total_pig, na.rm = T),
            Total_bird = sum(total_bird, na.rm = T)) %>% 
  mutate(PerUnder5 = (Under5 / total_pop) *100,
         PerOlder75 = (Older75 / total_pop)*100,
         pop_density = total_pop / Area_sqmi,
         PerProducer = (total_prod / total_pop)*100)
write_csv(pop3, path = "AggANDPpl_Census_3ZCTA_Agg.csv")

#### COMBINING ABOVE DATA SETS ####
#Load in flu data

flu_all <- read_csv(flu_all, path = "FluDiag&SHL_3ZCTA_YrIdx_Aggregated.csv")

#Census Data 
pop <- read_csv(file = "AggANDPpl_Census_3ZCTA_Agg.csv")
#This has 26 unique zip3's where the flu data only has 25
# this is because flu data did not have anything from zip3 = 560. 
# Will continue to only left_join for now to not intorudcue NAs

flu_pop <- left_join(flu_all, pop, by = "zip3")
sum(is.na(flu_pop)) # 0 missing data

#Weather Data
weather_3zip <- read_csv(file = "AllWeather_3zipStudySeasonIdxAgg.csv")
sum(is.na(weather_3zip)) #400. this is mainly from measurements in index year 15. will be ok becuase i don't have flu data for then. 

flood_3zip <- read_csv(file = "flood_3ZCTA_YearIdx_Agg.csv")
sum(is.na(flood_3zip)) #0 missing data

unique(flood_3zip$zip3)
unique(weather_3zip$zip3)

weather_flood_3zip <- left_join(weather_3zip, flood_3zip, by = c("zip3", "study_year_idx"))

sum(is.na(weather_flood_3zip)) # 2544 missing data points. 
#this has a good amount of missing data because flooding i do not have flooding data for each zip3 for each study year index. 


#Mergining flood, weather, population, and flu data together
All_data <- left_join(flu_pop, weather_flood_3zip, by = c("zip3", "study_year_idx"))
sum(is.na(All_data)) #1744. this is likely because the flood data does not start until year_index 3

#Filtering to year_index 3 as this is when we had data for weather, flooding and influenza records
All_data_filt <- All_data %>% filter(study_year_idx > 2)
sum(is.na(All_data_filt)) #944 missing data points

#double checking that all the missing is the flood data
sum(is.na(All_data_filt$totalAvgAbove2)) # 59 missing 
sum(is.na(All_data_filt$sd_avgAvgAbove)) # 59 missing
sum(is.na(All_data_filt$countAvgAbove))  #59 missing
sum(is.na(All_data_filt$totalMaxAbove2)) # 59 missing


All_data_filt <- All_data_filt %>% 
  mutate(vaccrate = total_vacc / total_pop * 10000,
         AsmAtckrate = total_AsmAtck / total_pop * 10000,
         Asmarate = total_Asma / total_pop * 10000)
# going to stop there. there are 10 columns in the flood_3zip dataset. this is where all the data is coming from
# going to save with NAs for now. beucase I can easily exclude on analysis sheet. can also do a sensativity analysis by setting them to zero.
write_csv(All_data_filt, path = "Flood_Flu_Weather_FloodFluSeason_3ZCTA&YrIdx_Missingflood.csv")

#### UNIVARIATE ANALYSES ####
model_data <- read_csv(file = "Flood_Flu_Weather_FloodFluSeason_3ZCTA&YrIdx_Missingflood.csv")
model_data_noNA <- na.exclude(model_data)

#### COR TESTS FOR ENV FACTORS ####
mean(model_data_noNA$flood.avgTAVG)
sd(model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flood.avgTAVG)

mean(model_data_noNA$flu.avgTAVG)
sd(model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flu.avgTAVG)

mean(model_data_noNA$flood.avgTMAX)
sd(model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flood.avgTMAX)

mean(model_data_noNA$flu.avgTMAX)
sd(model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flu.avgTMAX)

mean(model_data_noNA$flood.avgTMIN)
sd(model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flood.avgTMIN)

mean(model_data_noNA$flu.avgTMIN)
sd(model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flu.avgTMIN)

mean(model_data_noNA$flood.avgRhum)
sd(model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flood.avgRhum)

mean(model_data_noNA$flu.avgRhum)
sd(model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flu.avgRhum)

mean(model_data_noNA$flood.avgAhum)
sd(model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flood.avgAhum)

mean(model_data_noNA$flu.avgAhum)
sd(model_data_noNA$flu.avgAhum)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flu.avgAhum)

mean(model_data_noNA$totalAvgAbove2)
sd(model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$totalAvgAbove2)


mean(model_data_noNA$total_AvgDaysConsec)


#### COR TEST FOR OTHER COVARIATES
mean(model_data_noNA$vaccrate)
sd(model_data_noNA$vaccrate)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$vaccrate)

mean(model_data_noNA$Asmarate)
sd(model_data_noNA$Asmarate)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$Asmarate)

mean(model_data_noNA$AsmAtckrate)
sd(model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$AsmAtckrate)

mean(model_data_noNA$pop_density)
sd(model_data_noNA$pop_density)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$pop_density)

mean(model_data_noNA$PerOlder75)
sd(model_data_noNA$PerOlder75)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$PerOlder75)

mean(model_data_noNA$PerUnder5)
sd(model_data_noNA$PerUnder5)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$PerUnder5)

mean(model_data_noNA$PerProducer)
sd(model_data_noNA$PerProducer)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$PerProducer)


#### CORRELATION BETWEEN COVARIATES ####
#flood.TAVG
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$vaccrate)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$Asmarate)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$pop_density)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$flood.avgTAVG, model_data_noNA$PerProducer)

#flu.TAVG
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$vaccrate)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$Asmarate)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$pop_density)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$flu.avgTAVG, model_data_noNA$PerProducer)

#flood.TMAX
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$vaccrate)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$Asmarate)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$pop_density)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$flood.avgTMAX, model_data_noNA$PerProducer)

#flu.TMAX
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$vaccrate)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$Asmarate)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$pop_density)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$flu.avgTMAX, model_data_noNA$PerProducer)

#flood.TMIN
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$vaccrate)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$Asmarate)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$pop_density)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$flood.avgTMIN, model_data_noNA$PerProducer)

#flu.TMIN
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$vaccrate)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$Asmarate)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$pop_density)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$flu.avgTMIN, model_data_noNA$PerProducer)

#flood.Rhum
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$vaccrate)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$Asmarate)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$pop_density)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$flood.avgRhum, model_data_noNA$PerProducer)

#flu.Rhum
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$vaccrate)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$Asmarate)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$pop_density)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$flu.avgRhum, model_data_noNA$PerProducer)

#flood.Ahum
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$vaccrate)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$Asmarate)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$pop_density)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$flood.avgAhum, model_data_noNA$PerProducer)

#flu.Ahum
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$vaccrate)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$Asmarate)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$pop_density)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$flu.avgAhum, model_data_noNA$PerProducer)

#flooding
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$vaccrate)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$Asmarate)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$pop_density)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$totalAvgAbove2, model_data_noNA$PerProducer)

#Vacc Rate
cor.test(model_data_noNA$vaccrate, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$vaccrate, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$vaccrate, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$vaccrate, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$vaccrate, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$vaccrate, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$vaccrate, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$vaccrate, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$vaccrate, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$vaccrate, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$vaccrate, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$vaccrate, model_data_noNA$vaccrate)
cor.test(model_data_noNA$vaccrate, model_data_noNA$Asmarate)
cor.test(model_data_noNA$vaccrate, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$vaccrate, model_data_noNA$pop_density)
cor.test(model_data_noNA$vaccrate, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$vaccrate, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$vaccrate, model_data_noNA$PerProducer)

#Asma Rate
cor.test(model_data_noNA$Asmarate, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$Asmarate, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$Asmarate, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$Asmarate, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$Asmarate, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$Asmarate, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$Asmarate, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$Asmarate, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$Asmarate, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$Asmarate, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$Asmarate, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$Asmarate, model_data_noNA$vaccrate)
cor.test(model_data_noNA$Asmarate, model_data_noNA$Asmarate)
cor.test(model_data_noNA$Asmarate, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$Asmarate, model_data_noNA$pop_density)
cor.test(model_data_noNA$Asmarate, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$Asmarate, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$Asmarate, model_data_noNA$PerProducer)

#AsmaAtck Rate
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$vaccrate)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$Asmarate)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$pop_density)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$AsmAtckrate, model_data_noNA$PerProducer)

#Pop_Density
cor.test(model_data_noNA$pop_density, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$pop_density, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$pop_density, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$pop_density, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$pop_density, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$pop_density, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$pop_density, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$pop_density, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$pop_density, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$pop_density, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$pop_density, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$pop_density, model_data_noNA$vaccrate)
cor.test(model_data_noNA$pop_density, model_data_noNA$Asmarate)
cor.test(model_data_noNA$pop_density, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$pop_density, model_data_noNA$pop_density)
cor.test(model_data_noNA$pop_density, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$pop_density, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$pop_density, model_data_noNA$PerProducer)

#PerOlder75
cor.test(model_data_noNA$PerOlder75, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$PerOlder75, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$vaccrate)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$Asmarate)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$pop_density)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$PerOlder75, model_data_noNA$PerProducer)

#PerUnder5
cor.test(model_data_noNA$PerUnder5, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$PerUnder5, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$vaccrate)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$Asmarate)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$pop_density)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$PerUnder5, model_data_noNA$PerProducer)

#PerProducer
cor.test(model_data_noNA$PerProducer, model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$PerProducer, model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$PerProducer, model_data_noNA$flood.avgTMAX)
cor.test(model_data_noNA$PerProducer, model_data_noNA$flu.avgTMAX)
cor.test(model_data_noNA$PerProducer, model_data_noNA$flood.avgTMIN)
cor.test(model_data_noNA$PerProducer, model_data_noNA$flu.avgTMIN)
cor.test(model_data_noNA$PerProducer, model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$PerProducer, model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$PerProducer, model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$PerProducer, model_data_noNA$flu.avgAhum)

cor.test(model_data_noNA$PerProducer, model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$PerProducer, model_data_noNA$vaccrate)
cor.test(model_data_noNA$PerProducer, model_data_noNA$Asmarate)
cor.test(model_data_noNA$PerProducer, model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$PerProducer, model_data_noNA$pop_density)
cor.test(model_data_noNA$PerProducer, model_data_noNA$PerOlder75)
cor.test(model_data_noNA$PerProducer, model_data_noNA$PerUnder5)
cor.test(model_data_noNA$PerProducer, model_data_noNA$PerProducer)



cor.test(model_data_noNA$total_ILI, model_data_noNA$flu.avgAhum)


#### VIR ANALYSIS ####
standarized <- function(x){
  (x-mean(x)) / sd(x)  
}

stnd_lm_data <- as.data.frame(lapply(model_data_noNA[,c(4, 10, 19, 20, 21, 22, 65, 66, 73, 74, 81, 82, 87, 103, 104, 105)], standarized))



lm.test <- lm(flu_diag_lab ~ PerUnder5 + PerOlder75 + pop_density + PerProducer +  flood.avgTAVG + flu.avgTAVG +
                flood.avgRhum + flu.avgRhum + flood.avgAhum + flu.avgAhum + totalAvgAbove2 + vaccrate +
                AsmAtckrate + Asmarate, data = stnd_lm_data)

lm.test2 <- lm(flu_diag_lab ~ pop_density + PerProducer + flu.avgAhum + totalAvgAbove2 + vaccrate +
                 AsmAtckrate, data = stnd_lm_data)

ili.lm.test <- lm(total_ILI ~ PerUnder5 + PerOlder75 + pop_density + PerProducer +  flood.avgTAVG + flu.avgTAVG +
                    flood.avgRhum + flu.avgRhum + flood.avgAhum + flu.avgAhum + totalAvgAbove2 + vaccrate +
                    AsmAtckrate + Asmarate, data = stnd_lm_data)

ili.lm.test2 <- lm(total_ILI ~ pop_density + PerProducer + flu.avgAhum + totalAvgAbove2 + vaccrate +
                     AsmAtckrate, data = stnd_lm_data)
library(car)
vif(lm.test)
vif(lm.test2)

vif(ili.lm.test)
vif(ili.lm.test2)

#### MULTIVARIATE BAYSIAN CAR MODELS ####
iowa.shp <- st_read("Iowa_Shapefile/Merged510IowaZCTA_proj.shp")
zips <- unique(iowa.shp$Zip_3)
rownames(iowa.shp) <- zips


adj_matrx <- poly2nb(iowa.shp, row.names = zips)
colnames(adj_matrx2) <- zips
rownames(adj_matrx2) <- zips

NumCells=length(adj_matrx)
num=sapply(adj_matrx, length)
adj=unlist(adj_matrx)
sumNumNeigh=length(unlist(adj_matrx))


#### LOAD IN MODEL DATA ####
model_data <- read_csv(file = "Flood_Flu_Weather_FloodFluSeason_3ZCTA&YrIdx_Missingflood.csv")
#This data set has missing data from flood records (590)
sum(is.na(model_data)) #Yup still 944
# sum(is.na(model_data$totalAvgAbove2)) #59 
#going to set these to zero for one aspect of the sensitivity analysis
#Going to exclude this data for main analyses

model_data_noNA <- na.exclude(model_data)

# model_data[is.na(model_data)] <- 0
# sum(is.na(model_data))

#### MODELS #####
flucomplete_model_pois <- function(){
  #Model we are interested in  
  for(i in 1:regions){
    Y[i] ~ dpois(mu[i])
    log(mu[i]) <- (alpha + vaccrate*x1[i] + AsmAtckrate*x2[i] + 
                     flood*x4[i] + 
                     
                     Flu.avgAhum*x6[i] + 
                     
                     pop_density*x7[i] + PerProducer*x10[i] + log(n[i]) + 
                     
                     U[i])
    
  }
  #CAR Distribution + adjency matrix
  U[1:regions] ~ car.normal(adj[], weights[], num[], tau.u)
  for(j in 1:sumNumNeigh){weights[j] <- 1}
  
  ##Priors 
  alpha ~ dnorm(0,0.001)
  vaccrate ~ dnorm(0, 0.001)
  AsmAtckrate ~ dnorm(0, 0.001)
  
  flood ~ dnorm(0, 0.001)
  
  Flu.avgAhum ~ dnorm(0, 0.001)
  
  pop_density ~ dnorm(0, 0.001)
  PerProducer ~ dnorm(0, 0.001)
  
  tau.u ~ dgamma(0.1, 0.1)
  sigma <- sqrt(1/tau.u)
}



#Data to be used in the model
flucomplete_linedata_pois <- list(regions = NumCells,
                                  sumNumNeigh = sumNumNeigh,
                                  adj = adj,
                                  num = num,
                                  Y = model_data_noNA$flu_diag_lab,     #total_flu_diag for flu diagnosis only
                                  x1 = model_data_noNA$vaccrate,
                                  x2 = model_data_noNA$AsmAtckrate,
                                  
                                  x4 = model_data_noNA$totalAvgAbove2,
                                  
                                  x6 = model_data_noNA$flu.avgAhum,
                                  
                                  x7 = model_data_noNA$pop_density,
                                  x10 = model_data_noNA$PerProducer,
                                  n = model_data_noNA$total_pop)




flucomplete_lineinits_pois <- function() {
  list(alpha = 1,  vaccrate = 0, AsmAtckrate = 0, 
       
       flood = 0,
       
       Flu.avgAhum = 0, 
       
       pop_density = 0, PerProducer = 0,
       tau.u = 1, U = c(0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0))
}


set.seed(123)
lineout_flucomplete_nocoda <- bugs(data = flucomplete_linedata_pois, 
                                   inits = flucomplete_lineinits_pois, 
                                   parameters.to.save = c("alpha", "vaccrate", "AsmAtckrate", 
                                                          "flood",
                                                          
                                                          "Flu.avgAhum", 
                                                          
                                                          "pop_density", "PerProducer",
                                                          "sigma"), 
                                   model.file = flucomplete_model_pois, 
                                   n.iter = 2500000,
                                   n.burnin = 200000,
                                   n.chains = 3)
save(lineout_flucomplete_nocoda, file = "Aim 1 Analysis/October NH Models/Flu_CompleteRecords.rda", compress = "bzip2")
# load("Aim 1 Analysis/October NH Models/Flu_CompleteRecords.rda")
exp(lineout_flucomplete_nocoda$summary)
gelman.diag(lineout_flucomplete_nocoda)
#prop > 1
colnames(lineout_flucomplete_nocoda$sims.matrix)
mean(round(exp(lineout_flucomplete_nocoda$sims.matrix[,7]), digits = 2) > 1)

traceplot(as.mcmc.list(lineout_flucomplete_nocoda))


# ILI complete Records
ILIcomplete_model_pois <- function(){
  #Model we are interested in  
  for(i in 1:regions){
    Y[i] ~ dpois(mu[i])
    log(mu[i]) <- alpha + vaccrate*x1[i] + AsmAtckrate*x2[i] + 
      flood*x4[i] + 
      
      Flu.avgAhum*x6[i] + 
      
      pop_density*x7[i] + PerProducer*x10[i] + log(n[i]) + 
      
      U[i]
    
  }
  #CAR Distribution + adjency matrix
  U[1:regions] ~ car.normal(adj[], weights[], num[], tau.u)
  for(j in 1:sumNumNeigh){weights[j] <- 1}
  
  ##Priors 
  alpha ~ dnorm(0,0.001)
  vaccrate ~ dnorm(0, 0.001)
  AsmAtckrate ~ dnorm(0, 0.001)
  
  flood ~ dnorm(0, 0.001)
  
  Flu.avgAhum ~ dnorm(0, 0.001)
  
  pop_density ~ dnorm(0, 0.001)
  PerProducer ~ dnorm(0, 0.001)
  
  tau.u ~ dgamma(0.1, 0.1)
  sigma <- sqrt(1/tau.u)
}



#Data to be used in the model
ILIcomplete_linedata_pois <- list(regions = NumCells,
                                  sumNumNeigh = sumNumNeigh,
                                  adj = adj,
                                  num = num,
                                  Y = model_data_noNA$total_ILI,     #total_flu_diag for flu diagnosis only
                                  x1 = model_data_noNA$vaccrate,
                                  x2 = model_data_noNA$AsmAtckrate,
                                  
                                  x4 = model_data_noNA$totalAvgAbove2,
                                  
                                  x6 = model_data_noNA$flu.avgAhum,
                                  
                                  x7 = model_data_noNA$pop_density,
                                  x10 = model_data_noNA$PerProducer,
                                  n = model_data_noNA$total_pop)




ILIcomplete_lineinits_pois <- function() {
  list(alpha = 1,  vaccrate = 0, AsmAtckrate = 0, 
       
       flood = 0,
       
       Flu.avgAhum = 0,
       
       pop_density = 0, PerProducer = 0,
       tau.u = 1, U = c(0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0))
}

set.seed(123)
lineout_ILIcomplete_nocoda <- bugs(data = ILIcomplete_linedata_pois, 
                                   inits = ILIcomplete_lineinits_pois, 
                                   parameters.to.save = c("alpha", "vaccrate", "AsmAtckrate", 
                                                          "flood",
                                                          
                                                          "Flu.avgAhum", 
                                                          
                                                          "pop_density", "PerProducer", 
                                                          "sigma"), 
                                   model.file = ILIcomplete_model_pois, 
                                   n.iter = 3500000,  #needed to increase for model convergence
                                   n.burnin = 300000, #needed to increase for model convergence
                                   n.chains = 3)
save(lineout_ILIcomplete_nocoda, file = "Aim 1 Analysis/October NH Models/ILI_CompleteRecords.rda", compress = "bzip2")
# load("Aim 1 Analysis/October NH Models/ILI_CompleteRecords.rda")
exp(lineout_ILIcomplete_nocoda$summary)
gelman.diag(lineout_ILIcomplete_nocoda)

colnames(lineout_ILIcomplete_nocoda$sims.matrix)
mean(round(exp(lineout_ILIcomplete_nocoda$sims.matrix[,7]), digits = 2) > 1)

traceplot(as.mcmc.list(lineout_ILIcomplete_nocoda))

#### NO 09 ####
model_data <- read_csv(file = "Aim 1 Analysis/Flood_Flu_Weather_newseason_3ZCTA&YrIdx_Missingflood.csv")
#get rid of missing records
model_data_noNA <- na.exclude(model_data)
#subset out 09 data
model_data_noNA <- model_data_noNA %>% filter(study_year_idx != 5)
flucomplete_model_pois <- function(){
  #Model we are interested in  
  for(i in 1:regions){
    Y[i] ~ dpois(mu[i])
    log(mu[i]) <- alpha + vaccrate*x1[i] + AsmAtckrate*x2[i] +
      flood*x4[i] + 
      
      Flu.avgAhum*x6[i] +
      
      pop_density*x7[i] + PerProducer*x10[i] + log(n[i]) +
      
      U[i]
    
  }
  #CAR Distribution + adjency matrix
  U[1:regions] ~ car.normal(adj[], weights[], num[], tau.u)
  for(j in 1:sumNumNeigh){weights[j] <- 1}
  
  ##Priors 
  alpha ~ dnorm(0,0.001)
  vaccrate ~ dnorm(0, 0.001)
  AsmAtckrate ~ dnorm(0, 0.001)
  
  flood ~ dnorm(0, 0.001)
  
  Flu.avgAhum ~ dnorm(0, 0.001)
  
  pop_density ~ dnorm(0, 0.001)
  PerProducer ~ dnorm(0, 0.001)
  
  tau.u ~ dgamma(0.1, 0.1)
  sigma <- sqrt(1/tau.u)
}



#Data to be used in the model
flucomplete_linedata_pois <- list(regions = NumCells,
                                  sumNumNeigh = sumNumNeigh,
                                  adj = adj,
                                  num = num,
                                  Y = model_data_noNA$flu_diag_lab,     #total_flu_diag for flu diagnosis only
                                  x1 = model_data_noNA$vaccrate,
                                  x2 = model_data_noNA$AsmAtckrate,
                                  
                                  x4 = model_data_noNA$totalAvgAbove2,
                                  
                                  x6 = model_data_noNA$flu.avgAhum,
                                  
                                  x7 = model_data_noNA$pop_density,
                                  x10 = model_data_noNA$PerProducer,
                                  n = model_data_noNA$total_pop)




flucomplete_lineinits_pois <- function() {
  list(alpha = 1,  vaccrate = 0, AsmAtckrate = 0, 
       
       flood = 0,
       
       Flu.avgAhum = 0,
       
       pop_density = 0, PerProducer = 0,
       tau.u = 1, U = c(0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0))
}


set.seed(123)
lineout_flucompleteNo09_nocoda <- bugs(data = flucomplete_linedata_pois, 
                                       inits = flucomplete_lineinits_pois, 
                                       parameters.to.save = c("alpha", "vaccrate", "AsmAtckrate",
                                                              "flood",
                                                              
                                                              "Flu.avgAhum",
                                                              
                                                              "pop_density", "PerProducer", 
                                                              "sigma"), 
                                       model.file = flucomplete_model_pois, 
                                       n.iter = 2500000,
                                       n.burnin = 200000,
                                       n.chains = 3)
save(lineout_flucompleteNo09_nocoda, file = "Aim 1 Analysis/October NH Models/Flu_CompleteRecordsNo09.rda", compress = "bzip2")
# load("Aim 1 Analysis/October NH Models/Flu_CompleteRecordsNo09.rda")
exp(lineout_flucompleteNo09_nocoda$summary)
gelman.diag(lineout_flucompleteNo09_nocoda)

colnames(lineout_flucompleteNo09_nocoda$sims.matrix)
mean(round(exp(lineout_flucompleteNo09_nocoda$sims.matrix[,6]), digits = 2) > 1)


traceplot(as.mcmc.list(lineout_flucompleteNo09_nocoda))

#### CONSEC DAYS FLOODING ####
model_data <- read_csv(file = "Aim 1 Analysis/Flood_Flu_Weather_newseason_3ZCTA&YrIdx_Missingflood.csv")
#get rid of missing records
model_data_noNA <- na.exclude(model_data)
flucomplete_model_pois <- function(){
  #Model we are interested in  
  for(i in 1:regions){
    Y[i] ~ dpois(mu[i])
    log(mu[i]) <- alpha + vaccrate*x1[i] + AsmAtckrate*x2[i] +
      flood*x4[i] + 
      
      Flu.avgAhum*x6[i] +
      
      pop_density*x7[i] + PerProducer*x10[i] + log(n[i]) +
      
      U[i]
    
  }
  #CAR Distribution + adjency matrix
  U[1:regions] ~ car.normal(adj[], weights[], num[], tau.u)
  for(j in 1:sumNumNeigh){weights[j] <- 1}
  
  ##Priors 
  alpha ~ dnorm(0,0.001)
  vaccrate ~ dnorm(0, 0.001)
  AsmAtckrate ~ dnorm(0, 0.001)
  
  flood ~ dnorm(0, 0.001)
  
  Flu.avgAhum ~ dnorm(0, 0.001)
  
  pop_density ~ dnorm(0, 0.001)
  PerProducer ~ dnorm(0, 0.001)
  
  tau.u ~ dgamma(0.1, 0.1)
  sigma <- sqrt(1/tau.u)
}



#Data to be used in the model
flucomplete_linedata_pois <- list(regions = NumCells,
                                  sumNumNeigh = sumNumNeigh,
                                  adj = adj,
                                  num = num,
                                  Y = model_data_noNA$flu_diag_lab,     #total_flu_diag for flu diagnosis only
                                  x1 = model_data_noNA$vaccrate,
                                  x2 = model_data_noNA$AsmAtckrate,
                                  
                                  x4 = model_data_noNA$total_AvgDaysConsec,
                                  
                                  x6 = model_data_noNA$flu.avgAhum,
                                  
                                  x7 = model_data_noNA$pop_density,
                                  x10 = model_data_noNA$PerProducer,
                                  n = model_data_noNA$total_pop)




flucomplete_lineinits_pois <- function() {
  list(alpha = 1,  vaccrate = 0, AsmAtckrate = 0, Asmarate = 0, 
       
       flood = 0,
       
       Flu.avgAhum = 0,
       
       pop_density = 0, PerProducer = 0,
       tau.u = 1, U = c(0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0))
}



set.seed(123)
lineout_flucompleteConsec_nocoda <- bugs(data = flucomplete_linedata_pois, 
                                         inits = flucomplete_lineinits_pois, 
                                         parameters.to.save = c("alpha", "vaccrate", "AsmAtckrate",
                                                                "flood",
                                                                
                                                                "Flu.avgAhum",
                                                                
                                                                "pop_density", "PerProducer",
                                                                "sigma"), 
                                         model.file = flucomplete_model_pois, 
                                         n.iter = 2500000,
                                         n.burnin = 200000,
                                         n.chains = 3)
save(lineout_flucompleteConsec_nocoda, file = "Aim 1 Analysis/October NH Models/Flu_CompleteRecordsConsec.rda", compress = "bzip2")
# load("Aim 1 Analysis/October NH Models/Flu_CompleteRecordsConsec.rda")
exp(lineout_flucompleteConsec_nocoda$summary)
gelman.diag(lineout_flucompleteConsec_nocoda)

traceplot(as.mcmc.list(lineout_flucompleteConsec_nocoda))

colnames(lineout_flucompleteConsec_nocoda$sims.matrix)
mean(round(exp(lineout_flucompleteConsec_nocoda$sims.matrix[,7]), digits = 2) > 1)

#### SENSATIVITY ANALYSES ####
#### SETTING MISSING FLOOD TO 0 ####
#Set missing flood to 0
model_data <- read_csv(file = "Aim 1 Analysis/Flood_Flu_Weather_newseason_3ZCTA&YrIdx_Missingflood.csv")
model_data$totalAvgAbove2[is.na(model_data$totalAvgAbove2)] <- 0

flucomplete_model_pois <- function(){
  #Model we are interested in  
  for(i in 1:regions){
    Y[i] ~ dpois(mu[i])
    log(mu[i]) <- alpha + vaccrate*x1[i] + AsmAtckrate*x2[i] +
      flood*x4[i] + 
      
      Flu.avgAhum*x6[i] +
      
      pop_density*x7[i] + PerProducer*x10[i] + log(n[i]) +
      
      U[i]
    
  }
  #CAR Distribution + adjency matrix
  U[1:regions] ~ car.normal(adj[], weights[], num[], tau.u)
  for(j in 1:sumNumNeigh){weights[j] <- 1}
  
  ##Priors 
  alpha ~ dnorm(0,0.001)
  vaccrate ~ dnorm(0, 0.001)
  AsmAtckrate ~ dnorm(0, 0.001)
  
  flood ~ dnorm(0, 0.001)
  
  Flu.avgAhum ~ dnorm(0, 0.001)
  
  pop_density ~ dnorm(0, 0.001)
  PerProducer ~ dnorm(0, 0.001)
  
  tau.u ~ dgamma(0.1, 0.1)
  sigma <- sqrt(1/tau.u)
}



#Data to be used in the model
flucomplete_linedata_pois <- list(regions = NumCells,
                                  sumNumNeigh = sumNumNeigh,
                                  adj = adj,
                                  num = num,
                                  Y = model_data$flu_diag_lab,     #total_flu_diag for flu diagnosis only
                                  x1 = model_data$vaccrate,
                                  x2 = model_data$AsmAtckrate,
                                  
                                  x4 = model_data$totalAvgAbove2,
                                  
                                  x6 = model_data$flu.avgAhum,
                                  
                                  x7 = model_data$pop_density,
                                  x10 = model_data$PerProducer,
                                  n = model_data$total_pop)




flucomplete_lineinits_pois <- function() {
  list(alpha = 1,  vaccrate = 0, AsmAtckrate = 0, 
       
       flood = 0,
       
       Flu.avgAhum = 0, 
       
       pop_density = 0, PerProducer = 0,
       tau.u = 1, U = c(0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0))
}


set.seed(123)
lineout_flucompleteTo0_nocoda <- bugs(data = flucomplete_linedata_pois, 
                                      inits = flucomplete_lineinits_pois, 
                                      parameters.to.save = c("alpha", "vaccrate", "AsmAtckrate",
                                                             "flood",
                                                             
                                                             "Flu.avgAhum",
                                                             
                                                             "pop_density", "PerProducer", 
                                                             "sigma"), 
                                      model.file = flucomplete_model_pois, 
                                      n.iter = 2500000,
                                      n.burnin = 200000,
                                      n.chains = 3)
save(lineout_flucompleteTo0_nocoda, file = "Aim 1 Analysis/October NH Models/Flu_MissingTo0.rda", compress = "bzip2")
# load("Aim 1 Analysis/October NH Models/Flu_MissingTo0.rda")
exp(lineout_flucompleteTo0_nocoda$summary)
gelman.diag(lineout_flucompleteTo0_nocoda)

colnames(lineout_flucompleteTo0_nocoda$sims.matrix)
mean(round(exp(lineout_flucompleteTo0_nocoda$sims.matrix[,7]), digits = 2) > 1)

traceplot(as.mcmc.list(lineout_flucompleteTo0_nocoda))

#### USING FLOOD IMPUTE ####
#Set missing flood to impute value
model_data <- read_csv(file = "Flood_Flu_Weather_newseason_3ZCTA&YrIdx_Missingflood.csv")
flood_impute <- read_csv("Flood Impute.csv")

model_data2 <- model_data %>% 
  left_join(., flood_impute, by = c("zip3", "study_year_idx"))

a <- model_data %>% group_by(zip3) %>% summarise(flood = sum(totalAvgAbove2),
                                                 flood_noNA = sum(totalAvgAbove2, na.rm = T))

flucomplete_model_pois <- function(){
  #Model we are interested in  
  for(i in 1:regions){
    Y[i] ~ dpois(mu[i])
    log(mu[i]) <- alpha + vaccrate*x1[i] + AsmAtckrate*x2[i] +
      flood*x4[i] + 
      
      Flu.avgAhum*x6[i] +
      
      pop_density*x7[i] + PerProducer*x10[i] + log(n[i]) +
      
      U[i]

    
  }
  #CAR Distribution + adjency matrix
  U[1:regions] ~ car.normal(adj[], weights[], num[], tau.u)
  for(j in 1:sumNumNeigh){weights[j] <- 1}
  
  ##Priors 
  alpha ~ dnorm(0,0.001)
  vaccrate ~ dnorm(0, 0.001)
  AsmAtckrate ~ dnorm(0, 0.001)
  Asmarate ~ dnorm(0, 0.001)
  
  flood ~ dnorm(0, 0.001)
  
  Flu.avgAhum ~ dnorm(0, 0.001)
  
  pop_density ~ dnorm(0, 0.001)
  PerProducer ~ dnorm(0, 0.001)
  
  tau.u ~ dgamma(0.1, 0.1)
  sigma <- sqrt(1/tau.u)
}



#Data to be used in the model
flucomplete_linedata_pois <- list(regions = NumCells,
                                  sumNumNeigh = sumNumNeigh,
                                  adj = adj,
                                  num = num,
                                  Y = model_data2$flu_diag_lab,     
                                  x1 = model_data2$vaccrate,
                                  x2 = model_data2$AsmAtckrate,
                                  
                                  x4 = model_data2$AVG_imputed,
                                  
                                  x6 = model_data2$flu.avgAhum,
                                  
                                  x7 = model_data2$pop_density,

                                  x10 = model_data2$PerProducer,
                                  n = model_data2$total_pop)




flucomplete_lineinits_pois <- function() {
  list(alpha = 1,  vaccrate = 0, AsmAtckrate = 0,
       
       flood = 0,
       
       Flu.avgAhum = 0,
       
       pop_density = 0, PerProducer = 0,
       tau.u = 1, U = c(0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0))
}


set.seed(123)
lineout_flucompleteImpute_nocoda <- bugs(data = flucomplete_linedata_pois, 
                                         inits = flucomplete_lineinits_pois, 
                                         parameters.to.save = c("alpha", "vaccrate", "AsmAtckrate",
                                                                "flood",
                                                                
                                                                "Flu.avgAhum",
                                                                
                                                                "pop_density", "PerProducer", 
                                                                "sigma"), 
                                         model.file = flucomplete_model_pois, 
                                         n.iter = 2500000,
                                         n.burnin = 200000,
                                         n.chains = 3)
save(lineout_flucompleteImpute_nocoda, file = "Aim 1 Analysis/October NH Models/Flu_FloodImpute.rda", compress = "bzip2")
# load("Aim 1 Analysis/October NH Models/Flu_FloodImpute.rda")
exp(lineout_flucompleteImpute_nocoda$summary)
gelman.diag(lineout_flucompleteImpute_nocoda)

colnames(lineout_flucompleteImpute_nocoda$sims.matrix)
mean(round(exp(lineout_flucompleteImpute_nocoda$sims.matrix[,4]), digits = 2) > 1)

traceplot(as.mcmc.list(lineout_flucompleteImpute_nocoda))



#### TABLE 1 Summary Table of Influenza Diagnosis and USGS stream gauge reporting ####
diagcodes_3agg %>% group_by(study_year_idx) %>%
  summarise(total_flu = sum(total_flu_diag),
            total_visit = sum(total_pat))

flood2 <- flood %>% #flood is in code line 359
  mutate(zip3 = substr(Zip_Clean, 1, 3),
         study_year_idx = sapply(MonthYear, study_year_cat))
stations_use <- flood2 %>% 
  filter(study_year_idx > 2 &
           study_year_idx < 13 &
           !is.na(AvgAbove))

flood_sum <- stations_use %>% 
  group_by(study_year_idx) %>% 
  summarise(total_gauges = n_distinct(Station_ID),
            total_report = n_distinct(Station_ID[AvgAbove >= 1]),
            total_flood = sum(AvgAbove, na.rm = T))


#### TABLE 2 Bivariate Correlations between influenza diagnoses and potential covariates ####
#Average Temp Flood season
mean(model_data_noNA$flood.avgTAVG)
sd(model_data_noNA$flood.avgTAVG)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flood.avgTAVG)
#Average Temp Flu season
mean(model_data_noNA$flu.avgTAVG)
sd(model_data_noNA$flu.avgTAVG)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flu.avgTAVG)
#Average Relative Humidity Flood season
mean(model_data_noNA$flood.avgRhum)
sd(model_data_noNA$flood.avgRhum)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flood.avgRhum)
#Average Relative Humidity Flu season
mean(model_data_noNA$flu.avgRhum)
sd(model_data_noNA$flu.avgRhum)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flu.avgRhum)
#Average Absolute Humidity Flood season
mean(model_data_noNA$flood.avgAhum)
sd(model_data_noNA$flood.avgAhum)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flood.avgAhum)
#Average Absolute Humidity Flu season
mean(model_data_noNA$flu.avgAhum)
sd(model_data_noNA$flu.avgAhum)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$flu.avgAhum)
#Average3 Days above flood stage (flooding)
mean(model_data_noNA$totalAvgAbove2)
sd(model_data_noNA$totalAvgAbove2)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$totalAvgAbove2)

#Average Vaccination Rate
mean(model_data_noNA$vaccrate)
sd(model_data_noNA$vaccrate)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$vaccrate)

#Average Asthma Rates
mean(model_data_noNA$Asmarate)
sd(model_data_noNA$Asmarate)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$Asmarate)

#Average Asthma Attack Rates
mean(model_data_noNA$AsmAtckrate)
sd(model_data_noNA$AsmAtckrate)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$AsmAtckrate)

#Population Density
mean(model_data_noNA$pop_density)
sd(model_data_noNA$pop_density)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$pop_density)

#Percent pop Older than 75
mean(model_data_noNA$PerOlder75)
sd(model_data_noNA$PerOlder75)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$PerOlder75)

#Percent pop Younger than 5
mean(model_data_noNA$PerUnder5)
sd(model_data_noNA$PerUnder5)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$PerUnder5)

#Percent pop in Animal Production
mean(model_data_noNA$PerProducer)
sd(model_data_noNA$PerProducer)
cor.test(model_data_noNA$flu_diag_lab, model_data_noNA$PerProducer)





