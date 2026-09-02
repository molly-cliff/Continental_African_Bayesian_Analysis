library(dplyr)
library(INLA)
library(sf)
library(spdep)
library(sp)
library(caret)
library(pROC)
library(tidyr)
library(broom)
library(corrplot)
library(car)
library(ggplot2)
library(reshape2)

data_join <-function() {
  library(dplyr)
  dataset_files <- list(
    rainfall = "final_rainfall_data.csv",
    meridonal_zonal = "final_meridonal_zonal_data.csv",
    aod = "final_aod_data.csv",
    humidity = "final_humidity_data.csv",
    windspeed = "final_windspeed_data.csv",
    temp = "final_temp_data.csv"
  )
  
  # Load outbreak data
  spatiotemporaloutbreaks <- read.csv("full-outbreak-matched.csv", stringsAsFactors = FALSE)
  
  # Check columns in outbreak data
  print("Column names in outbreak data:")
  print(names(spatiotemporaloutbreaks))

  
  
  
  
  
  # Start with FALSE for all rows
bf1$vaccine <- 0







#4 campaigns in Nigeria- first December 2011
bf1$vaccine[
  bf1$ADMN1_country %in% c("Bauchi Nigeria",
"Gombe Nigeria",
"Jigawa Nigeria",
"Katsina Nigeria",
"Zamfara Nigeria")& 
  bf1$year_month2 >= 201112
] <- 1


#4 campaigns in Nigeria- second December 2012


bf1$vaccine[
  bf1$ADMN1_country %in% c("Kano Nigeria",
"Sokoto Nigeria",
"Borno Nigeria",
"Yobe Nigeria")& 
  bf1$year_month2 >= 201212
] <- 1


#4 campaigns in Nigeria- third November 2013
bf1$vaccine[
  bf1$ADMN1_country %in% c("Adamawa Nigeria",
"Federal Capital Territory Nigeria",
"Kadun, Nigeria",
"Kebbi Nigeria",
"Nasarawa Nigeria",
"Niger Nigeria",
"Plateau Nigeria",
"Taraba Nigeria")& 
  bf1$year_month2 >= 201311
] <- 1


#4 campaigns in Nigeria- forth October 2014

bf1$vaccine[
  bf1$ADMN1_country %in% c("Anambra Nigeria",
"Benue Nigeria",
"Ebonyi Nigeria",
"Enugu Nigeria",
"Cross River Nigeria",
"Kogi Nigeria",
"Kwara Nigeria",
"Imo Nigeria",
"Oyo Nigeria")& 
  bf1$year_month2 >= 201410

] <- 1


#Cameroon May 2013- Northern districts vaccination
bf1$vaccine[
  bf1$ADMN1_country %in% c("Nord Cameroon", "Est Extrême-Nord Cameroon")& 
    bf1$year_month2 >= 201305
] <- 1


#2 campaigns in South Sudan- first March 2015


bf1$vaccine[
  bf1$ADMN1_country %in% c("Central Equatoria South Sudan",
"Eastern Equatoria South Sudan" ,
"Lakes South Sudan",
"North Bahr-al-Ghazal South Sudan",
"Warap South Sudan",
"West Bahr-al-Ghazal South Sudan")& 
  bf1$year_month2 >= 201503
] <- 1


#2 campaigns in South Sudan- second May 2018




bf1$vaccine[
  bf1$ADMN1_country %in% c("Jungoli South Sudan",
                           "Unity  South Sudan" ,
                           "Lakes South Sudan",
                           "Upper Nile South Sudan",
                           "West Equatoria South Sudan")& 
    bf1$year_month2 >= 201805
] <- 1





#Limited campaign carried out in high risk regions
#Phase 1, scheduled for November 2011, will cover the health regions of N'Djamena, Mayo Kebbi Est, and Chari Baguirmi. It targets a population of 1,809,657 people and will require 2,113,918 doses of MenAfriVac.


bf1$vaccine[
  bf1$district_country %in% c("Loug Chari Chad", "Baguirmi Chad", "N'Djamena Chad", "Mayo-Boneye Chad", 
                              "Mont Illi Chad", "Kabbia Chad") & 
    bf1$year_month2 >= 201111
] <- 1

bf1$vaccine[
  bf1$ADMN1_country %in% c(
    "Mayo-Kebbi Est Chad",
    "Chari-Baguirmi Chad") & 
    bf1$year_month2 >= 201111
] <- 1



#Phase 2, scheduled for December 2011, will cover the health regions of Guéra, Logone Oriental, Mandoul, Tandjilé, and Moyen Chari. It targets a population of 2,455,601 people and will require 2,897,609 doses of MenAfriVac.

bf1$vaccine[
  bf1$ADMN1_country %in% c(
    "Guéra Chad",
    "Logone Oriental Chad",
    "Mandoul Chad",
    "Tandjilé Chad",
    "Moyen-Chari Chad"
    ) & 
    bf1$year_month2 >= 201112
] <- 1



#Phase 3, scheduled for November 2012, will cover the health regions of Hadjer Lamis, Lac, Logone Occidental, and Mayo Kebbi Ouest. It targets a population of 1,761,188 people and will require 2,078,201 doses of MenAfriVac.


bf1$vaccine[
  bf1$ADMN1_country %in% c(
    
    "Hadjer-Lamis Chad",
    "Lac Chad",
    "Logone Occidental Chad",
    "Mayo-Kebbi Ouest Chad"
  ) & 
    bf1$year_month2 >= 20211
] <- 1




#Phase 4, scheduled for December 2012, will cover the health regions of Batha, Borkou, Tibesti, Ennedi, Kanem, Bahr El-Ghazal, Ouaddaï, Sila, Salamat, and Wadi Fira. It targets a population of 2,614,112 people and will require 3,084,652 doses of MenAfriVac (see Figure 1 below).

bf1$vaccine[
  bf1$ADMN1_country %in% c(
    "Batha Chad",
    "Borkou Chad",
    "Tibesti Chad",
    "Ennedi Est Chad",
    "Ennedi Ouest Chad",
    "Kanem Chad",
    "Barh el Ghazel Chad",
    "Ouaddaï Chad",
    "Salamat Chad",
    "Sila Chad",
    "Wadi Fira Chad"
    
  ) & 
    bf1$year_month2 >= 20211
] <- 1





#Northern high risk areas of togo vaccinated in November 2014

bf1$vaccine[
  bf1$ADMN1_country %in% c( "Centre Togo", "Kara Togo","Plateaux Togo","Savanes Togo") & 
    bf1$year_month2 >= 201411
] <- 1



#Two capaigns in guinea, one in June 2014 and other in August 2015

bf1$vaccine[
  bf1$district_country %in% c( "Siguiri Guinea", 
                               "Mandiana Guinea") & 
    bf1$year_month2 >= 201406
] <- 1
    
    
bf1$vaccine[
  bf1$district_country %in% c( "Gaoual Guinea",
"Koundara Guinea",
"Dabola Guinea",
"Dinguiraye Guinea",
"Faranah Guinea",
"Kissidougou Guinea",
"Kankan, Guinea",
"Kérouané Guinea",
"Kouroussa Guinea",
"Koubia Guinea",
"Labé Guinea",
"Lélouma Guinea",
"Mali Guinea",
"Tougué Guinea",
"Beyla Guinea") & 
  bf1$year_month2 >= 201508
] <- 1    
    
#Three campaigns, country wide in Ethiopia moving West to East in Oct 2013, Oct 2014 and Oct 2015


bf1$vaccine[
  bf1$district_country == "Debub Gondar Ethiopia",
"Debub Wollo Ethiopia",
"Mirab Gojjam Ethiopia",
"Misraq Gojjam Ethiopia",
"Bahir Dar Special Zone Ethiopia",
"Semen Gondar Ethiopia",
"Semen Wello Ethiopia",
"Agew Awi Ethiopia",
"Agnuak Majang Ethiopia",
"Nuer Ethiopia",
"Keffa Ethiopia",
"Misraq Wellega Ethiopia",
"Bench Maj Ethiopia",
"Debub Omo Ethiopia",
"Sheka Ethiopia",
"Mirab Welega Ethiopia",
"Kelem Wellega Ethiopia",
"Jimma Ethiopia",
"Asosa Ethiopia",
"Kemashi Ethiopia",
"Metekel Ethiopia",
"Debubawi Ethiopia",
"Mehakelegnaw Ethiopia",
"Mi'irabawi Ethiopia",
"Misraqawi Ethiopia",
"Semien Mi'irabaw Ethiopia" & bf1$year_month2 >= 201310
] <- 1


bf1$vaccine[
  bf1$ADMN1_country %in% c( "Addis Abeba Ethiopia",
                            "Oromia Ethiopia",
                            "Somali Southern Nations, Nationalities Ethiopia") & 
    bf1$year_month2 >= 201411
] <- 1


bf1$vaccine[
  bf1$COUNTRY %in% c("Ethiopia") & 
    bf1$year_month2 >= 201510
] <- 1

#Niger, campaign countrywide December 2010-Jan 2011
bf1$vaccine[
  bf1$district_country == "Niger Filingué" & bf1$year_month2 >= 201009
] <- 1

bf1$vaccine[
  bf1$district_country %in% c( "Téra Niger", "Tillabéry Niger", "Ouallam Niger", "Say Niger",
                               "Niamey Niger", "Kollo Niger", "Dosso Niger", "Boboye Niger") & 
    bf1$year_month2 >= 201012
] <- 1

bf1$vaccine[
  bf1$COUNTRY %in% c("Niger") & 
    bf1$year_month2 >= 201101
] <- 1


#Mali pilot in September 2010, followed by camapigns December 2010 and November 2011
#details of campaign one not sure where
bf1$vaccine[
  bf1$district_country %in% c( "Dioïla  Mali", "Fanna Mali") & 
    bf1$year_month2 >= 201009
] <- 1


bf1$vaccine[
  bf1$ADMN1_country %in% c( "Koulikoro  Mali") & 
    bf1$year_month2 >= 201012
] <- 1

bf1$vaccine[
  bf1$ADMN1_country %in% c( "Kayes Mali",
                            "Sikasso Mali",
                            "Mopti Mali",
                            "Timbuktu Mali",
                            "Gao Mali",
                            "Kidal Mali") & 
    bf1$year_month2 >= 201111
] <- 1



# Mauritania, high risk counties in OCtober 2014, but no detail of which ones

bf1$vaccine[bf1$COUNTRY %in% c( "Mauritania"
)& bf1$year_month2 >= 201410
] <- 1

# Benin high risk counties in Nov 2012

bf1$vaccine[
  bf1$district_country %in% c( "Banikoara Benin", "Gogounou Benin", "Kandi Benin", "Karimama Benin","Malanville Benin",
      "Ségbana Benin",
                               "Boukoumbé Benin", "Cobly Benin", "Kérou  Benin", "Kouandé Benin", "Matéri Benin",
                               "Natitingou Benin", "Péhunco Benin", "Tanguiéta Benin", "Toucountouna Benin", 
                               "Bembèrèkè Benin", "Kalalé Benin", "Nikki Benin", "N'Dali Benin", "Parakou Benin",
                               "Pèrèrè Benin", "Sinendé Benin", "Tchaourou Benin",  "Bantè Benin",
                               "Dassa-Zoumè Benin", "Glazoué Benin", "Ouèssè Benin", "Savalou Benin", "Savè Benin",
                               "Bassila Benin", "Copargo Benin", "Djougou Urbain Benin", "Djougou Rural Benin",
                               "Ouaké Benin") & 
    bf1$year_month2 >= 201211
] <- 1





# CID, high risk counties in December 2014

bf1$vaccine[
  bf1$district_country %in% c("Gbeke, Côte d'Ivoire",
"Marahoué, Côte d'Ivoire",
"Bounkani, Côte d'Ivoire",
"Bagoué, Côte d'Ivoire",
"Hambol, Côte d'Ivoire",
"Tchologo, Côte d'Ivoire",
"Poro, Côte d'Ivoire",
"Béré, Côte d'Ivoire",
"Gontougo, Côte d'Ivoire",
"Bafing, Côte d'Ivoire",
"Kabadougou, Côte d'Ivoire",
"Worodougou, Côte d'Ivoire",
"Folon, Côte d'Ivoire") & 
    bf1$year_month2 >= 201412
] <- 1




# CAR, high risk counties in March 2017 and other regions in May 2017
bf1$vaccine[
  bf1$ADMN1_country %in% c(
    "Mambéré-Kadéï Central African Republic",
    "Nana-Mambéré Central African Republic",
    "Sangha-Mbaéré Central African Republic",
    "Ouham Central African Republic",
    "Ouham-Pendé Central African Republic",
    "Kémo Central African Republic",
    "Ouaka Central African Republic",
    "Nana-Grébizi Central African Republic",
    "Haute-Kotto Central African Republic",
    "Vakaga Central African Republic",
    "Bamingui-Bangoran Central African Republic" ) &  bf1$year_month2 >= 201703] <- 1

bf1$vaccine[
  bf1$district_country %in% c(
"M'Baïki Central African Republic",
"Boda Central African Republic",
"Bimbo Central African Republic",
"Mobay, Central African Republic",
"Alindao Central African Republic",
"Kembé Central African Republic",
"Bangassou Central African Republic" ) &  bf1$year_month2 >= 201705] <- 1

bf1$vaccine[
  bf1$ADMN1_country %in% c(
"Bangui Central African Republic",
"Haut-Mbomou Central African Republic" ) &  bf1$year_month2 >= 201703] <- 1



# Kenya, 5 high risk counties in June 2019

bf1$vaccine[
  bf1$ADMN1_country %in% c(
    "Mandera, Kenya", "Marsabit, Kenya","Turkana, Kenya","Wajir, Kenya","West Pokot, Kenya"
  ) &  bf1$year_month2 >= 201906] <- 1




# DRC, high risk areas March 2016
bf1$vaccine[
  bf1$district_country %in% c(
"Buta Democratic Republic of the Congo",
"Buta (ville) Democratic Republic of the Congo",
"Isiro Democratic Republic of the Congo",
"Aru Democratic Republic of the Congo",
"Bunia Democratic Republic of the Congo",
"Kisangani Democratic Republic of the Congo",
"Butembo Democratic Republic of the Congo",
"Goma Democratic Republic of the Congo",
"Bukavu Democratic Republic of the Congo",
"Uvira Democratic Republic of the Congo") & 
  bf1$year_month2 >= 201603
] <- 1







# Ghana, high risk areas in October 2012, rest of country in July 2016
bf1$vaccine[
  bf1$district_country %in% c(  "Bole Ghana","Bunkpurugu Nakpanduri Ghana","Central Gonja Ghana","East Gonja Ghana",
                                "Chereponi Ghana","East Mamprusi Ghana","Gushegu Ghana","Karaga Ghana",
                                "Kpandai Ghana" ,    "Nanumba North Ghana",   "Nanumba South Ghana", "Saboba Ghana", 
                                "Savelugu Ghana","Sawla-Tuna-Kalba Ghana",  "Tamale Ghana","Sagnerigu Ghana","Tolon Ghana" ,
                                "Kumbungu Ghana","West Gonja Ghana" ,"West Mamprusi Ghana",  "Yendi Ghana" ,
                                "Zabzugu Ghana","Bawku Ghana" ,"Bawku West Ghana","Bolgatanga Ghana", 
                                "Bongo Ghana","Builsa South Ghana" ,"Builsa North Ghana","Garu Ghana",
                                "Kasena Nankana East Ghana","Kasena Nankana West Ghana","Talensi Ghana" ,"Bolgatanga Ghana", 
                                "Bongo Ghana" ,"Builsa South Ghana" ,"Builsa North Ghana", "Garu Ghana" ,
                                "Kasena Nankana East Ghana" ,"Kasena Nankana West Ghana","Talensi Ghana" ) & 
    bf1$year_month2 >= 201210
] <- 1

bf1$vaccine[bf1$COUNTRY %in% c("Ghana"
)& bf1$year_month2 >= 201608
] <- 1



# Gambia, whole country in November 2013

bf1$vaccine[bf1$COUNTRY %in% c("Gambia"
)& bf1$year_month2 >= 201311
] <- 1



# Sudan, roll out in two phases, first sub selection of ADMN1 districts in October 2012 adn rest in September 2013

bf1$vaccine[
  bf1$ADMN1_country %in% c("Khartoum Sudan", "Al Jazirah Sudan", "Blue Nile Sudan", "Al Qadarif Sudan", 
                           "West Darfur Sudan", "North Darfur Sudan", "South Darfur Sudan", "Central Darfur Sudan", "East Darfur Sudan"
  ) &  bf1$year_month2 >= 201210] <- 1

bf1$vaccine[
  bf1$country %in% c("Sudan") &  bf1$year_month2 >= 201309] <- 1

# Burkina Faso, 1 test distrct in September 2010 and then whole country from december 2010
bf1$vaccine[
  bf1$district_country %in% c("Sanmatenga Burkina Faso") & 
    bf1$year_month2 >= 201009
] <- 1

bf1$vaccine[
  bf1$COUNTRY == "Burkina Faso" & bf1$year_month2 >= 201012
] <- 1


#Burundi 1 country wide vacination programme December 2018
bf1$vaccine[
  bf1$COUNTRY == "Burundi" & bf1$year_month2 >= 201812
] <- 1



#Senegal country wide in November 2012
bf1$vaccine[bf1$COUNTRY %in% c("Senegal"
)& bf1$year_month2 >= 201211
] <- 1


#Guinea Bissau country wide in June 2016

bf1$vaccine[bf1$COUNTRY %in% c("Guinea Bissau" 
)& bf1$year_month2 >= 201606
] <- 1

#Eritrea, 2 full country vaccination programmes in Nov and December 2017
#add admn1 country in there as well as new variable

bf1$vaccine[
  bf1$ADMN1_country %in% c("Debubawi Keyih Bahri Eritrea",
                              "Semenawi Keyih Bahri Eritrea","Maekel Eritrea"
                              ) &  bf1$year_month2 >= 201711] <- 1

bf1$vaccine[
  bf1$ADMN1_country %in% c("Anseba Eritrea",
                              "Debub Eritrea","Gash Barka Eritrea"
  ) &  bf1$year_month2 >= 201712] <- 1


#Uganda vaccination, 39 high risk districts in Jan 2017- no details of where they are so assuming full country
bf1$vaccine[bf1$COUNTRY %in% c("Uganda")& bf1$year_month2 >= 201701
] <- 1


print(names(bf1))

# Start with outbreak summary as the base
merged_data <- bf1

# Merge each environmental dataset
for (name in names(dataset_files)) {
  env_data <- read.csv(dataset_files[[name]], stringsAsFactors = FALSE)
  
  # If the env dataset has year/month/district_country, rebuild the code correctly
  if (all(c("year", "month", "district_country") %in% names(env_data))) {
    env_data$year_month <- paste(env_data$year, sprintf("%02d", env_data$month), sep = " ")
    env_data$code <- paste(env_data$year_month, env_data$district_country, sep = " ")
  }
  
  # Deduplicate by code
  if ("code" %in% names(env_data)) {
    env_data <- env_data %>% distinct(code, .keep_all = TRUE)
  } else {
    warning(paste("Column 'code' not found in", name, "dataset. Deduplication skipped."))
  }
  
  # Left join onto merged_data
  merged_data <- merged_data %>%
    left_join(env_data, by = "code")
}

# Final column names
print("Final column names after merging:")
print(names(merged_data))
print(table(merged_data$outbreak_2, merged_data$country))

# Select & rename useful columns
myvars <- c("year", "month", "name_2", "outbreak_2", "district_country", "country", "code",
            "rainfall", "eastward_wind", "north_wind", "aod", "humidity", "windspeed", "temp")

# Keep only columns that exist (protects against missing columns)
existing_vars <- intersect(myvars, names(merged_data))
merged_data <- merged_data[, existing_vars]

# Data structure and summary
print("---- Structure of merged data ----")
str(merged_data)

print("---- Summary statistics ----")
print(summary(merged_data))

print("---- Missing values per column ----")
na_counts <- sapply(merged_data, function(x) sum(is.na(x)))
print(na_counts)

print("---- Proportion of missing values per column ----")
na_proportion <- round(na_counts / nrow(merged_data), 3)
print(na_proportion)

# Summary for numeric columns
numeric_cols <- sapply(merged_data, is.numeric)
print("---- Means and ranges for numeric columns ----")
print(summary(merged_data[, numeric_cols]))

# Unique counts for character columns
categorical_cols <- sapply(merged_data, is.character)
print("---- Unique values for character columns ----")
unique_counts <- sapply(merged_data[, categorical_cols, drop = FALSE], function(x) length(unique(x)))
print(unique_counts)

# Example rows
print("---- Example rows ----")
print(head(merged_data, 10))


landcover_data <- read.csv("final_landcover_data.csv", stringsAsFactors = FALSE)

# Create merge key in landcover data
landcover_data$year_district <- paste(landcover_data$year, landcover_data$district_country, sep = " ")

# Create same merge key in merged_data
merged_data$year_district <- paste(merged_data$year, merged_data$district_country, sep = " ")

# Deduplicate landcover data (optional, in case there are duplicates)
landcover_data <- landcover_data %>% distinct(year_district, .keep_all = TRUE)

# Merge landcover onto merged_data
merged_data <- merged_data %>%
  left_join(landcover_data, by = "year_district")

# Optional: Remove the temporary merge key
merged_data$year_district <- NULL
merged_data<- na.omit(merged_data)
# Check final structure
print("---- Structure after merging landcover ----")
str(merged_data)

# Data structure and summary
print("---- Structure of merged data ----")
str(merged_data)

print("---- Summary statistics ----")
print(summary(merged_data))

print("---- Missing values per column ----")
na_counts <- sapply(merged_data, function(x) sum(is.na(x)))
print(na_counts)

print("---- Proportion of missing values per column ----")
na_proportion <- round(na_counts / nrow(merged_data), 3)
print(na_proportion)

# Summary for numeric columns
numeric_cols <- sapply(merged_data, is.numeric)
print("---- Means and ranges for numeric columns ----")
print(summary(merged_data[, numeric_cols]))

# Unique counts for character columns
categorical_cols <- sapply(merged_data, is.character)
print("---- Unique values for character columns ----")
unique_counts <- sapply(merged_data[, categorical_cols, drop = FALSE], function(x) length(unique(x)))
print(unique_counts)

# Example rows
print("---- Example rows ----")
print(head(merged_data, 10))

write.csv(merged_data, "joined_outbreak_vaccine_data.csv", row.names = FALSE)
}