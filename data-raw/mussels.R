## code to prepare NOAA Gulfwatch data on PAHs in Casco Bay blue mussels (Mytilus edulis)

##  Data derived from NOOA Gulfwatch files, available here:
## https://gulfofmaine.org/public/gulfwatch-contaminants-monitoring/findings/
## For additional data on pre-processing, see CBEP archive on shellfish toxics.


library(dplyr)

mussels <- read.csv('data-raw/gulfwatch_pah.csv') %>%
  select(-Original_Code, -Sample_Type)

usethis::use_data(mussels, overwrite = TRUE)
