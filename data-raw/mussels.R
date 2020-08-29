## code to prepare NOAA Gulfwatch data on PAHs in Casco Bay blue mussels (Mytilus edulis)

##  Data derived from NOOA Gulfwatch files, available here:
## https://gulfofmaine.org/public/gulfwatch-contaminants-monitoring/findings/
## For additional data on pre-processing, see CBEP archive on shellfish toxics.


library(dplyr)
library(readr)

mussels <- read_csv('data-raw/gulfwatch_pah.csv',
                    col_types = cols(
                      Original_Code = col_character(),
                      Site = col_character(),
                      Year = col_integer(),
                      Sample_Type = col_character(),
                      PAH = col_character(),
                      Flag = col_logical(),
                      Concentration = col_double())) %>%
  select(-Original_Code, -Sample_Type)

usethis::use_data(mussels, overwrite = TRUE)
