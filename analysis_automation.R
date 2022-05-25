library(here)
source(
  "https://github.com/tmieno2/OnFarmExperiments/blob/master/Functions/prepare.R?raw=TRUE",
  local = TRUE
)
source(here("./Code/functions.R"))

field_year_ls <- jsonlite::fromJSON(
  file.path(
    here("Data"),
    "field_parameter.json"
  ),
  flatten = TRUE
) %>%
  data.table() %>%
  .[, field_year := paste(farm, field, year, sep = "_")] %>%
  .[ec == "received", ] %>%
  .[, field_year]

ffy <- "Hord_F98_2019"

data_sf <- st_read(here("Data", ffy, "analysis_data.shp")) 

analysis_make_report(ffy, rerun = TRUE)
render(here("Data", ffy, "analysis_report_exp.Rmd"))


