library(shiny)
library(stringr)
library(plotly)
library(tidyr)
library(shinyWidgets)

current_dir <- getwd()
data_path <- file.path(current_dir, "../data/gapminder_clean.csv")

df <- read.csv(data_path) %>% 
  select(-1) %>% 
  drop_na() %>%
  # Remove ending . first, then replace any . or more with _.
  rename_with(.fn = function(col) gsub("\\.+", "_", gsub("\\.$", "", col)))

sel_cols <- df %>% 
  select(where(is.numeric), -Year, -pop) %>% 
  names() 

sel_continent <- unique(df$continent)

sel_years <- sort(unique(df$Year))
