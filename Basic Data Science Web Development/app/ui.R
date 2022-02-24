optionsUI <- function(id) {
  ns <- NS(id)
  tagList(
    selectInput("select_cont", "Continent:", 
                choices = sel_continent, 
                selected = sel_continent[1], 
                multiple = TRUE),
    fluidRow(
      column(12, 
             actionButton("sel_all_cntry", "Select All"), 
             actionButton("reset_cntry", "Reset"))
    ), 
    hr(),
    h5(strong("Choose Columns:")),
    fluidRow(
      column(3, switchInput(inputId = "lin_scale_x",
                            label = "X-axis",
                            onLabel = "Linear",
                            offLabel = "Log")),
      column(9, selectInput("select_col_1", 
                            label = NULL, 
                            choices = sel_cols, 
                            width = "100%"))
    ),
    fluidRow(
      column(3, switchInput(inputId = "lin_scale_y",
                            label = "Y-axis",
                            onLabel = "Linear",
                            offLabel = "Log")),
      column(9, selectInput("select_col_2", 
                            label = NULL, 
                            choices = sel_cols, 
                            width = "100%"))
    ),
    hr(),
    sliderTextInput(inputId = "slider_year", 
                    label = "Years:",
                    choices = sel_years,
                    selected = c(sel_years[1], sel_years[length(sel_years)]),
                    grid = TRUE)
  )
}

plotUI <- function(id) {
  ns <- NS(id)
  tagList(
    plotlyOutput("plot", height = 700)
  )
}

ui <- fluidPage(
  # linear to log scale
  # change data type (col of df)
  # change year shown 
  
  titlePanel("GapMinder"),
  sidebarLayout(
    sidebarPanel(optionsUI("options")),
    mainPanel(plotUI("plot"))
  )
  
)