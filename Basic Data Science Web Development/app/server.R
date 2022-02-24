server <- function(input, output, session) {
  
  filt_data <- reactive({
    sel_col_1 <- as.symbol(input$select_col_1)
    sel_col_2 <- as.symbol(input$select_col_2)
    df %>% 
      group_by(continent) %>%
      select(Country_Name, Year, pop, !!sel_col_1, !!sel_col_2) %>%
      # Filter years
      filter(between(Year, input$slider_year[1], input$slider_year[2])) %>%
      # Filter country
      filter(continent %in% input$select_cont) %>%
      # Cond mutate based on if log/linear
      mutate(!!sel_col_1 := case_when(
         isTRUE(input$lin_scale_x) ~ !!sel_col_1,
         isFALSE(input$lin_scale_x) ~ log(!!sel_col_1))
         ) %>%
      mutate(!!sel_col_2 := case_when(
        isTRUE(input$lin_scale_y) ~ !!sel_col_2,
        isFALSE(input$lin_scale_y) ~ log(!!sel_col_2))
      )
  })

  output$plot <- renderPlotly({
    col_x <- as.symbol(input$select_col_1)
    col_y <- as.symbol(input$select_col_2)
    
    filt_df <- filt_data() %>%
      mutate(X = !!col_x, Y = !!col_y)
    
    plot_ly(data = filt_df, type = "scatter", mode = "markers",
            x = as.formula(paste0("~", input$select_col_1)), 
            y = as.formula(paste0("~", input$select_col_2)), 
            color = ~continent, 
            hoverinfo = "text",
            size = ~pop,
            text = ~ paste("Continent:", continent, 
                           "<br>Country:", Country_Name,
                           "<br>Year:", Year,
                           sprintf("<br>%s:", input$select_col_1), X,
                           sprintf("<br>%s:", input$select_col_2), Y)) %>%
      layout(title = sprintf("%s vs. %s", input$select_col_1, input$select_col_2))
  })
  
  # Select all
  observe({
    if (input$sel_all_cntry) updateSelectInput(session, "select_cont", selected = sel_continent)
  })
  
  # Reset
  observe({
    if (input$reset_cntry) updateSelectInput(session, "select_cont", selected = sel_continent[1])
  })
}
