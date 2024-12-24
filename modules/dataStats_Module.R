# Data Stats Module

# Module UI function
dataStatsUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(12,
             box(width = 12, 
                 solidHeader = TRUE, 
                 status = "teal",
                 title = "Basic stats results",
                 collapsible=TRUE,
                 gt_output(ns("my_gt_table"))
             )
             
      )
    )
  )
}




# Module Server function
dataStatsServer <- function(id, filter_data){
  stopifnot(is.reactive(filter_data))
  
  moduleServer(id, function(input, output, session){
    observe({
      # Build a gtsummary table
      output$my_gt_table <- render_gt ({
        data_stats = (filter_data() %>% 
                        tbl_summary(
                          by = SAMPLE,
                          include = c("DATABASE", "GENE", "COVERAGE_PCT","IDENTITY_PCT"),
                          statistic = list(
                            all_continuous() ~ "{mean} ({sd})",
                            all_categorical() ~ "{n} / {N} ({p}%)"
                          ),
                          digits = all_continuous() ~ 2,
                          label = c(COVERAGE_PCT ~ "COVERAGE",
                                    IDENTITY_PCT ~ "IDENTITY"
                          ),
                          sort = list(everything() ~ "frequency"),
                          missing_text = "Missing",
                          missing = "no"
                        )) 
        
        add_p(data_stats, simulate.p.value = TRUE) %>% 
          add_overall() %>%
          as_gt() %>%
          tab_header(md("**Table 1. Summary of Samples groups and genetics elements profile**"))
        
      })
      
      
    })
    
  })
  
}