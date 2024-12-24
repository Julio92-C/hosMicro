

# Import a CSV File Module

# Module UI function
csvFileUI <- function(id) {
  # `NS(id)` returns a namespace function, which was save as `ns` and will
  # invoke later.
  ns <- NS(id)
  
  
  tagList(
    fluidRow(
      column(3,
             box(
               width = 12,
               collapsible=TRUE,
               status = "black",
               solidHeader = TRUE,
               background = "gray",
               title = "Input Data",
               fileInput(ns("file"), "Choose CSV File",
                         multiple = TRUE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
               tags$hr(),
               checkboxInput(ns('header'), 'Header', TRUE),
               radioButtons(ns('sep'), 'Separator',
                            c(Comma=',',
                              Semicolon=';',
                              Tab='\t'),
                            'Comma'),
               radioButtons(ns('quote'), 'Quote',
                            c(None='',
                              'Double Quote'='"',
                              'Single Quote'="'"),
                            'Double Quote'),
               sliderInput(ns('obs_num'),
                           label = "Number of Observations",
                           value = 12, min = 1, max = 100)
             ),
             
             ),
      
      column(9,
             box(width = 12,
                 solidHeader = TRUE,
                 status = "black",
                 DT::dataTableOutput(ns("data"))
             )  
         )
    )
   
  )
  
}


# Module server function
csvFileServer <- function(id, obs_num) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {
      
      # The selected file, if any
      userFile <- reactive({
        # If no file is selected, don't do anything
        validate(need(input$file, message = FALSE))
        input$file
      })
      
      # The user's data, parsed into a data frame
      dataframe <- reactive({
         read.csv(userFile()$datapath, header=input$header, sep=input$sep, quote=input$quote)
        
      })
      
      # We can run observers in here if we want to
      observe({
        msg <- sprintf("File %s was uploaded", userFile()$name)
        cat(msg, "\n")
        
        print(colnames(dataframe()))
      })
      
      
      
      # Render dataTable
      output$data = DT::renderDataTable({
        DT::datatable(dataframe(),
                      options = list(pageLength = input$obs_num, 
                                     scrollX = TRUE))
      })
      
      # Return the reactive that yields the data frame
      return(dataframe)
      
    }
    
 
  )    
}