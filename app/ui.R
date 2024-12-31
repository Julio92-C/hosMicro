#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinythemes)
library(plotly)
library(gt)
library(shinydashboard)
library(shinydashboardPlus)

theme = shinythemes::shinytheme("united")  # <--- Specify theme here

# Global variables
source("modules/input_csvFile_Module.R", local = T)
source("modules/mdr_dataPlot_Module.R", local = T)
source("modules/dataStats_Module.R", local = T)



# Define UI for application that draws a histogram
dashboardPage(
  skin = "blue-light",
  dashboardHeader(title = "hosMicro"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Contact", tabName = "contact", icon = icon("envelope")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    )
  ),
  
  dashboardBody(
    tags$style(".content-wrapper {overflow-y: hidden !important; }"),
    tabItems(
      tabItem(tabName = "home",
              h1("Welcome to hosMicro"),
              h2("This is a web application to explore the hospital microbiome"),
              h3("Please select a tab to begin")
      ),
      
      tabItem(tabName = "dashboard",
              tabBox(
                # title = "Welcome to dashboad section",
                width = 12,
                side = "right", 
                height = "1500px",
                selected = "hosMicro",
                
                # <--- MDR-hvkp start Session--->
                tabPanel("hosMicro", 
                         tabsetPanel(
                           tabPanel("Description",
                                    fluidRow(
                                      column(6,
                                             br(),
                                             tags$img(src = "AMR_Pathogens.png", alt = "AMR-Pathogens", width = "100%", height= "auto")
                                             
                                      ),
                                      
                                      column(6,
                                             br(),
                                             box(
                                               width = 12,
                                               status = "teal",
                                               solidHeader = TRUE,
                                               title = "Genomic diversity and population structure analysis of Hospital microbiome",
                                               tags$p("The spread of antimicrobial-resistant bacteria globally is a pressing issue that has captured my attention."),
                                               tags$ul(
                                                 tags$li("It's concerning to see how certain bacteria are becoming resistant to the drugs we use to treat them, leading to the rise of superbugs that are difficult to control."),
                                                 tags$li("In 2017, WHO priority pathogens: A list of antibiotic-resistant bacteria assessed to be of highest priority for new antibiotic development."),
                                                 tags$li("In 2019, nearly 1.3 million deaths including 140 thousand newborns were caused by AMR. This is expected to rise to 10 million deaths by 2050. "),
                                                 tags$li("By understanding the patterns and drivers of antimicrobial resistance, we can develop strategies to prevent its spread and improve patient outcomes.")
                                               ),
                                               tags$p("Therefore, this web-based application has been developed to facilitate the exploration of
                                              genomic diversity and population structure analysis of the hospital microbiome."),
                                               
                                               tags$p("The application is designed
                                              to provide a user-friendly interface for researchers and clinicians to analyze
                                              and visualize the metagenomic data generated from the hospital environment."),
                                               
                                               tags$p("It offers a range of tools
                                              and features to help users identify taxonomy, track the spread of the pathogen,
                                              and understand their AMR, VFs and MGEs profiles.")
                                             )
                                      )
                                    )
                                    
                                    
                           ),
                           tabPanel("Data",
                                    br(),
                                    # Import csvFileUI Module for MDR-hvKp dataset
                                    csvFileUI("MDR_hvKp_dataSet")
                                    
                           ),
                           tabPanel("MetaData",
                                    br(),
                                    # Import csvFileUI Module for MDR-hvKp dataset
                                    csvFileUI("MetadataSet")
                                    
                           ),
                           
                           tabPanel("Plot",
                                    br(),
                                    # Import MDR dataPlot Module UI
                                    MDRdataPlotUI("MDR_hvKp_dataPlot")
                                    
                           ),
                           
                           tabPanel("Stats Analysis",
                                    br(),
                                    # Import MDR dataPlot Module UI
                                    dataStatsUI("hosMicro_dataStat")
                                    
                           )
                         )
                         
                         
                )
                # <--- MDR-hvkp End--->
                
                
                
              ),
              # <--- Project Section End --->
              
              
              # <---- Dashboard end section ---->
      ),
      
      
      # <--- Start Contact Section ----->
      tabItem(tabName = "contact",
              p("Name: Julio C. Ortega Cambara"),
              p("Email: 32104617@student.uwl.ac.uk"),
              p("PhD C. Bioinformatics"),
              p("School of Biomedical Sciences"),
              p("University of West London"),
              p("St Mary's Rd, London W5 5RF")
              
              
              ),
      
      # <--- Start About Section ----->
            tabItem(tabName = "about",
              h1("About hosMicro"),
              p("This is a web application to explore the hospital microbiome, a critical tool for understanding the complex microbial ecosystems within healthcare settings. By providing detailed insights into the types and behaviors of microorganisms present in hospitals, this application aids in identifying potential sources of infections and tracking the spread of antibiotic-resistant bacteria. Such knowledge is invaluable for developing targeted infection control strategies, improving patient outcomes, and ensuring a safer hospital environment. Moreover, the application facilitates ongoing research and collaboration among healthcare professionals, fostering a proactive approach to managing hospital-acquired infections and enhancing overall public health."),
              
      )
    )
  )
)
