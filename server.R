#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(DT)
library(plotly)
library(ggplot2)
library(dplyr)
library(forcats)
library(thematic)
library(gtsummary)
library(paletteer)
library(scales)
library(reshape2)
library(heatmaply)
library(tidyr)
library(VennDiagram)
library(tidyverse)
library(stringr)
library(ggsignif) 
library(vegan)
library(circlize)
library(pheatmap)

# Define server logic required 
# Define server logic required to draw a table
shinyServer(function(input, output) {
  thematic::thematic_shiny()
  
  # Input csvFileServer Module for MDR-hvKp dataSet
    dataframe <- csvFileServer("MDR_hvKp_dataSet", obs_num = obs_num)
    
    # Input csvFileServer Module for MDR-hvKp dataSet
    metadata <- csvFileServer("MetadataSet", obs_num = obs_num)
  # 
  # # Input MDRdataPlotServer Module for MDR-hvKp dataSet
  filter_data <- MDRdataPlotServer("MDR_hvKp_dataPlot", dataframe = dataframe, metadata = metadata)
  # 
  # # input MDRdataStatsServer Module for MDR-hvKp dataSet
  dataStatsServer("hosMicro_dataStat", filter_data = filter_data)
  
  
  
})
