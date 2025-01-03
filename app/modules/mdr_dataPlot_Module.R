# MDR-hvKp dataPlot Module


# Module UI function
MDRdataPlotUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(3,
             box(width = 12, 
                 solidHeader = TRUE, 
                 collapsible=TRUE,
                 status = "black",
                 title = "Choose a Database:",
                 background = "gray",
                 
                 selectInput(ns("DATABASE"),
                             label = "Database",
                             choices = c("card",
                                         "vfdb",
                                         "plasmidfinder"
                             ),
                             selected = "card",
                             multiple = TRUE
                 ),
                 
                 h4("Filter data-set by:"),
                 sliderInput(ns("NAME"),
                             label = "Taxa count",
                             value = 20, min = 0, max = 100
                             
                 ),
                 sliderInput(ns("GENE"),
                             label = "Gene count",
                             value = 10, min = 0, max = 100
                             
                 ),
                 
                 sliderInput(ns("COVERAGE_PCT"),
                             label = "Gene Coverage (%)",
                             value= 90, min = 80, max = 100

                             ),
                 sliderInput(ns("IDENTITY_PCT"),
                             label = "Gene Identity (%)",
                             value = 90, min = 80, max = 100

                 ),
                 
                 
                 
                 fluidRow(
                   column(6,
                          actionButton(ns("plot_button"), icon = icon("chart-bar"), "Plot", 
                                       style = "
                                       width: 100%;
                                       background-color: #428bca;
                                       color: white;
                                       border: none;
                                       border-radius: 5px;
                                       padding: 10px 20px;
                                       font-size: 16px;
                                       font-weight: bold;"
                          )
                   ),
                   column(6,
                          downloadButton(ns("download"), "Report",
                                         style = "
                                         width: 100%;
                                         background-color: #428bca;
                                         color: white;
                                         border: none;
                                         border-radius: 5px;
                                         padding: 10px 20px;
                                         font-size: 16px;
                                         font-weight: bold;"
                          )  
                   )
                 )
                 
                 
                 
                 
             )
      ),
      
      column(9,
             fluidRow(
               
               infoBoxOutput(ns("general_stats"), width = 8),
               
               valueBoxOutput(ns("Drugs_count"), width = 4)
               
               
             ),
             
             fluidRow(
               valueBoxOutput(ns("ARG_count"), width = 4),
               
               valueBoxOutput(ns("VFs_count"), width = 4),
               
               valueBoxOutput(ns("MGEs_count"), width = 4)
               
               
               
               
             ),
             
             fluidRow(
               
               tabBox(width = 12,
                      height = 700,
                      
                      tabPanel("Taxa Distribution",
                               fluidRow(
                                 
                                   plotlyOutput(ns("freq_plot"), height = 500)
                                 
                               ),
                               
                               fluidRow(
                                 column(4,
                                        
                                 ),
                                 column(4,
                                       selectInput(ns("mainVar_freq"),
                                                   label = "Choose a genotype:",
                                                   choices = c("NAME", "GENE"),
                                                   selected = "NAME",
                                                   multiple = FALSE,
                                       )
                                 ),
                                 
                                 column(4,
                                   
                                 )
                               ),
                               
                               fluidRow(
                                 column(7,
                                        plotOutput(ns("taxaCount_violinPlot"), height = 500),
                                 ),
                                 column(5,
                                        plotOutput(ns("Air_vs_Surface"), height = 400, width = 400),
                                        
                                        
                                 ),
                                 
                               ),
                               
                               
                      ),
                      
                      
                      tabPanel("Relative Abundance",
                               fluidRow(
                                 plotlyOutput(ns("Relative_abundance"), height = 500)
                               )
                               
                      ),
                      
                      
                      
                      tabPanel("ARG/VFs/MGEs profile",
                               fluidRow(
                                 plotOutput(ns("ARG_profile"), height = 500),
                               ),
                               
                               fluidRow(
                                 plotOutput(ns("GeneticsElements_sample"), height = 500)
                               )
                               
                               

                      ),
                      
                      tabPanel("ARG/VFs/MGEs co-distribution",
                               fluidRow(
                                 plotlyOutput(ns("PCoA_Analysis"), height = 500),
                               ),
                               
                               fluidRow(
                                 column(6,
                                        plotOutput(ns("AirSurf_violinPlot"), height = 500)  
                                 ),
                                 column(6,
                                        plotOutput(ns("geneticElements_occurrence"), height = 500, width = 500)


                                 )
                                 
                               )
                                 
                               
                      ),
                       
                      tabPanel("Air microbiome",
                               fluidRow(
                                 plotOutput(ns("Air_microbiome"), height = 900),
                               ),
                               
                               fluidRow(
                                 column(3,
                                        
                                 ),
                                 column(6,
                                        box(width = 12, 
                                            solidHeader = TRUE, 
                                            collapsible=TRUE,
                                            status = "black",
                                            title = "Dynamic Row Control for Heatmap:",
                                            background = "gray",
                                            # Slider input for selecting a range of values
                                            sliderInput(ns("row_count"), 
                                                "Select a range:",
                                                 min = 1, 
                                                 max = 226, 
                                                 value = c(1, 50))
                                    
                                        )
                                        
                                        
                                 ),
                                 column(3,
                                        
                                 )
                               ),
                               
                               
                      ),
                      
                      tabPanel("Key Pathogens/AMR",
                               fluidRow(plotOutput(ns("Key_Pathoges_AMR"), height = 600)),
                               
                          )
                     )
                      
               )
               
             )
      )
      
    )
  
  
}

# Module server function
MDRdataPlotServer <- function(id, dataframe, metadata){
  stopifnot(is.reactive(dataframe))
  stopifnot(is.reactive(metadata))
  
  moduleServer(id, function(input, output, session){
    
    # Filter data by user input
    filter_data <- eventReactive(input$plot_button, {
      validate(
        need(dataframe(), "Please input a data-sets as a csv file"),
        need(input$DATABASE, "Please select a database")
      )
      
      # Group the data by name and count the occurrences
      name_counts <- dataframe() %>% 
        group_by(NAME) %>%
        summarize(count = n())
      
      gene_counts <- dataframe() %>% 
        group_by(GENE) %>%
        summarize(count = n())
      
      # print(gene_counts)


      filter_data <- dataframe() %>%
        filter(NAME %in% name_counts$NAME[name_counts$count > input$NAME],
               GENE %in% gene_counts$GENE[gene_counts$count > input$GENE],
               COVERAGE_PCT >= input$COVERAGE_PCT,
               IDENTITY_PCT >= input$IDENTITY_PCT,
               DATABASE %in% input$DATABASE)
      
    })
    
    # Download plots to a HTML file
    output$download <- downloadHandler(
      filename = function() {
        paste("report-", Sys.Date(), ".html", sep="")
      },


      content = function(file) {
        # Code to generate the report goes here
        
        # Set up parameters to pass to Rmd document
        params <- list(filter_data = filter_data(),
                       dataframe = dataframe(),
                       metadata = metadata(),
                       mainVar_freq = input$mainVar_freq,
                       row_count = input$row_count)
        

        rmarkdown::render("reports/shinyReport.Rmd", output_file = file, 
                          params = params, 
                          envir = new.env(parent = globalenv())
                          )
      }
    )
    
    # Display an infoBox with general stats
    output$general_stats <- renderInfoBox({
      infoBox(
        title = "General stats",
        value = paste("No. of Reads:", length(unique(filter_data()$SEQUENCE_ID)), "|",
                      "Samples:", length(unique(filter_data()$SAMPLE)), "|",
                      "Taxa:", length(unique(filter_data()$NAME)), "|",
                      "Genes:", length(unique(filter_data()$GENE)), "|",
                      "Resistance:", length(unique(filter_data()$RESISTANCE)) 
                      
                      ),
        subtitle = "Taxonomy classification carried out using Kraken2. Profiling of ARGs, VFs, and MGEs with Abricate",
        color = "aqua",
        icon = icon("chart-line"),
        fill = TRUE
      )
    })
    
    output$ARG_count <- renderValueBox({
      arg_data <- filter_data() %>%
        filter(DATABASE == "card") %>%
        distinct(NAME, GENE) %>%
        group_by(GENE) %>%
        arrange(GENE)
      
      total_ARG_count <- length(unique(arg_data$GENE))
     
      
      valueBox(
        value = paste0(total_ARG_count, "/",
                       round(100 * total_ARG_count/length(unique(filter_data()$NAME)), 2), "%"),
        subtitle = "Total Antimacrobial Resistance genes", icon("dna"),
        color = "purple"
      )
    })
    
    
    
    output$VFs_count <- renderValueBox({
      vfs_data <- filter_data() %>%
        filter(DATABASE == "vfdb") %>%
        distinct(NAME, GENE) %>%
        group_by(GENE) %>%
        arrange(GENE)
      
      total_VFs_count <- length(unique(vfs_data$GENE))
      
      
      valueBox(
        value = paste0(total_VFs_count, "/",
                       round(100 * total_VFs_count/length(unique(filter_data()$NAME)), 2), "%"),
        subtitle = "Total Virulence Factors", icon("dna"),
        color = "purple"
      )
    })
    
    output$MGEs_count <- renderValueBox({
      mges_data <- filter_data() %>%
        filter(DATABASE == "plasmidfinder") %>%
        distinct(NAME, GENE) %>%
        group_by(GENE) %>%
        arrange(GENE)
      
      total_mges_count <- length(unique(mges_data$GENE))
      
      
      valueBox(
        value = paste0(total_mges_count, "/",
                       round(100 * total_mges_count/length(unique(filter_data()$NAME)), 2), "%"),
        subtitle = "Total Mobile Genetic Elements", icon("dna"),
        color = "purple"
      )
    })
    
    output$Drugs_count <- renderValueBox({
      drugs_data <- filter_data() %>%
        filter(DATABASE == "card") %>%
        separate_longer_delim(RESISTANCE, delim = ";") %>%
        distinct(GENE, RESISTANCE) %>%
        group_by(GENE) %>%
        arrange(GENE) 
         
      # drugs_class <- drugs_data$RESISTANCE
      # print(drugs_class)
      
      total_drugs_count <- length(unique(drugs_data$RESISTANCE))
      
      valueBox(
        value = paste0(total_drugs_count, "/",
                       round(100 * total_drugs_count/length(unique(filter_data()$NAME)), 2), "%"),
        subtitle = "Total Drug classes", icon("capsules"),
        color = "orange"
      )
    })
    
    # Color palettes
    palt1 <- paletteer_d("palettesForR::Cranes", n = 200)
    palt1 <- as.character(palt1)
    
    
    
    # Render ST type frequency plot based on the filter data sets
    output$freq_plot <- renderPlotly({
      # Plot ST type frequency
      a <- ggplot(filter_data(), aes(x= fct_infreq(.data[[input$mainVar_freq]]),
                                     ))
      a <- a + geom_bar(fill = "skyblue")
      a <- a + theme_classic()
      a <- a + geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5)  # Add count labels on top of the bars
      a <- a + labs(fill = "Legend")
      a <- a + theme(legend.position = "Bottom")
      a <- a + theme(legend.justification = "center")
      a <- a + theme(axis.text = element_text(size = 10.5))
      a <- a + theme(plot.title = element_text(hjust = 0.5))
      a <- a + theme(axis.text.x = element_text(angle = 45))
      a <- a + labs(x="", y="Count")
      a
      
      # Plot Model type frequency with PLOTLY ####
      ST_frequency_plot <- plotly::ggplotly(a)
      ST_frequency_plot
      
    })
    
    
    # Render Air vs Surface venn Diagram based on the taxa count
    output$Air_vs_Surface <- renderPlot({
      # Pivot the Sample values into columns
      df_wide <- dataframe() %>%
        arrange(SAMPLE) %>%
        pivot_wider(names_from = SAMPLE, values_from = c(TAXID)) %>%
        group_by(NAME) %>%
        summarise(across(16:75, ~ paste(., collapse = ", "))) %>%
        mutate(across(2:61, ~ str_replace_all(., "(NA,|,NA|,NA,|NA|,| )", "")))
      
      colnames(df_wide)
      
      # Assign the values of col1 as row names
      # Air dataset
      air_dataset <- df_wide %>%
        unite("Air", 2:28, sep = ",") %>%
        select(NAME, Air) %>%
        mutate(across(2, ~ str_replace_all(., "(,| )", "")))
      
      
      # Surface dataset
      Surface_dataset <- df_wide %>%
        unite("Surface", 29:61, sep = ",") %>%
        select(NAME, Surface) %>%
        mutate(across(2, ~ str_replace_all(., "(,| )", "")))
      
      # Merge Air and Surface dataset
      air_Sur_dataset <- merge(air_dataset, Surface_dataset, by="NAME", all = TRUE)
      
      # Assign the values of col1 as row names
      air_Sur_dataset <- air_Sur_dataset %>%
        remove_rownames %>% 
        column_to_rownames(var="NAME") %>%
        mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))
      
      
      # Convert from a dataframe to a matrix
      sample2arg_data_matrix <- data.matrix(air_Sur_dataset, rownames.force = NA)
      
      # Count NA in data sets
      sum(is.na(sample2arg_data_matrix))
      
      # Replace all NA values with 0 values
      sample2arg_data_matrix <- sample2arg_data_matrix %>% replace(is.na(.), 0)
      
      # Convert to a binary matrix
      sample2arg_data_Bmatrix <- as.matrix((sample2arg_data_matrix > 0) + 0)
      
      # Create a list of sets
      sets <- apply(sample2arg_data_Bmatrix, 2, function(col) which(col == 1))
      names(sets) <- colnames(sample2arg_data_Bmatrix)
      # print(sets)
      # class(sets)
      
      ### Plot the Venn diagram
      venn.plot <- venn.diagram(
        x = sets,
        category.names = c("Air", "Surface"),
        fill = c(Air="#04b5d9", Surface="#aeb6b8"), 
        alpha = 0.3,
        height = 300, 
        width = 300,
        filename = NULL,
        cat.cex = 1.5,  # Adjust this value to change the font size
        cex = 2       # Adjust this value to change the font size of the numbers
      )
      
      grid.newpage()
      grid.draw(venn.plot)
      
      
      
    })
    
    
    
     # Render Virulence-Reistance heatmap plot based on the filter data sets
     output$Relative_abundance <- renderPlotly({
       tryCatch({
       
       # Filter by Species and counts greater than 30
       df_species <- filter_data() %>%
         group_by(NAME) %>%
         arrange(NAME) %>%
         group_by(SAMPLE, NAME) %>%
         summarise(Count = n(), .groups = "drop") %>%
         group_by(SAMPLE) %>%
         mutate(Percentage = Count / sum(Count) * 100)
       
         # Annotation col colors for taxa
         # Description variable
         taxa <- as.factor(filter_data()$NAME)
         print(taxa)
         
         # Generate a color palette based on the number of levels in gene_family
         t <- length(levels(taxa))
         print(t)
         t_palette <- paletteer_d("palettesForR::Cranes", n = t)
         t_palette
         
         
         
         #  stacked barplot where each bar represents a sample, 
         #  and the segments within each bar represent the relative abundance of different taxa
          ra <- ggplot(df_species, aes(x = factor(SAMPLE), y = Percentage, fill = NAME)) +
           geom_bar(stat = "identity") +
           scale_fill_manual(values = t_palette) +
           labs(x = "Air and Surface sample groups", y = "Percentage", fill = "NAME") +
           theme_classic() +
           theme(legend.position="top") +
           theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
           theme(text = element_text(size = 16))
         
         # Plot Coverage vs seq_length with PLOTLY ####
         Relative_abundace_plotly <- plotly::ggplotly(ra)
         Relative_abundace_plotly
         
         
         
       }, error = function(e) {
         # Handle the error
         # You can display an error message or take any other appropriate action
         showNotification("An error occurred while generating the plot. Number of requested colors greater than this palette can offer which is 256.",
                          type = "error")

         # Return an empty plot
         return(NULL)
       })
       
    })
     
     # Render ARG, VFs and MGEs distribution / Air and Surface samples groups
     output$taxaCount_violinPlot <- renderPlot({
       tryCatch({
         # Filter the data set by AMR, VFs and MGEs
         abri_kraken2_filtered <- dataframe() %>%
           filter(grepl("card|vfdb|plasmidfinder", DATABASE)) %>%
           arrange(NAME) %>%
           group_by(SAMPLE, GENE) %>%
           summarise(Gene_count = n(), .groups = 'drop') %>%
           ungroup() %>% 
           mutate(Gene_log_count = log(Gene_count + 1))  # Adding 1 to avoid log(0)
         
         # Rename samples
         abri_kraken2_filtered[1:536, 1] <- "Air"
         abri_kraken2_filtered[537:1669, 1] <- "Surface"
         
         
         # Calculate p-value for significance
         p_value <- t.test(Gene_log_count ~ SAMPLE, data = abri_kraken2_filtered)$p.value
         
         t_test_result <- t.test(Gene_log_count ~ SAMPLE, data = abri_kraken2_filtered)
         # print(t_test_result)
         
         
         # Generate the violin plot with statistical significance for AMR
         TaxaCount <- ggplot(abri_kraken2_filtered, aes(x = SAMPLE, y = Gene_log_count, fill = SAMPLE)) +
           geom_violin(trim = FALSE, scale = "width", width = 0.6) +
           #geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
           scale_fill_manual(values = c("#04b5d9", "#aeb6b8")) +
           labs(title = "",
                x = "Sample groups",
                y = "Log(Taxa Count)") +
           theme_classic() +
           theme(legend.position = "none") +
           theme(text = element_text(size = 20)) +
           #geom_signif(comparisons = list(c("Air", "Surface")), 
           #            map_signif_level = TRUE, 
           #            textsize = 16) +
           # annotate("text", x = 1.5, y = max(abri_kraken2_filtered$Gene_log_count) + 1.5, 
           #         label = paste("p =", format(p_value, digits = 2)), 
           #         size = 4, color = "black")+
           theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 14)) +
           theme(axis.text.y = element_text(size = 16))
         
         
         # Plot Air vs Surface violin plot
         TaxaCount
         
         
         
       }, error = function(e) {
         # Handle the error
         # You can display an error message or take any other appropriate action
         showNotification("An error occurred while generating the plot. Adjust the level of each variable on the left to the minimum value.",
                          type = "error")
         
         # Return an empty plot
         return(NULL)
       })
       
       
     })
     
     # ARG Profile
     output$ARG_profile <- renderPlot({
       Gene_count <- table(filter_data()$DATABASE)
       
       # PLot the Gene count by database 
       pie_chart <- pie(Gene_count, labels = paste0(names(Gene_count), " / ", Gene_count, " / ", round(100 * Gene_count/sum(Gene_count), 2), "%"),
                        # main="Gene count by database",
                        col = c("#b02d02", "#5c0404", "#c47a02"),
                        cex = 1.3, # Adjust label size
                        radius = 1, # adjust the size
                        xlim = c(-1.5, 1.5) # Adjust x-axis limits to create more space for labels
       )
       
       # Plot pie chart
       pie_chart
                        
     })
     
     # ARG, VFs and MGEs by samples
     output$GeneticsElements_sample <- renderPlot({
       # Count occurrences
       data_summary <- filter_data() %>%
         group_by(SAMPLE, DATABASE) %>%
         summarise(Count = n())
       
       # Calculate percentages
       data_summary <- data_summary %>%
         group_by(SAMPLE) %>%
         mutate(Percentage = Count / sum(Count) * 100)
       
       
       
       #  stacked barplot where each bar represents a sample, 
       #  and the segments within each bar represent the counts of different genes
       stacked_barplot <- ggplot(data_summary, aes(x = factor(SAMPLE), y = Percentage, fill = DATABASE)) +
         geom_bar(stat = "identity") +
         scale_fill_manual(values = c("#b02d02", "#5c0404", "#c47a02")) +
         labs(x = "Sample", y = "Percentage", fill = "DATABASE") +
         theme_classic() +
         theme(legend.position="top") +
         theme(text = element_text(size = 20)) +
         theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) +
         theme(axis.text.y = element_text(size = 14))
       
       # Plot the graph
       stacked_barplot
       
     })
    
     # Genetic elements co-occurrence on Air and Surface samples groups
     output$geneticElements_occurrence <- renderPlot ({
       
       # Pivot the Sample values into columns
       df_wide_ann_rows <- dataframe() %>%
         arrange(SAMPLE) %>%
         pivot_wider(names_from = DATABASE, values_from = c(GENE)) %>%
         group_by(NAME) %>%
         summarise(across(15:17, ~ paste(., collapse = ", "))) %>%
         mutate(across(2:4, ~ str_replace_all(., "(NA,|,NA|,NA,|NA| )", "")))
       
       
       
       # Assign the values of Sample as row names
       df_wide_ann_rows <- df_wide_ann_rows %>%
         na.omit()  %>%
         remove_rownames %>%
         column_to_rownames(var="NAME") %>%
         mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))
       
       
       # Convert from a dataframe to a matrix
       sample2arg_data_matrix <- data.matrix(df_wide_ann_rows, rownames.force = NA)
       
       # Count NA in data sets
       sum(is.na(sample2arg_data_matrix))
       
       # Replace all NA values with 0 values
       sample2arg_data_matrix <- sample2arg_data_matrix %>% replace(is.na(.), 0)
       
       # Convert to a binary matrix
       sample2arg_data_Bmatrix <- as.matrix((sample2arg_data_matrix > 0) + 0)
       
       # Create a list of sets
       sets <- apply(sample2arg_data_Bmatrix, 2, function(col) which(col == 1))
       names(sets) <- colnames(sample2arg_data_Bmatrix)
       # print(sets)
       # class(sets)
       
      ### Plot the Venn diagram
       venn.plot <- venn.diagram(
         x = sets,
         category.names = c("VFs", "ARGs", "MGEs"),
         fill = c(card="#b02d02", plasmidfinder="#5c0404", vfdb="#c47a02"), 
         alpha = 0.5,
         height = 50, width = 50,
         filename = NULL,
         cat.cex = 1.5,  # Adjust this value to change the font size
         cex = 2       # Adjust this value to change the font size of the numbers
       )
       grid.newpage()
       grid.draw(venn.plot)
       
       
     })
     
      
     # Render PCoA based on the ARG, VFs and MGEs profile
     output$PCoA_Analysis <- renderPlotly({
       suppressWarnings({
       # Annotation rows
       # Pivot the Genes values into columns
       df_wide_ann_rows <- dataframe() %>%
         arrange(GENE) %>%
         pivot_wider(names_from = GENE, values_from = c(NAME)) %>%
         group_by(SAMPLE) %>%
         summarise(across(16:292, ~ paste(., collapse = ", "))) %>%
         mutate(across(2:277, ~ str_replace_all(., "(NA,|,NA|,NA,|NA| )", "")))
       
       colnames(df_wide_ann_rows)
       
       # Assign the values of Sample as row names and empty as rows as NA value
       df_wide_ann_rows <- df_wide_ann_rows %>%
         na.omit()  %>%
         remove_rownames %>%
         column_to_rownames(var="SAMPLE") %>%
         mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))
       
       
       # Convert from a dataframe to a matrix
       sample2arg_data_matrix <- data.matrix(df_wide_ann_rows, rownames.force = NA)
       
       # Count NA in data sets
       sum(is.na(sample2arg_data_matrix))
       
       # Replace all NA values with 0 values
       sample2arg_data_matrix <- sample2arg_data_matrix %>% replace(is.na(.), 0)
       
       # Calculate Bray-Curtis distances
       bray_curtis_dist <- vegdist(sample2arg_data_matrix, method = "bray")
       
       # Perform PCoA
       pcoa_result <- cmdscale(bray_curtis_dist, eig = TRUE, k = 2)
       
       # Calculate percentage of variance explained
       eig_values <- pcoa_result$eig
       var_explained <- eig_values / sum(eig_values) * 100
       
       # Create a data frame for plotting
       pcoa_df <- data.frame(SAMPLE = rownames(sample2arg_data_matrix),
                             PC1 = pcoa_result$points[, 1],
                             PC2 = pcoa_result$points[, 2])
       
       
       # Merge the metadata and sample dataset
       ann_pcoa <- inner_join(pcoa_df, metadata(), by="SAMPLE")
       
       # Annotation col colors for Location
       # Description variable
       location_col <- as.factor(ann_pcoa$LOCATION)
       
       # Generate a color palette based on the number of levels in gene_family
       l <- length(levels(location_col))
       l_palette <- paletteer_d("ggsci::default_igv", n = l)
       l_palette
       
       
       
       # Plot the results
       pcaoa_plot <- ggplot(ann_pcoa, aes(x = PC1, y = PC2, label = SAMPLE, colour = LOCATION)) +
         geom_point(size = 3) +
         #geom_text(vjust = 1.5, hjust = 1.5) +
         #stat_ellipse(geom = "polygon", aes(fill = location), alpha = 0.2) +
         stat_ellipse(lwd = 0.8) + # Change line width
         scale_color_manual(values = l_palette) +
         #scale_fill_manual(values = l_palette) +
         labs(title = "",
              color = "LOCATION",
              x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
              y = paste0("PC2 (", round(var_explained[2], 2), "%)")
         ) +
         theme_classic() +
         theme(text = element_text(size = 16))
       
       # ggplot graph
       # pcaoa_plot
       
       # Plot ST vs No. Resistance Genes with PLOTLY ####
       PCoA_plotly <- plotly::ggplotly(pcaoa_plot)
       PCoA_plotly
       
       })

     })
     
     # Render ARG, VFs and MGEs distribution / Air and Surface samples groups
     output$AirSurf_violinPlot <- renderPlot({
       tryCatch({
         # Filter the data set by AMR, VFs and MGEs
         abri_kraken2_filtered <- dataframe() %>%
           filter(grepl("card", DATABASE)) %>%
           arrange(NAME) %>%
           group_by(SAMPLE, GENE) %>%
           summarise(Gene_count = n(), .groups = 'drop') %>%
           ungroup() %>% 
           mutate(Gene_log_count = log(Gene_count + 1))  # Adding 1 to avoid log(0)
         
         # Rename samples
         abri_kraken2_filtered[1:93, 1] <- "Air"
         abri_kraken2_filtered[94:342, 1] <- "Surface"
         
         
         # Calculate p-value for significance
         p_value <- t.test(Gene_log_count ~ SAMPLE, data = abri_kraken2_filtered)$p.value
         
         t_test_result <- t.test(Gene_log_count ~ SAMPLE, data = abri_kraken2_filtered)
         print(t_test_result)
         
         
         # Generate the violin plot with statistical significance for AMR
         AirSurf <- ggplot(abri_kraken2_filtered, aes(x = SAMPLE, y = Gene_log_count, fill = SAMPLE)) +
           geom_violin(trim = FALSE, scale = "width", width = 0.6) +
           #geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
           scale_fill_manual(values = c("#04b5d9", "#aeb6b8")) +
           labs(title = "",
                x = "Sample groups",
                y = "Log(Number of ARGs)") +
           theme_classic() +
           theme(legend.position = "none") +
           theme(text = element_text(size = 20)) +
           geom_signif(comparisons = list(c("Air", "Surface")), 
                       map_signif_level = TRUE, 
                       textsize = 16) +
           annotate("text", x = 1.5, y = max(abri_kraken2_filtered$Gene_log_count) + 1.5, 
                    label = paste("p =", format(p_value, digits = 2)), 
                    size = 4, color = "black")+
           theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 16)) +
           theme(axis.text.y = element_text(size = 16))
         
         
         # Plot Air vs Surface violin plot
         AirSurf
         
         
         
       }, error = function(e) {
         # Handle the error
         # You can display an error message or take any other appropriate action
         showNotification("An error occurred while generating the plot. Adjust the level of each variable on the left to the minimum value.",
                          type = "error")
         
         # Return an empty plot
         return(NULL)
       })
       
       
     })
     
       #  # Key Pathogens-AMR co-occurrence
      output$Key_Pathoges_AMR <- renderPlot({
          # Filter and mutate the dataset
          abri_kraken2_filtered <- dataframe() %>%
            filter(grepl("card", DATABASE)) %>%
            mutate(GENE = case_when(
              grepl("vanR_gene_in_vanA_cluster", GENE) ~ "vanRA",
              grepl("vanH_gene_in_vanA_cluster", GENE) ~ "vanHA",
              grepl("vanS_gene_in_vanA_cluster", GENE) ~ "vanSA",
              TRUE ~ GENE
            ),
            NAME = ifelse(grepl("Rahnella aquatilis CIP 78.65 = ATCC 33071", NAME), "Rahnella aquatilis", NAME)) %>%
            arrange(NAME) %>%
            mutate(RESISTANCE = ifelse(nchar(RESISTANCE) > 20, "Multi-drug", RESISTANCE))
          
          # Define modify_strings function
          modify_strings <- function(df, column_name, target_strings, replacement_strings) {
            if (length(target_strings) != length(replacement_strings)) {
              stop("target_strings and replacement_strings must be of the same length")
            }
            for (i in seq_along(target_strings)) {
              df <- df %>%
                mutate(!!sym(column_name) := str_replace_all(!!sym(column_name), target_strings[i], replacement_strings[i]))
            }
            return(df)
          }
          
          # List of target and replacement strings
          target_strings <- c("cephalosporin;penam", "carbapenem;penam", "Cephalosporin;Penem")
          replacement_strings <- rep("Beta-lactam", length(target_strings))
          
          # Call Modify String function
          abri_kraken2_filtered <- modify_strings(abri_kraken2_filtered, "RESISTANCE", target_strings, replacement_strings)
          
          # Clean the RESISTANCE strings
          abri_kraken2_filtered <- abri_kraken2_filtered %>%
            mutate(RESISTANCE = str_replace_all(RESISTANCE, "_", " ") %>% str_to_title())
          
          # Filter by species
          df <- abri_kraken2_filtered %>%
            filter(grepl("Pseudomonas aeruginosa|Enterococcus faecium|Staphylococcus aureus|Acinetobacter baumannii|Klebsiella pneumoniae|Klebsiella quasipneumoniae|Acinetobacter lwoffii|Enterobacter cloacae complex|Enterobacter hormaechei|Escherichia coli|Staphylococcus aureus|Rahnella aquatilis", NAME))
          
          # Create a connection matrix
          categories <- unique(c(sort(df$SAMPLE), sort(df$NAME), sort(df$GENE), sort(df$RESISTANCE)))
          mat <- matrix(0, nrow = length(categories), ncol = length(categories), dimnames = list(categories, categories))
          
          for (i in 1:nrow(df)) {
            mat[df$SAMPLE[i], df$NAME[i]] <- 1
            mat[df$NAME[i], df$GENE[i]] <- 1
            mat[df$GENE[i], df$RESISTANCE[i]] <- 1
          }
          
          # Set colors for categories
          color_sectors <- c(
            paletteer_d("ggsci::default_igv", n = length(unique(df$SAMPLE))),
            paletteer_d("ggsci::default_igv", n = length(unique(df$NAME))),
            paletteer_d("ggsci::default_igv", n = length(unique(df$GENE))),
            paletteer_d("ggsci::default_igv", n = length(unique(df$RESISTANCE)))
          )
          
          ## reset the graphic parameters and internal variables
          circos.clear()
          
          print("Printing Chord Diagram...")
          
          # Graph parameters
          circos.par(track.height = 0.1, start.degree = 120, gap.degree = 2, canvas.xlim = c(-1, 1), canvas.ylim = c(-1, 1), circle.margin = c(1, 1), unit.circle.segments = 500)
          
          # Create the chord diagram
          chordDiagram(mat, transparency = 0.5, annotationTrack = "grid", scale = FALSE, directional = 1, diffHeight = mm_h(3), grid.col = color_sectors, preAllocateTracks = list(track.height = 0.1, unit.circle.segments = 100, start.degree = 90, scale = TRUE))
          
          # Add colors to the sectors
          circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
            sector.index = get.cell.meta.data("sector.index")
            circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.8), cex = 1)
          }, bg.border = NA)

      })
      
      
      # Air microbiome heatmap
      output$Air_microbiome <- renderPlot({
        # Pivot the Sample values into columns
        df_wide <- dataframe() %>%
          arrange(SAMPLE) %>%
          pivot_wider(names_from = SAMPLE, values_from = TAXID) %>%
          group_by(NAME) %>%
          summarise(across(16:42, ~ paste(., collapse = ", "))) %>%
          mutate(across(2:28, ~ str_replace_all(., "(NA,|,NA|,NA,|NA|,| )", "")))
        
        # Assign the values of col1 as row names and handle empty strings
        processed_mgs2arg_data <- df_wide %>%
          remove_rownames() %>% 
          column_to_rownames(var = "NAME") %>%
          mutate(across(everything(), ~ na_if(., "")))
        
        # Remove rows where all values are NA
        sample2arg_data <- processed_mgs2arg_data %>%
          filter(rowSums(is.na(.)) < ncol(.))
        
        # Convert from a dataframe to a matrix and replace NA with 0
        sample2arg_data_matrix <- data.matrix(sample2arg_data, rownames.force = NA) %>%
          replace(is.na(.), 0)
        
        # Count NA in data sets
        sum(is.na(sample2arg_data_matrix))
        
        # Convert to a binary matrix
        sample2arg_data_Bmatrix <- as.matrix((sample2arg_data_matrix > 0) + 0)
        
        # Order the matrix by row names
        sample2arg_data_Bmatrix <- sample2arg_data_Bmatrix[order(rownames(sample2arg_data_Bmatrix)), ]
        
        # Rearranging a Binary Matrix by Row Count
        rearrange_matrix <- function(mat) {
          # Calculate the row sums (count of '1's)
          ordered_indices <- order(rowSums(mat), decreasing = TRUE)
          
          # Rearrange the matrix
          return(mat[ordered_indices, ])
        }
        
        # Rearranging the binary matrix
        sorted_matrix <- rearrange_matrix(sample2arg_data_Bmatrix)
        
        # Modify ordering of the clusters using clustering callback option
        callback = function(hc, mat){
          sv = svd(t(mat))$v[,1]
          dend = reorder(as.dendrogram(hc), wts = sv)
          as.hclust(dend)
        }
        
        # Annotations col names
        # Transpose sample2arg_data
        sample2arg_data_transposed <- as.data.frame(t(sample2arg_data))
        
        # Convert rownames into column and select sample
        sample2arg_data_transposed <- sample2arg_data_transposed %>%
          mutate(SAMPLE = rownames(.)) %>%
          select(SAMPLE)
        
        # Covert to a dataframe and trim whitespace
        MetadataLocations <- data.frame(metadata())
        MetadataLocations$sample <- str_trim(MetadataLocations$SAMPLE)
        
        # Merge the metadata and sample dataset
        ann_col <- inner_join(sample2arg_data_transposed, MetadataLocations, by = "SAMPLE") %>%
          na.omit() %>%
          remove_rownames() %>%
          column_to_rownames(var = "SAMPLE") 
        
        # Generate a color palette based on the number of levels in location
        d_palette <- paletteer_d("ggsci::category10_d3", n = nlevels(as.factor(ann_col$LOCATION)))
        
        # Create a named dataframe for the col colors
        df4 <- ann_col %>%
          distinct(LOCATION) %>%
          mutate(color = d_palette)
        
        # Create a named vector for rows
        location_colors <- setNames(as.character(df4$color), df4$LOCATION)
        
        # Annotation rows
        # Pivot the Sample values into columns
        df_wide_ann_rows <- dataframe() %>%
          arrange(SAMPLE) %>%
          pivot_wider(names_from = DATABASE, values_from = c(GENE)) %>%
          group_by(NAME) %>%
          summarise(across(15:17, ~ paste(., collapse = ", "))) %>%
          mutate(across(2:4, ~ str_replace_all(., "(NA,|,NA|,NA,|NA| )", ""))) %>%
          na.omit() %>%
          remove_rownames() %>%
          column_to_rownames(var = "NAME") %>%
          mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))
        
        # Classify functions
        classify <- function(input_string, type) {
          if (is.na(input_string)) return(NA)
          elements <- unlist(strsplit(input_string, ","))
          unique_elements <- length(unique(elements))
          
          if (type == "virulence") {
            return(ifelse(unique_elements == 1, "Virulent", "Hypervirulent"))
          } else if (type == "amr") {
            return(ifelse(unique_elements == 1, "Drug-resistant", "Multidrug-resistant"))
          } else if (type == "MGE") {
            return(ifelse(unique_elements == 1, "Single MGE", "Multi-MGEs"))
          }
        }
        
        # Apply classification functions
        ann_rows <- df_wide_ann_rows %>%
          mutate(VFs = sapply(vfdb, classify, type = "virulence"),
                 AMR = sapply(card, classify, type = "amr"),
                 MGEs = sapply(plasmidfinder, classify, type = "MGE")) %>%
          select(VFs, AMR, MGEs)
        
        # Resistance category variable
        Resistance_category <- dataframe() %>%
          distinct(NAME, GENE, RESISTANCE) %>%
          mutate(RESISTANCE = ifelse(nchar(RESISTANCE) > 20, "Multi-drug", RESISTANCE))
        
        # Modify strings function
        modify_strings <- function(df, column_name, target_strings, replacement_strings) {
          if (length(target_strings) != length(replacement_strings)) {
            stop("target_strings and replacement_strings must be of the same length")
          }
          
          for (i in seq_along(target_strings)) {
            df <- df %>%
              mutate(!!sym(column_name) := str_replace_all(!!sym(column_name), target_strings[i], replacement_strings[i]))
          }
          
          return(df)
        }
        
        # Call Modify String function
        Resistance_category <- modify_strings(Resistance_category, "RESISTANCE", 
                                              c("cephalosporin;penam", "carbapenem;penam", "Cephalosporin;Penem"), 
                                              c("Beta-lactam", "Beta-lactam", "Beta-lactam"))
        
        # Clean string function
        clean_string <- function(x) {
          str_to_title(str_replace_all(x, "_", " "))
        }
        
        # Apply the function to the RESISTANCE column
        Resistance_category <- Resistance_category %>%
          mutate(RESISTANCE = sapply(RESISTANCE, clean_string)) %>%
          distinct(NAME, .keep_all = TRUE) %>%
          column_to_rownames(var = "NAME") %>%
          select(RESISTANCE)
        
        # Merge all the annotation rows
        ann_row_all <- merge(ann_rows, Resistance_category, by = 0) %>%
          rename(DRUG = RESISTANCE) %>%
          column_to_rownames(var = "Row.names") %>%
          select(DRUG, AMR, VFs, MGEs)
        
        
        # Generate a color palette for drug categories
        r_palette <- paletteer_d("ggsci::default_igv", n = length(unique(ann_row_all$DRUG))-1)
        
        # Create a named dataframe for the Resistance category
        ann_row_drug_col <- ann_row_all %>%
          select(DRUG) %>%
          na.omit() %>%
          distinct(DRUG) %>%
          mutate(color = r_palette)
        
        # Create a named vector for rows
        drug_colors <- setNames(as.character(ann_row_drug_col$color), ann_row_drug_col$DRUG)
        
        # Annotation colors customization
        ann_colors <- list(
          LOCATION = location_colors,
          DRUG = drug_colors,
          AMR = c(`Drug-resistant` = "#a02d02", `Multidrug-resistant` = "#522501"),
          MGEs = c(`Multi-MGEs` = "#5c0404", `Single MGE` = "#fcd2d2"),
          VFs = c(`Virulent` = "#fce5d2", `Hypervirulent` = "#c47a02")
        )
        
        
        # Create heatmap using pheatmap package ##
       # pheatmap(sorted_matrix[input$row_count[1]:input$row_count[2], ], display_numbers = FALSE, cluster_cols = TRUE, cluster_rows = FALSE,
       #                           scale = "none", 
       #                           clustering_callback = callback,  
       #                           border_color = "NA", color = c("#CCCCCCFF", "#666666FF"),
       #                           legend_breaks = c(0, 1),
       #                           legend_labels = c("Absent", "Present"),
       #                           annotation_row = ann_row_all,
       #                           annotation_col = ann_col[4],
       #                           show_rownames = TRUE,
       #                           # cutree_cols = 15,
       #                           # cutree_rows = 20,
       #                           annotation_colors = ann_colors,
       #                           fontsize_row = 14,  # Adjust this value for row names
       #                           fontsize_col = 14,   # Adjust this value for column names
       #                           fontsize = 14  # Adjust this value for the legend text
       #                           
       #                           
       #  )
        
       # Attempt to create the heatmap
       tryCatch({
         pheatmap(sorted_matrix[input$row_count[1]:input$row_count[2], ], display_numbers = FALSE, cluster_cols = TRUE, cluster_rows = FALSE,
                  scale = "none", 
                  clustering_callback = callback,  
                  border_color = "NA", color = c("#CCCCCCFF", "#666666FF"),
                  legend_breaks = c(0, 1),
                  legend_labels = c("Absent", "Present"),
                  annotation_row = ann_row_all,
                  annotation_col = ann_col[4],
                  show_rownames = TRUE,
                  # cutree_cols = 15,
                  # cutree_rows = 20,
                  annotation_colors = ann_colors,
                  fontsize_row = 14,  # Adjust this value for row names
                  fontsize_col = 14,   # Adjust this value for column names
                  fontsize = 14  # Adjust this value for the legend text
                  
                  
         )
       }, error = function(e) {
         # Handle error: display a message in the console
         message("An error occurred while generating the heatmap: ", e$message)
         
       })
       
       
        # Plot Air microbiome pheatmap
        # heatmap_plot
        
      })

    
    # Return the reactive that yields the data frame
    return(filter_data)
  }
  
  )
}