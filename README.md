# hosMicro
Welcome to our R Shiny application! This interactive and user-friendly web application is designed to explore the hospital microbiome and its genetic makeup in depth. Our tool provides a comprehensive platform for researchers, clinicians, and bioinformaticians to analyze and visualize metagenomics data with ease.

![](app/www/AMR_Pathogens.png)

## Key Features:
- Hospital Microbiome: Gain insights into the microbiome specific to hospital settings and its impact on patient health.
- Pathogen Identification: Identify and study various pathogens present in hospital environments.
- Antimicrobial Resistance: Investigate antimicrobial resistance patterns and their implications for public health.
- Virulence Factors: Analyze virulence factors to understand the mechanisms of pathogenicity.
- Mobile Genetic Elements: Explore the role of mobile genetic elements in the spread of resistance and virulence factors.

We hope this application will be a valuable resource for your research and help advance our understanding of the hospital microbiome and its genetic components.

## Quick start
If docker is available, pull the image and run the Shiny app using the following commands:
```Sh
docker pull julio92ont/hosmicro:1.1.2

docker run --rm -p 3838:3838 -v /root/shiny_save julio92ont/hosmicro:1.1.2 R -e "shiny::runApp('/root/shiny_save', host='0.0.0.0', port=3838)"
```

## Installation
If installing from the source, an  R version >= 4.1.0 with the corresponding packages listed below is required. Most of these are easy to install on a linux-based system:
```Sh
### Clone the repository and move to the hosMicro directory
git clone https://github.com/Julio92-C/hosMicro.git
cd hosMicro

### Install packages and dependencies
sudo apt update
sudo apt install r-base
R --version
R -e "install.packages(pkgs=c('shiny','shinydashboard', 'shinydashboardPlus', 'DT', 'dplyr', 'plotly', 'readr', 'ggplot2', 'scales', 'forcats','thematic', 'gtsummary', 'paletteer', 'reshape2', 'tidyr', 'VennDiagram', 'tidyverse', 'stringr', 'ggsignif', 'vegan', 'circlize', 'pheatmap'), repos='https://cran.rstudio.com/')"

### Run the shiny app
R -e "shiny::runApp('.', host='0.0.0.0', port=3838)"
```

## Input
In order to successfully run hosMicro, two specific input CSV files are required:

1. Merged Report CSV: This file should combine the kraken2 report with the Abricate summary report. Its purpose is to effectively link taxonomy with antimicrobial resistance (AMR), virulence factors (VFs), and mobile genetic elements (MGEs) profiles. This integration is crucial for an in-depth analysis of the hospital microbiome.

2. Project Metadata CSV: The second file must contain essential metadata about the project. This typically includes details such as sample identifiers, collection dates, locations, and other relevant descriptors that provide context to the metagenomics data.

For reference and guidance, two sample CSV files demonstrating the required format and content are available in the data repository.


## Tutorial
If you'd like to get to know the app functionality better, feel free to check out the tutorial linked below! It's a great way to learn about all its features.

[![App tutorial.](app/www/hosMicro_shinnyApp.png)](https://www.youtube.com/watch?v=9njf0_LXSOI)

## Acknowledgments:
We want to express our gratitude to our colleague, Dr. C. Oscar Previtali, for his support in curating the essential dataset for this application. We also want to thank our supervisors, Prof. Hermine V. Mkrtchyan, Dr. Piotr Cuber, and Dr. Raju Misra, for their invaluable contributions and guidance throughout the project. Your expertise has been crucial in enhancing our understanding and implementation of the metagenomics bioinformatics pipeline. We are truly thankful for everything you have done.

# Author:
Julio C. Ortega Cambara 

PhD C. Bioinformatics

School of Biomedical Sciences

University of West London

St Mary's Rd, London W5 5RF
  
