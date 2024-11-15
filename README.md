# Banff Automation System


# **Overview**
Banff Automation System aims to offer a fully consistent and comprehensive Banff classification regarding the kidney allograft rejection diagnosis. This might be able to reduce the intra- and inter- observer variability, provide educational values, offer the consistent therapeutic patient management. The system has been validated with its reclassification ability in multicenter cohorts and 2 clinical trials of adult and pediatric patients from Europe and North America. The system is available as an application at https://transplant-prediction-system.shinyapps.io/Banff_automation/. The application is saved at the shinyapps.io by RStudio, https://www.shinyapps.io/. The code algorithm of the application is shown here. The code is written in RStudio.

# **System requirement**
Hardware requirements:
The Banff Automation System does not require more than a standard computer.

# Software requirements

**OS Requirements**

The system is online. There is no OS requirement.

**Browser requirement**

The application has been tested on the following browsers:

-	Google Chrome: 102.0.5005.61 (x86_64)
-	Firefox: version 72.0.2 (64-bit)
- Safari: version 15.5
- Opera: version 88.0.4412.27 (x86_64)

**R package dependencies**

If you would like to try the algorithm on your own machine, you need to install few packages.
R version 3.5.1 (Feather Spray), RStudio version 1.4.1106

Packages:

-	Shiny
-	Shinythemes
-	tidyverse
-	dplyr
-	collapsibleTree
-	shinycssloaders
-	shinyWidgets
-	rmarkdown
-	knitr
-	readxl
-	kableExtra
-	stringr
-	shinyjs
-	writexl

# **Installation guide**
Please download:

-	app.R
-	test1.Rmd
-	colorLOCAL_updated.rds
-	banff_df3.rds
-	test_biopsies.xlsm
-	banff_ref_total.bibtex
-   test_biopsies.xlsm
-   test_biopsies.xlsx


Install R and RStudio to run the code.
Install all the required packages if needed. If R and RStudio are already installed, the whole installation time would not be greater than 10 minutes.

app.R file is the main code file for the Banff Automation System.
code_explanation.docx file is to help comprehend the app.R code structure block by block.
test1.Rmd file is to create a PDF report integrated in app.R file.
colorLOCAL_updated and banff_df3.rds files are for coloring scheme data for the visualization tree-diagram.
banff_ref_total.bibtex is the bibliography file integrated in the PDF generation.
test_biopsies.xlsx and .xlsm files are the minimal data to reproduce the system. 


**Instructions for use:**
When the app.R file is opened in the RStudio, you can click “Run App” to run the Banff Automation System on-machine. If you would like to input your own biopsy data to test the system, you might want to create a similar Excel file as “test_biopsies.xlsm”. Please use the same column names as the Excel file has.
The first “Previsualization” tab allows a single biopsy for a rejection diagnosis. The second tab allows either single or multiple biopsies as Excel input. You can get the output results in either PDF or Excel format. 
“test_biopsies.xlsm” file is a sample file, which you might want to try first before a usage.


For the complete reproducibility, you should do as follows:
1. Download all files.
2. Create a folder then put all the files into the folder.
3. Create a child folder in it and name it "www".
4. Please put test_biopsies.xlsx and test_biopsies.xlsm files in the "www" folder.
5. Install all necessary libraries.
6. Run app.R file in R Studio. - In case of difficulty of understanding of the code structure, please refer to code_explanation.docx file.
7. Run the Shiny application. This will open a separated screen or a web browser. --> This is the replication of Banff Automation System (https://transplant-prediction-system.shinyapps.io/Banff_automation/).
8. Go to the 2nd panel of the application.
9. Click "Browse" button and import the downloaded test_biopsies.xlsm or test_biopsies.xlsx file.
10. Select either PDF report or Excel report.
11. Finally, select the Generate Report button. The report will be downloaded.


**For the figures of the manuscript (the sankey diagrams and survival figure)**:
Please download code_sankey.R, code_survival.R, banff_discrepancy_adult.csv, banff_discrepancy_pediatric.csv, and Figure5_data.csv files.
The code_sankey.R and code_survival.R files include all source codes from importing required libraries to sankey diagrams and survival figure generations, respectively.
The three CSV files excluded patient information and included deidentified patient IDs to anonymize the data.

# License
This project is covered under the Creative Commons Attribute 4.0 International License.
