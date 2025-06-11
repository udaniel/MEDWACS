# Machineborne Early Diabetic Warning And Control System (MEDWACS)


# **Overview**
MEDWACS aims to offer a robust early warning system for prediabetes or diabetes, yet, with only with simple parameters. The system is available as an application at https://dtu-quantitative-sustainability-assessment.shinyapps.io/MEDWACS/. The application is saved at the shinyapps.io by RStudio, https://www.shinyapps.io/. The code to reproduce the Figures is available here. The code is written in RStudio.

# **System requirement**
Hardware requirements:
MEDWACS does not require more than a standard computer.

**R package dependencies**

If you would like to reproduce the figures on your own machine, you need to install few packages.
R version 4.3.2 (Feather Spray), RStudio version 2024.9.0.375.

Packages:

-	tidyverse (version 2.0.0)
-	ggplot2 (version 3.5.1)
-	ggExtra (version 0.10.1)
-	parallel (version 4.3.2)
-	doParallel (version 1.0.17)
-	caret (version 6.0-94)
-	Boruta (version 8.0.0)
-	xgboost (version 1.7.8.1)
-	pROC (version 1.18.5)
-	yardstick (version 1.3.1)
-	rsample (version 1.2.1)
-	kernelshap (version 0.7.0)
-	shapviz (version 0.9.4)
-	patchwork (version 1.3.0) 
-	formattable (version 0.2.1)

# **Instructions guide**
Install R and RStudio to run the code.
Install all the required packages if needed. If R and RStudio are already installed, the whole installation time would not be greater than 10 minutes.

Please download all files in the repository.
nhanes_all_prep.rds file has the dataset for the final analysis.
metrics.R file includes metrics to help bootstrap.
warning_system_all.rds file is MEDWACS.

Follow code_Figures.R file to reproduce the figures.
