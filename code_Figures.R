## Figures made on  R version 4.3.2 ##

##### Install packages #####
# Skip this if you already have installed them
install.packages(c("tidyverse", "ggplot2", "ggExtra", "parallel", "doParallel", "caret", "Boruta",
                   "xgboost", "pROC", "yardstick",  "rsample", "kernelshap", "shapviz", "patchwork",
                   "formattable", "dcurves"))


#### Load basic packages and datasets ####

library(tidyverse)
library(ggplot2)
library(ggExtra)
library(formattable)
library(parallel); library(doParallel)
library(dcurves)


nhanes_all_prep <- read_rds("nhanes_all_prep.rds")
nhanes_all_prep_weights <- read_rds("nhanes_all_prep_weights.rds")
dictionary_nhanes <- read.csv("dictionary_nhanes.csv") %>% as_tibble()
source("metrics.R") # Metrics to perform bootstrap



#### Figure 1. Fasting plasma glucose vs glycohemoglobin ####
nhanes_all_prep %>% 
    mutate(quadrants = ifelse(plasma_fasting >= 100 & glycohemoglobin >= 5.7, "1st",
                              ifelse(plasma_fasting < 100 & glycohemoglobin >= 5.7, "2nd",
                                     ifelse(plasma_fasting < 100 & glycohemoglobin < 5.7, "3rd", "4th")))) %>% 
    pull(quadrants) %>% table() -> n_quadrants

nhanes_all_prep %>% 
    mutate(quadrants = ifelse(plasma_fasting >= 100 & glycohemoglobin >= 5.7, "1st",
                              ifelse(plasma_fasting < 100 & glycohemoglobin >= 5.7, "2nd",
                                     ifelse(plasma_fasting < 100 & glycohemoglobin < 5.7, "3rd", "4th")))) %>% 
    pull(quadrants) %>% table() %>% prop.table() -> percent_quadrants

quadrant_1st <- paste0("n=", comma(n_quadrants[1], format = "d"), ", ", 
                       formatC(signif(percent_quadrants[1] * 100, 3), digits = 3, format = "fg", flag = "#"), "%")
quadrant_2nd <- paste0("n=", comma(n_quadrants[2], format = "d"), ", ",
                       formatC(signif(percent_quadrants[2] * 100, 3), digits = 3, format = "fg", flag = "#"), "%")
quadrant_3rd <- paste0("n=", comma(n_quadrants[3], format = "d"), ", ",
                       formatC(signif(percent_quadrants[3] * 100, 3), digits = 3, format = "fg", flag = "#"), "%")
quadrant_4th <- paste0("n=", comma(n_quadrants[4], format = "d"), ", ",
                       formatC(signif(percent_quadrants[4] * 100, 3), digits = 3, format = "fg", flag = "#"), "%")
quadrant_all <- "Total of 17,458 participants"


# Weights
nhanes_all_prep_weights %>%
    mutate(quadrants = ifelse(plasma_fasting >= 100 & glycohemoglobin >= 5.7, "1st",
                              ifelse(plasma_fasting < 100 & glycohemoglobin >= 5.7, "2nd",
                                     ifelse(plasma_fasting < 100 & glycohemoglobin < 5.7, "3rd", "4th")))) %>%
    select(quadrants, adjusted_mec_weights) %>%
    group_by(quadrants) %>%
    summarise(sum = sum(adjusted_mec_weights)) %>%
    mutate(prop = sum / sum(sum) * 100) %>% 
    mutate(prop_fix = formatC(signif(prop, 3), digits = 3, format = "fg", flag = "#")) -> weighted_distribution

quadrant_1st <- paste0(quadrant_1st, "\n", "(weighted, ", weighted_distribution$prop_fix[1], "%)")
quadrant_2nd <- paste0(quadrant_2nd, "\n", "(weighted, ", weighted_distribution$prop_fix[2], "%)")
quadrant_3rd <- paste0(quadrant_3rd, "\n", "(weighted, ", weighted_distribution$prop_fix[3], "%)")
quadrant_4th <- paste0(quadrant_4th, "\n", "(weighted, ", weighted_distribution$prop_fix[4], "%)")




nhanes_all_prep %>% 
    ggplot(aes(x = plasma_fasting, y = glycohemoglobin)) +
    geom_point() +
    geom_hline(yintercept = 5.7, color = "cornflowerblue", linewidth = 1.5, linetype = "dotted") +
    geom_vline(xintercept = 100, color = "cornflowerblue", linewidth = 1.5, linetype = "dotted") + 
    theme_bw() +
    xlab("Fasting Plasma Glucose (mg/dL)") +
    ylab("Glycohemoglobin (%)") +
    annotate("text",
             x = c(570, 40, 40, 570, 550), 
             y = c(14, 14, 2, 2, 20), 
             label = c(quadrant_1st, quadrant_2nd, quadrant_3rd, quadrant_4th, quadrant_all), 
             size = c(4, 4, 4, 4, 5)) + 
    scale_x_continuous(breaks = seq(0, 800, 100), limits = c(0, NA)) +
    scale_y_continuous(breaks = seq(0, 20, 2), limits = c(0, NA)) +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 15)) -> fpg_vs_glyco

fpg_vs_glyco
ggMarginal(fpg_vs_glyco, type = "density", xparams = list(fill = "deepskyblue", color = "white"),
           yparams = list(fill = "tomato", color = "white")) -> figure1

figure1




#### Figure 2. Boruta selected exposome parameters ####
# Dichotomize the outcome
nhanes_all_prep %>% 
    mutate(binarize_diab = ifelse(plasma_fasting >= 100 | glycohemoglobin >= 5.7, "X1", "X0"),
           binarize_diab = factor(binarize_diab, levels = c("X1", "X0"))) %>% 
    select(-plasma_fasting,
           -glycohemoglobin) -> nhanes_all_prep

nhanes_all_prep$binarize_diab %>% table()
nhanes_all_prep$binarize_diab %>% table() %>% prop.table


# Data split
library(caret)
set.seed(123)
ind <- createDataPartition(nhanes_all_prep$binarize_diab, p = 0.7, list = F)
all_train_data <- nhanes_all_prep[ind, ]
all_test_data <- nhanes_all_prep[-ind, ]

# Boruta
library(Boruta)
library(xgboost)

# setting multicores
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

set.seed(123)
all_boruta_result <- Boruta(y = all_train_data$binarize_diab,
                            x = all_train_data %>% select(-SEQN, -SEQN_new, -SDDSRVYR, -binarize_diab),
                            maxRuns = 2000,
                            doTrace = 1,
                            getImp = getImpXgboost)

all_boruta_result_fixed <- TentativeRoughFix(all_boruta_result)

all_confirmed_features <- names(all_boruta_result_fixed$finalDecision)[all_boruta_result_fixed$finalDecision == "Confirmed"]

all_boruta_result_fixed$ImpHistory %>% 
    as_tibble() %>%
    select(all_of(all_confirmed_features)) %>% 
    gather(key = "parameters", value = "importance") -> all_boruta_result_detail

# Edit the parameters' names from code to description 
all_old_names_boruta <- unique(all_boruta_result_detail$parameters)

dictionary_nhanes %>% 
    filter(variable_codename_use %in% all_old_names_boruta) %>%
    distinct(variable_codename_use, .keep_all = T) %>% 
    mutate(variable_description_use = ifelse(variable_codename_use == "RIAGENDR", "Gender of the participant",
                                             ifelse(variable_codename_use == "RIDAGEYR", "Age in years of the participant at the time of screening",
                                                    ifelse(variable_codename_use == "BMXBMI", "Body Mass Index (kg/m²)", variable_description_use)))) -> all_new_names_boruta


all_boruta_result_detail %>% 
    mutate(parameters = all_new_names_boruta$variable_description_use[match(parameters, all_new_names_boruta$variable_codename_use)],
           parameters = fct_reorder(parameters, importance)) -> all_boruta_result_detail

all_boruta_result_detail %>% 
    ggplot(aes(x = parameters, y = importance)) +
    geom_violin(aes(color = parameters, fill = parameters), position = position_dodge(), alpha = 0.7, scale = "width") +
    xlab("Exposome") + ylab("Importance") +
    scale_y_continuous(breaks = seq(-2, 1, 0.1)) +
    scale_colour_viridis_d(option = "turbo") +
    scale_fill_viridis_d(option = "turbo") +
    coord_flip() + 
    theme_bw() +
    theme(
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"),
        legend.position = "none"
    ) -> figure2

figure2




#### Figure 3. MEDWACS' discrimination and calibration performance ####
warning_system_all <- read_rds("warning_system_all.rds") # MEDWACS
warning_system_all

# Find Youden's cutoff for binary classification and metrics
library(pROC); library(yardstick)

all_train_data_tmp <- all_train_data %>% rownames_to_column()


warning_system_all$pred %>% 
    as_tibble() %>% 
    mutate(RIDAGEYR = all_train_data_tmp$RIDAGEYR[match(rowIndex, all_train_data_tmp$rowname)]) %>% 
    filter(RIDAGEYR < 65) -> to_youden_young

warning_system_all$pred %>% 
    as_tibble() %>% 
    mutate(RIDAGEYR = all_train_data_tmp$RIDAGEYR[match(rowIndex, all_train_data_tmp$rowname)]) %>% 
    filter(RIDAGEYR >= 65) -> to_youden_old


to_youden_young %>% roc(obs, X1, levels = c("X1", "X0")) -> tmp_youden_obj_young_all
to_youden_old %>% roc(obs, X1, levels = c("X1", "X0")) -> tmp_youden_obj_old_all


warning_threshold_young_all <- coords(tmp_youden_obj_young_all, x = "best", best.method = "youden")
warning_threshold_old_all <- coords(tmp_youden_obj_old_all, x = "best", best.method = "youden")

warning_thresholds_all <- c(warning_threshold_young_all$threshold, warning_threshold_old_all$threshold)

# External validation
NHANES_2023_final <- read_rds("NHANES_2023_final.rds") # NHANES 2021-2023
knhanes_2023_final <- read_rds("knhanes_2023_final.rds") # KNHANES 2023

## KNHANES dataset is missing upper leg length. This should be imputed with linear regression; height is also included for better prediction accuracy
all_train_data %>% 
    select(RIDAGEYR, RIAGENDR, VNAVEBPXSY, BMXARMC, BMXWAIST, BMXLEG, BMXBMI, BMXHT) -> all_train_forBMXLEG_KNHANES


# Make Upper Leg Length regression
set.seed(123)
model_ctrl_upper_leg <- trainControl(method = "repeatedcv",
                                     repeats = 3,
                                     number = 10,
                                     allowParallel = T,
                                     savePredictions = "final",
                                     index = createMultiFolds(all_train_forBMXLEG_KNHANES$BMXLEG, k = 10, times = 3))


all_final_features <- warning_system_all$trainingData %>% select(-.outcome) %>% names()


warning_system_leg <-
    caret::train(
        BMXLEG ~ RIDAGEYR + RIAGENDR + VNAVEBPXSY + BMXARMC + BMXWAIST + BMXBMI + BMXHT,
        data = all_train_forBMXLEG_KNHANES %>% select(all_of(all_final_features), BMXHT),
        trControl = model_ctrl_upper_leg,
        preProc = c("center", "scale"),
        method = "lm",
        metric = "MAE"
    )

warning_system_leg
warning_system_leg %>% summary()

# Impute the upper leg length
knhanes_2023_final %>% 
    mutate(BMXLEG = predict(warning_system_leg, knhanes_2023_final)) -> knhanes_2023_final


## Prepare the predicted dataset (internal validation [test set])
all_test_data_wPred <- 
    all_test_data %>%
    mutate(predict_test = predict(warning_system_all, all_test_data, type = "prob")$X1,
           pred_raw = predict(warning_system_all, all_test_data),
           threshold = ifelse(RIDAGEYR < 65, warning_thresholds_all[1], warning_thresholds_all[2]),
           pred_youden = ifelse(predict_test >= threshold, "X1", "X0"),
           pred_youden = factor(pred_youden, levels = c("X1", "X0")))



#### Figure 3A. Discrimination of MEDWACS with bootstrap ####
library(rsample)
set.seed(27)
boots_system_test_all <- bootstraps(all_test_data_wPred, times = 1000, apparent = TRUE)
result_system_test_all <- boots_system_test_all %>% mutate(result = map(splits, multi_boot_warning))

result_detail_all <- 
    bind_rows(result_system_test_all$result[1:1000]) %>% 
    dplyr::rename(metric = term,
                  performance = estimate) %>% 
    mutate(metric = factor(metric, levels = c("roc_auc_perc",
                                              "accuracy_perc", "f_meas_perc",
                                              "sens_perc", "spec_perc",
                                              "ppv_perc", "npv_perc",
                                              "brier_class_perc"),
                           labels = c("100 X ROCAUC", 
                                      "Accuracy (%)", "F1 Score (%)",
                                      "Sensitivity (%)", "Specificity (%)",
                                      "PPV (%)", "NPV (%)",
                                      "1 - Brier Score (%)")))


## Prepare the predicted dataset (NHANES 2021-2023)
NHANES_2023_final_wPred <- 
    NHANES_2023_final %>%
    mutate(predict_test = predict(warning_system_all, NHANES_2023_final, type = "prob")$X1,
           pred_raw = predict(warning_system_all, NHANES_2023_final),
           threshold = ifelse(RIDAGEYR < 65, warning_thresholds_all[1], warning_thresholds_all[2]),
           pred_youden = ifelse(predict_test >= threshold, "X1", "X0"),
           pred_youden = factor(pred_youden, levels = c("X1", "X0")))




set.seed(27)
boots_system_external_NHANES23_all <- bootstraps(NHANES_2023_final_wPred, times = 1000, apparent = TRUE)
result_system_external_NHANES23_all <- boots_system_external_NHANES23_all %>% mutate(result = map(splits, multi_boot_warning))


result_detail_external_NHANES23_all <- 
    bind_rows(result_system_external_NHANES23_all$result[1:1000]) %>% 
    dplyr::rename(metric = term,
                  performance = estimate) %>% 
    mutate(metric = factor(metric, levels = c("roc_auc_perc",
                                              "accuracy_perc", "f_meas_perc",
                                              "sens_perc", "spec_perc",
                                              "ppv_perc", "npv_perc",
                                              "brier_class_perc"),
                           labels = c("100 X ROCAUC", 
                                      "Accuracy (%)", "F1 Score (%)",
                                      "Sensitivity (%)", "Specificity (%)",
                                      "PPV (%)", "NPV (%)",
                                      "1 - Brier Score (%)")))

## Prepare the predicted dataset (KNHANES 2023)
KNHANES_2023_final_wPred <- 
    knhanes_2023_final %>%
    mutate(predict_test = predict(warning_system_all, knhanes_2023_final, type = "prob")$X1,
           pred_raw = predict(warning_system_all, knhanes_2023_final),
           threshold = ifelse(RIDAGEYR < 65, warning_thresholds_all[1], warning_thresholds_all[2]),
           pred_youden = ifelse(predict_test >= threshold, "X1", "X0"),
           pred_youden = factor(pred_youden, levels = c("X1", "X0")))


set.seed(27)
boots_system_external_KNHANES23_all <- bootstraps(KNHANES_2023_final_wPred, times = 1000, apparent = TRUE)
result_system_external_KNHANES23_all <- boots_system_external_KNHANES23_all %>% mutate(result = map(splits, multi_boot_warning))


result_detail_external_KNHANES23_all <- 
    bind_rows(result_system_external_KNHANES23_all$result[1:1000]) %>% 
    dplyr::rename(metric = term,
                  performance = estimate) %>% 
    mutate(metric = factor(metric, levels = c("roc_auc_perc",
                                              "accuracy_perc", "f_meas_perc",
                                              "sens_perc", "spec_perc",
                                              "ppv_perc", "npv_perc",
                                              "brier_class_perc"),
                           labels = c("100 X ROCAUC", 
                                      "Accuracy (%)", "F1 Score (%)",
                                      "Sensitivity (%)", "Specificity (%)",
                                      "PPV (%)", "NPV (%)",
                                      "1 - Brier Score (%)")))

result_detail_all %>%
    mutate(validation = "Internal\n(test set)") %>% 
    bind_rows(result_detail_external_NHANES23_all %>% mutate(validation = "External #1\n(NHANES 2021-2023)")) %>% 
    bind_rows(result_detail_external_KNHANES23_all %>% mutate(validation = "External #2\n(KNHANES 2023)")) %>% 
    mutate(validation = factor(validation, levels = c("Internal\n(test set)", "External #1\n(NHANES 2021-2023)", "External #2\n(KNHANES 2023)"))) -> result_detail_int_ext_all


result_detail_int_ext_all %>%
    ggplot(aes(x = metric, y = performance)) +
    geom_violin(aes(color = validation, fill = validation), position = position_dodge(0.4), alpha = 0.7) +
    geom_boxplot(aes(color = validation, fill = validation), position = position_dodge(0.4), width = 0.1, color = "black", alpha = 0.2, show.legend = F) + 
    xlab("Metric") + ylab("Performance (%)") +
    scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
    scale_colour_viridis_d(option = "turbo") +
    scale_fill_viridis_d(option = "turbo") +
    guides(color = "none",
           fill = guide_legend(title = "Validation")) + 
    theme_bw() +
    theme(
        axis.text.x = element_text(size = 17),
        axis.title.x = element_text(size = 22, vjust = 0, face = "bold"),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_text(size = 22, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 17)
    ) -> figure3A_discrimination

figure3A_discrimination



#### Figure 3B. Calibration of MEDWACS ####
calPlotData_test_all <- calibration(binarize_diab ~ predict_test, all_test_data_wPred)
calPlotData_test_all <- calPlotData_test_all$data %>% as_tibble()
calPlotData_test_all %>% 
    mutate(across(c(Percent:Upper, midpoint), ~ .x / 100)) -> calPlotData_test_all

calPlotData_NHANES2023_all <- calibration(binarize_diab ~ predict_test, NHANES_2023_final_wPred)
calPlotData_NHANES2023_all <- calPlotData_NHANES2023_all$data %>% as_tibble()
calPlotData_NHANES2023_all %>% 
    mutate(across(c(Percent:Upper, midpoint), ~ .x / 100)) -> calPlotData_NHANES2023_all

calPlotData_KNHANES2023_all <- calibration(binarize_diab ~ predict_test, KNHANES_2023_final_wPred)
calPlotData_KNHANES2023_all <- calPlotData_KNHANES2023_all$data %>% as_tibble()
calPlotData_KNHANES2023_all %>% 
    mutate(across(c(Percent:Upper, midpoint), ~ .x / 100)) -> calPlotData_KNHANES2023_all

calPlotData_test_all %>% 
    mutate(validation = "Internal\n(test set)") %>% 
    bind_rows(calPlotData_NHANES2023_all %>% mutate(validation = "External #1\n(NHANES 2021-2023)")) %>% 
    bind_rows(calPlotData_KNHANES2023_all %>% mutate(validation = "External #2\n(KNHANES 2023)")) %>% 
    mutate(validation = factor(validation, levels = c("Internal\n(test set)", "External #1\n(NHANES 2021-2023)", "External #2\n(KNHANES 2023)"))) -> calPlotData_all_int_ext



ggplot(calPlotData_all_int_ext, aes(x = midpoint, y = Percent, color = validation)) +
    geom_smooth(se = F, method = "loess", linewidth = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    theme_bw() +
    labs(color = "Validation") + 
    scale_x_continuous(limits = c(-0.1, 1.1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(-0.1, 1.1), breaks = seq(0, 1, 0.2)) +
    scale_colour_viridis_d(option = "turbo") +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20, face = "bold"),
          plot.title = element_text(size = 20, face = "bold"),
          legend.position = "top",
          legend.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15)) +
    xlab("Predicted probabilities") + ylab("Observed probabilities") -> figure3B_calibration

figure3B_calibration



#### Figure 4. MEDWACS interpretation ####
#### Figure 4A. MEDWACS interpretation with SHapley Additive exPlanations (SHAP) ####
library(kernelshap); library(shapviz)
# warning_system_all <- read_rds("warning_system_all.rds")
all_final_features <- warning_system_all$trainingData %>% select(-.outcome) %>% names()

# setting multicores
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

# sample 2000 participants from train data
set.seed(123)
all_sample_for_SHAP_system <- sample(1:nrow(all_train_data), 2000, replace = F)

all_train_data %>% 
    dplyr::select(all_of(all_final_features)) %>% 
    dplyr::slice(all_sample_for_SHAP_system) %>% 
    as.data.frame() -> all_train_data_for_SHAP_system

# sample 200 participants among the 2000 samples
set.seed(123)
all_sample_for_SHAP_system_bg <- sample(1:nrow(all_train_data_for_SHAP_system), 200, replace = F)

all_train_data_for_SHAP_system %>% 
    dplyr::slice(all_sample_for_SHAP_system_bg) %>% 
    as.data.frame() -> all_train_data_for_SHAP_system_bg


set.seed(1)
all_shap_result_system <- permshap(warning_system_all, 
                                   X = all_train_data_for_SHAP_system,
                                   bg_X = all_train_data_for_SHAP_system_bg,
                                   type = "prob",
                                   parallel = T,
                                   verbose = T)


all_sv_system <- shapviz(all_shap_result_system)

# Fixing the names for better readability in the figure
names(all_sv_system) <- c("prediabetes_diabetes", "Normal")
all_old_names_system <- colnames(all_sv_system$Normal)
dictionary_nhanes %>% 
    filter(variable_codename_use %in% all_old_names_system) %>%
    distinct(variable_codename_use, .keep_all = T) %>% 
    mutate(variable_description_use = ifelse(variable_codename_use == "RIDAGEYR", "Age in years of the participant at the time of screening",
                                             ifelse(variable_codename_use == "RIAGENDR", "Gender of the participant",
                                                    ifelse(variable_codename_use == "BMXBMI", "Body Mass Index (kg/m²)", variable_description_use)))) %>% 
    arrange(match(variable_codename_use, all_old_names_system)) %>% 
    pull(variable_description_use) -> all_new_names_system


colnames(all_sv_system$prediabetes_diabetes) <- all_new_names_system
colnames(all_sv_system$Normal) <- all_new_names_system


all_sv_system$Normal <- NULL
all_sv_imp_system <- sv_importance(all_sv_system, kind = "bee", max_display = 20, viridis_args = list(option = "H"))

all_sv_imp_system &
    theme_classic() &
    xlab("SHAP value (prediction probability)") &
    ggtitle("Prediabetes or Diabetes") &
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15, face = "bold", vjust = -1),
          legend.title = element_text(size = 15, face = "bold"),
          plot.title = element_text(size = 15, face = "bold")) -> figure4A_SHAP

figure4A_SHAP




#### Figure 4B. MEDWACS interpretation with combined-effect of parameters ####
library(patchwork)

# warning_system_all <- read_rds("warning_system_all.rds")
warning_system_bmi <- read_rds("warning_system_bmi.rds")
warning_system_bmi


# Set the value grids
RIDAGEYR_grid <- seq(from = quantile(nhanes_all_prep$RIDAGEYR, .05), 
                     to = quantile(nhanes_all_prep$RIDAGEYR, .95),
                     length.out = 50)
BMXWAIST_grid <- seq(from = quantile(nhanes_all_prep$BMXWAIST, .05), 
                     to = quantile(nhanes_all_prep$BMXWAIST, .95), 
                     length.out = 50)

BMXARMC_list <- quantile(nhanes_all_prep$BMXARMC, c(0.05, 0.50, 0.95))
BMXLEG_list <- quantile(nhanes_all_prep$BMXLEG, c(0.05, 0.50, 0.95))
female <- "2"
male <- "1"
low_BP <- 110
high_BP <- 140

expand_grid(RIDAGEYR = RIDAGEYR_grid,
            BMXWAIST = BMXWAIST_grid,
            RIAGENDR = 1:2,
            VNAVEBPXSY = c(low_BP, high_BP),
            BMXARMC = BMXARMC_list,
            BMXLEG = BMXLEG_list) %>% 
    mutate(RIAGENDR = factor(RIAGENDR, levels = 1:2)) %>% 
    mutate(BMXBMI = predict(warning_system_bmi, .)) %>% # BMXBMI is regressed
    mutate(predict_test = predict(warning_system_all, ., type = "prob")$X1) -> data_fig_simulation



plot_list_FLS <- list(); plot_list_FHS <- list(); plot_list_MLS <- list(); plot_list_MHS <- list()
quadrants <- c("FLS", "FHS", "MLS", "MHS")


for (quad in quadrants) {
    if (quad == "FLS") {
        data_fig_simulation %>% 
            filter(RIAGENDR == female & VNAVEBPXSY == low_BP) -> tmp_data
    } else if (quad == "FHS") {
        data_fig_simulation %>% 
            filter(RIAGENDR == female & VNAVEBPXSY == high_BP) -> tmp_data
    } else if (quad == "MLS") {
        data_fig_simulation %>% 
            filter(RIAGENDR == male & VNAVEBPXSY == low_BP) -> tmp_data
    } else {
        data_fig_simulation %>% 
            filter(RIAGENDR == male & VNAVEBPXSY == high_BP) -> tmp_data
    }
    
    num_cell <- 1
    
    
    for (i in 1:3) {
        
        if (i == 1) {
            arm <- BMXARMC_list[1]
        } else if (i == 2) {
            arm <- BMXARMC_list[2]
        } else {
            arm <- BMXARMC_list[3]
        }
        
        
        for (j in 1:3) {
            
            if (j == 1) {
                leg <- BMXLEG_list[1]
            } else if (j == 2) {
                leg <- BMXLEG_list[2]
            } else {
                leg <- BMXLEG_list[3]
            }
            
            tmp_data %>%
                filter(BMXARMC == arm) %>%
                filter(BMXLEG == leg) -> tmp_forfig
            
            
            
            tmp_forfig %>% 
                ggplot(aes(x = RIDAGEYR, y = BMXWAIST, fill = predict_test)) +
                geom_raster(interpolate = T) + 
                scale_fill_gradient2(low = "green", 
                                     mid = "yellow", 
                                     high = "red",
                                     midpoint = 0.5,
                                     breaks = seq(0, 1, 0.2),
                                     limits = c(0, 1),
                                     na.value = NA) +
                guides(fill = guide_colourbar(theme = theme(
                    legend.key.height = unit(10, "cm"),
                    legend.key.width = unit(1, "cm")
                ))) +
                xlab("Age of the participant") + 
                ylab("Waist circumference (cm)") + 
                labs(fill = "Probability") +
                scale_x_continuous(breaks = seq(30, 100, 20)) +
                scale_y_continuous(breaks = seq(80, 120, 20)) +
                theme_minimal() +
                theme(axis.title = element_text(size = 20, face = "bold"),
                      axis.text = element_text(size = 13),
                      legend.title = element_text(size = 17, face = "bold", margin = ggplot2::margin(b = 15)),
                      legend.text = element_text(size = 13),
                      panel.grid = element_blank()) -> p_tmp
            
            
            if (quad == "FLS") {
                plot_list_FLS[[num_cell]] <- p_tmp
            } else if (quad == "FHS") {
                plot_list_FHS[[num_cell]] <- p_tmp
            } else if (quad == "MLS") {
                plot_list_MLS[[num_cell]] <- p_tmp
            } else {
                plot_list_MHS[[num_cell]] <- p_tmp
            }
            
            num_cell <- num_cell + 1
        }
    }
    
}



((wrap_plots(plot_list_FLS, byrow = T, axis_titles = "collect") | (ggplot() + theme_void()) | wrap_plots(plot_list_MLS, byrow = T, axis_titles = "collect")) + plot_layout(widths = c(4, 0.3, 4))) /
    (ggplot() + theme_void()) / 
    ((wrap_plots(plot_list_FHS, byrow = T, axis_titles = "collect") | (ggplot() + theme_void()) | wrap_plots(plot_list_MHS, byrow  = T, axis_titles = "collect")) + plot_layout(widths = c(4, 0.3, 4))) +
    plot_layout(guides = "collect", heights = c(4, 0.3, 4)) &
    theme(legend.title = element_text(vjust = 2),
          legend.box.spacing = unit(1, "cm")) -> figure4B_combined

figure4B_combined


#### Figure 5. Decision curve analysis of MEDWACS ####

dca(binarize_diab ~ USPSTF_2021 + ADA_2022 + predict_test, 
    data = all_test_data_wPred %>% 
        mutate(binarize_diab = fct_relevel(binarize_diab, "X0"),
               USPSTF_2021 = ifelse(RIDAGEYR >= 35 & RIDAGEYR <= 70 & BMXBMI >= 25, 1, 0),
               ADA_2022 = ifelse(RIDAGEYR >= 35, 1, 0)), 
    thresholds = seq(0, 1, by = 0.01),
    label = list(USPSTF_2021 = "USPSTF (2021)",
                 ADA_2022 = "ADA (2022)",
                 predict_test = "MEDWACS")) %>%
    plot(smooth = TRUE) +
    ggplot2::labs(x = "Screen Threshold Probability") +
    scale_y_continuous(breaks = seq(0, 0.7, 0.1)) + 
    # geom_line(linewidth = 2) +
    # scale_color_discrete(name = "Screen", labels = c("Screen All", "Screen None", "USPSTF (2021)", "ADA (2022)", "MEDWACS")) +
    scale_colour_viridis_d(option = "turbo",
                           name = "Screen",
                           labels = c("Screen All", "Screen None", "USPSTF (2021)", "ADA (2022)", "MEDWACS")) +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 15, face = "bold")) -> dcurve_result

dcurve_result






