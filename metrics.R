#### cross-validation MAE measure ####
cross_validation_MAE <- function(ensemble_model) {
    
    library(janitor)
    # ensemble_model <- greedy_ensemble
    
    list_carry_ens <- ensemble_model$ens_model$resample %>% as_tibble() %>% mutate(across(RMSE:MAE, ~signif(.x, 3)))
    model_list <- ensemble_model$models
    n_model <- length(names(model_list))
    name_model <- names(model_list)
    list_carry <- list()
    for (name in name_model) {
        list_carry[[name]] <- model_list[[name]]$resample %>% as_tibble() %>% mutate(across(RMSE:MAE, ~signif(.x, 3)))
    }
    
    list_carry[["ensemble_model"]] <- list_carry_ens
    
    whole_df <- bind_rows(list_carry, .id = "model")
    whole_df %>% 
        group_by(model) %>% 
        summarise(mean_MAE = mean(MAE)) %>% pull(mean_MAE)
    
    whole_df %>% 
        group_by(model) %>% 
        summarise(mean_MAE = formatC(signif(mean(MAE), 3), digits = 3, format = "fg", flag = "#"),
                  sd_MAE = formatC(signif(sd(MAE), 3), digits = 3, format = "fg", flag = "#"),
                  mean_RMSE = formatC(signif(mean(RMSE), 3), digits = 3, format = "fg", flag = "#"),
                  sd_RMSE = formatC(signif(sd(RMSE), 3), digits = 3, format = "fg", flag = "#"),
                  mean_Rsquared = formatC(signif(mean(Rsquared), 3), digits = 3, format = "fg", flag = "#"),
                  sd_Rsquared = formatC(signif(sd(Rsquared), 3), digits = 3, format = "fg", flag = "#")) %>% 
        arrange(factor(model, levels = name_model)) %>% 
        mutate(together_MAE = paste0(mean_MAE, " (",
                                     sd_MAE, ")"),
               
               together_RMSE = paste0(mean_RMSE, " (",
                                      sd_RMSE, ")"),
               
               together_Rsquared = paste0(mean_Rsquared, " (",
                                          sd_Rsquared, ")")
        ) -> summarised_results
    
    whole_df %>% 
        mutate(MAE = formatC(MAE, digits = 3, format = "fg", flag = "#"),
               RMSE = formatC(RMSE, digits = 3, format = "fg", flag = "#"),
               Rsquared = formatC(Rsquared, digits = 3, format = "fg", flag = "#")) %>% 
        select(-RMSE, -Rsquared) %>% 
        pivot_wider(names_from = model, values_from = MAE) %>% 
        separate_wider_delim(Resample, delim = ".", names = c("fold", "repetition")) %>% 
        mutate(Resample = paste0(repetition, ".", fold)) %>% 
        relocate(Resample) %>% 
        select(-fold, -repetition) -> cv_MAE_results
    
    cv_MAE_results %>% 
        bind_rows(summarised_results %>% select(model, together_MAE) %>% t() %>% as_tibble() %>% row_to_names(row_number = 1)) -> cv_MAE_results_final
    
    whole_df %>% 
        mutate(MAE = formatC(MAE, digits = 3, format = "fg", flag = "#"),
               RMSE = formatC(RMSE, digits = 3, format = "fg", flag = "#"),
               Rsquared = formatC(Rsquared, digits = 3, format = "fg", flag = "#")) %>% 
        select(-MAE, -Rsquared) %>% 
        pivot_wider(names_from = model, values_from = RMSE) %>% 
        separate_wider_delim(Resample, delim = ".", names = c("fold", "repetition")) %>% 
        mutate(Resample = paste0(repetition, ".", fold)) %>% 
        relocate(Resample) %>% 
        select(-fold, -repetition) %>% 
        mutate_all(as.character) -> cv_RMSE_results
    
    cv_RMSE_results %>% 
        bind_rows(summarised_results %>% select(model, together_RMSE) %>% t() %>% as_tibble() %>% row_to_names(row_number = 1)) -> cv_RMSE_results_final
    
    whole_df %>% 
        mutate(MAE = formatC(MAE, digits = 3, format = "fg", flag = "#"),
               RMSE = formatC(RMSE, digits = 3, format = "fg", flag = "#"),
               Rsquared = formatC(Rsquared, digits = 3, format = "fg", flag = "#")) %>% 
        select(-MAE, -RMSE) %>% 
        pivot_wider(names_from = model, values_from = Rsquared) %>% 
        separate_wider_delim(Resample, delim = ".", names = c("fold", "repetition")) %>% 
        mutate(Resample = paste0(repetition, ".", fold)) %>% 
        relocate(Resample) %>% 
        select(-fold, -repetition) %>% 
        mutate_all(as.character) -> cv_RS_results
    
    cv_RS_results %>% 
        bind_rows(summarised_results %>% select(model, together_Rsquared) %>% t() %>% as_tibble() %>% row_to_names(row_number = 1)) -> cv_RS_results_final
    
    list(MAE = cv_MAE_results_final, RMSE = cv_RMSE_results_final, Rsquared = cv_RS_results_final)
}



#### cross-validation whole system (not including binary metrics) ####
cross_validation_list_all <- function(model_list, threshold) {
    
    library(janitor)
    # ensemble_model <- greedy_ensemble
    
    list_carry_ens <- model_list$ens_model$resample %>% as_tibble() %>% mutate(across(ROC, ~signif(.x, 3)))
    model_list <- model_list$models
    n_model <- length(names(model_list))
    name_model <- names(model_list)
    list_carry <- list()
    for (name in name_model) {
        list_carry[[name]] <- model_list[[name]]$resample %>% as_tibble() %>% mutate(across(ROC, ~signif(.x, 3)))
    }
    
    list_carry[["ensemble_model"]] <- list_carry_ens
    
    whole_df <- bind_rows(list_carry, .id = "model")
    
    whole_df %>% 
        group_by(model) %>% 
        summarise(mean_ROC = formatC(signif(mean(ROC), 3), digits = 3, format = "fg", flag = "#"),
                  sd_ROC = formatC(signif(sd(ROC), 3), digits = 3, format = "fg", flag = "#")) %>% 
        arrange(factor(model, levels = name_model)) %>% 
        mutate(together_ROC = paste0(mean_ROC, " (",
                                     sd_ROC, ")")) -> summarised_results
    
    whole_df %>% 
        select(-Sens, -Spec) %>% 
        mutate(ROC = formatC(ROC, digits = 3, format = "fg", flag = "#")) %>% 
        pivot_wider(names_from = model, values_from = ROC) %>% 
        separate_wider_delim(Resample, delim = ".", names = c("fold", "repetition")) %>% 
        mutate(Resample = paste0(repetition, ".", fold)) %>% 
        relocate(Resample) %>% 
        select(-fold, -repetition) -> cv_ROC_results
    
    cv_ROC_results %>% 
        bind_rows(summarised_results %>% select(model, together_ROC) %>% t() %>% as_tibble() %>% row_to_names(row_number = 1)) -> cv_ROC_results_final
    
    cv_ROC_results_final
}



#### cross-validation warning system measure ####
cross_validation_warning <- function(model, thresholds) {
    
    # model <- warning_system_all
    # thresholds <- warning_thresholds_all
    
    library(janitor)
    original_data <- 
        model$trainingData %>% 
        rowid_to_column()
    
    model$pred %>%
        as_tibble() %>% 
        mutate(RIDAGEYR = original_data$RIDAGEYR[match(rowIndex, original_data$rowid)]) %>% 
        mutate(threshold = ifelse(RIDAGEYR < 65, thresholds[1], thresholds[2]),
               pred_youden = ifelse(X1 >= threshold, "X1", "X0"),
               pred_youden = factor(pred_youden, levels = c("X1", "X0"))) %>% 
        group_by(Resample) -> tmp_data_pred
    
    
    # model$pred %>%
    #     as_tibble() %>% 
    #     mutate(pred_youden = ifelse(X1 >= threshold, "X1", "X0"),
    #            pred_youden = factor(pred_youden, levels = c("X1", "X0"))) %>% 
    #     group_by(Resample) -> tmp_data_pred
    
    tmp_data_pred %>% 
        roc_auc(truth = obs, X1) %>% 
        bind_rows(tmp_data_pred %>% accuracy(truth = obs, pred_youden),
                  tmp_data_pred %>% sens(truth = obs, pred_youden),
                  tmp_data_pred %>% specificity(truth = obs, pred_youden)) -> whole_df
    

    
    whole_df %>% 
        pivot_wider(names_from = .metric, values_from = .estimate) %>% 
        separate_wider_delim(Resample, delim = ".", names = c("fold", "repetition")) %>% 
        mutate(Resample = paste0(repetition, ".", fold)) %>% 
        relocate(Resample) %>% 
        select(-fold, -repetition, -.estimator) %>% 
        arrange(Resample) %>% 
        mutate(across(roc_auc:specificity, ~ round(.x, 3))) -> whole_df_prep
        
    whole_df_prep %>% 
        summarise(mean_roc = mean(roc_auc),
                  sd_roc = sd(roc_auc),
                  mean_acc = mean(accuracy),
                  sd_acc = sd(accuracy),
                  mean_sens = mean(sens),
                  sd_sens = sd(sens),
                  mean_spec = mean(specificity),
                  sd_spec = sd(specificity)) %>% 
        mutate(across(mean_roc:sd_spec, ~ formatC(signif(.x, 3), digits = 3, format = "fg", flag = "#")),
               roc_auc = paste0(mean_roc, " (", sd_roc, ")"),
               accuracy = paste0(mean_acc, " (", sd_acc, ")"),
               sens = paste0(mean_sens, " (", sd_sens, ")"),
               specificity = paste0(mean_spec, " (", sd_spec, ")")) %>% 
        select(-c(mean_roc:sd_spec)) -> summarised_results
    
    whole_df_prep %>% 
        mutate(across(roc_auc:specificity, ~ formatC(signif(.x, 3), digits = 3, format = "fg", flag = "#"))) %>% 
        bind_rows(summarised_results) -> final_results
    
    final_results
}



library(ggsci)
library(rsample)

continuous_boot_plasma <- function(splits) {
    
    validation_pred <- analysis(splits)
    results <- postResample(pred = validation_pred$predict_test, obs = validation_pred$plasma_fasting)
    
    tibble(term = c("MAE", "RMSE", "Rsquared"),
           estimate = c(results[names(results) == "MAE"],
                        results[names(results) == "RMSE"],
                        results[names(results) == "Rsquared"])
    )
}

continuous_boot_glyco <- function(splits) {
    
    validation_pred <- analysis(splits)
    results <- postResample(pred = validation_pred$predict_test, obs = validation_pred$glycohemoglobin)
    
    tibble(term = c("MAE", "RMSE", "Rsquared"),
           estimate = c(results[names(results) == "MAE"],
                        results[names(results) == "RMSE"],
                        results[names(results) == "Rsquared"])
    )
}


multi_boot <- function(splits) {
    
    validation_pred <- analysis(splits)
    results <- postResample(pred = validation_pred$predict_test, obs = validation_pred$plasma_fasting)
    hard_results <- 
        validation_pred %>% yardstick::accuracy(truth = binarize_diab, estimate = binarize_pred) %>% 
        bind_rows(validation_pred %>% yardstick::f_meas(truth = binarize_diab, estimate = binarize_pred),
                  validation_pred %>% yardstick::sens(truth = binarize_diab, estimate = binarize_pred),
                  validation_pred %>% yardstick::spec(truth = binarize_diab, estimate = binarize_pred),
                  validation_pred %>% yardstick::ppv(truth = binarize_diab, estimate = binarize_pred),
                  validation_pred %>% yardstick::npv(truth = binarize_diab, estimate = binarize_pred))
    
    
    
    
    tibble(term = c("MAE", "RMSE", "Rsquared_perc", paste0(hard_results$.metric, "_perc")),
           estimate = c(results[names(results) == "MAE"],
                        results[names(results) == "RMSE"],
                        results[names(results) == "Rsquared"] * 100,
                        hard_results$.estimate * 100)
    )
}


multi_boot_warning <- function(splits) {
    
    validation_pred <- analysis(splits)
    # validation_pred <- test_data_warningSystem_wPred
    total_results <- 
        validation_pred %>% yardstick::roc_auc(truth = binarize_diab, predict_test) %>% 
        bind_rows(validation_pred %>% yardstick::accuracy(truth = binarize_diab, estimate = pred_youden),
                  validation_pred %>% yardstick::f_meas(truth = binarize_diab, estimate = pred_youden),
                  validation_pred %>% yardstick::sens(truth = binarize_diab, estimate = pred_youden),
                  validation_pred %>% yardstick::spec(truth = binarize_diab, estimate = pred_youden),
                  validation_pred %>% yardstick::ppv(truth = binarize_diab, estimate = pred_youden),
                  validation_pred %>% yardstick::npv(truth = binarize_diab, estimate = pred_youden),
                  validation_pred %>% yardstick::brier_class(truth = binarize_diab, predict_test) %>% mutate(.estimate = 1 - .estimate))
    
    
    tibble(term = paste0(total_results$.metric, "_perc"),
           estimate = total_results$.estimate * 100)
}



ROC_boot_all <- function(splits) {
    
    validation_pred <- analysis(splits)
    total_results <- 
        validation_pred %>% 
        dplyr::select(binarize_diab:ensemble) %>% 
        tidyr::gather(key, value, rf:ensemble) %>% 
        dplyr::group_by(key) %>% 
        yardstick::roc_auc(truth = binarize_diab, value)
    
    tibble(term = total_results$key,
           estimate = total_results$.estimate * 100)
}
