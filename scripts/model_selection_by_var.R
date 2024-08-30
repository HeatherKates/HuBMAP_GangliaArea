analyze_dependent_variable <- function(dependent_variable) {
  
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(gridExtra)
  library(ggplot2)
  
  # Load and preprocess the data
  data <- read.csv("data/HuBMAP_ganglia_data.csv")
  data <- data %>% filter(!CaseID == "")
  colnames(data) <- gsub("\\.\\.","\\.", colnames(data))
  colnames(data) <- gsub("\\.\\.","\\.", colnames(data))
  
  # Summarize the data
  summary_table <- data %>%
    group_by(CaseID, Disease.Status) %>%
    summarize(
      n = n(),
      avg_value = mean(!!sym(dependent_variable)),
      unique_slide_age_count = n_distinct(Slide.Age)
    ) %>%
    ungroup()
  
  # Convert Slide.Age to a factor
  data$Slide.Age <- as.factor(data$Slide.Age)
  
  # Transform the dependent variable
  log_var_name <- paste0("log_", dependent_variable)
  data[[log_var_name]] <- log(data[[dependent_variable]] + 1)
  
  # Prepare the random effects
  random_effects <- c("Tissue.Region", "Donor.Age", "T1D_Duration_Binned", "Slide.Age", "Ganglia.Type", "CaseID")
  data$Tissue.Region <- as.factor(data$Tissue.Region)
  data$Ganglia.Type <- as.factor(data$Ganglia.Type)
  data$CaseID <- as.factor(data$CaseID)
  
  # Bin T1D.Duration into categories
  data$T1D_Duration_Binned <- gsub("^0$", "Onset", data$T1D.Duration)
  data$T1D_Duration_Binned[is.na(data$T1D_Duration_Binned)] <- "None"
  data$T1D_Duration_Binned[data$T1D.Duration > 0 & data$T1D.Duration < 11] <- "Short"
  data$T1D_Duration_Binned[data$T1D.Duration >= 11] <- "Long"
  data$T1D_Duration_Binned <- as.factor(data$T1D_Duration_Binned)
  
  # Filter out random effects with only one level
  random_effects <- random_effects[sapply(random_effects, function(effect) {
    if (is.factor(data[[effect]])) {
      length(levels(data[[effect]])) > 1
    } else {
      TRUE
    }
  })]
  
  # Create models using combinations of random effects
  models <- list()
  for (i in 1:length(random_effects)) {
    for (j in i:length(random_effects)) {
      if (i != j) {
        model_name <- paste0(random_effects[i], "_", random_effects[j])
        formula <- as.formula(paste(log_var_name, "~ Disease.Status + (1 |", random_effects[i], ") + (1 |", random_effects[j], ")"))
        models[[model_name]] <- lmer(formula, data = data)
      } else {
        model_name <- random_effects[i]
        formula <- as.formula(paste(log_var_name, "~ Disease.Status + (1 |", random_effects[i], ")"))
        models[[model_name]] <- lmer(formula, data = data)
      }
    }
  }
  
  # Add baseline models
  models$NoRandom <- lm(as.formula(paste(log_var_name, "~ Disease.Status")), data = data)
  models$CaseID_SlideAge <- lmer(as.formula(paste(log_var_name, "~ Disease.Status + (1 | CaseID) + (1 | Slide.Age)")), data = data)
  
  # Calculate AIC, BIC, logLik, and p-values
  aic_values <- sapply(models, AIC)
  bic_values <- sapply(models, BIC)
  logLik_values <- sapply(models, logLik)
  p_values <- sapply(models, function(model) {
    if (inherits(model, "lm")) {
      tidy(anova(model))$p.value[1]
    } else {
      tidy(anova(model))$p.value[1]
    }
  })
  
  # Rank the models
  aic_rank <- rank(aic_values, na.last = TRUE)
  bic_rank <- rank(bic_values, na.last = TRUE)
  logLik_rank <- rank(-logLik_values, na.last = TRUE)
  
  # Handle single effect models
  single_effect_models <- list()
  for (effect in random_effects) {
    formula <- as.formula(paste(log_var_name, "~", effect))
    single_effect_models[[effect]] <- lm(formula, data = data)
  }
  
  # Update results table with single effect models
  all_models <- c(models, single_effect_models)
  p_values <- sapply(all_models, function(model) {
    if (inherits(model, "lm")) {
      tidy(anova(model))$p.value[1]
    } else {
      tidy(anova(model))$p.value[1]
    }
  })
  
  # Extract model formulas and combine them into a single string
  model_formulas <- sapply(all_models, function(model) {
    paste(deparse(formula(model)), collapse = " ")
  })

  
  # Create results table
  results_table <- data.frame(
    `Dependent Variable` = rep(dependent_variable, length(all_models)),
    Model = model_formulas,
    Fixed_Effect = c(rep("Disease.Status", length(models)), random_effects),
    Random_Effects = c(sapply(names(models), function(name) {
      if (name == "NoRandom") {
        return("None")
      } else {
        return(gsub("_", ", ", name))
      }
    }), rep("None", length(single_effect_models))),
    Model_Fit_Rank = c(aic_rank, rep(NA, length(single_effect_models))),
    Model_Likelihood_Rank = c(logLik_rank, rep(NA, length(single_effect_models))),
    Significance_of_Fixed_Effect = round(p_values, 4)
  )
  
  # Write the results table to a CSV file
  write.csv(results_table, file = paste0("results/", dependent_variable, "_model_results.csv"), row.names = FALSE)
  
  # Generate and save plots
  effects <- c(random_effects, "Disease.Status")
  plot_list <- list()
  plots_per_page <- 4
  
  for (i in seq_along(effects)) {
    effect <- effects[i]
    if (is.numeric(data[[effect]])) {
      plot <- ggplot(data, aes_string(x = effect, y = log_var_name)) +
        geom_point() +
        geom_smooth(method = "lm") +
        labs(title = paste("Distribution of", dependent_variable, "by", effect), 
             x = effect, 
             y = paste("log(", log_var_name, " coverage)")) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    } else {
      plot <- ggplot(data, aes_string(x = effect, y = log_var_name, fill = effect)) +
        geom_boxplot() +
        labs(title = paste("Distribution of", dependent_variable, "by", effect), 
             x = effect, 
             y = paste("log(", log_var_name, " coverage)")) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_fill_discrete(name = effect)
    }
    
    plot_list[[i]] <- plot
    
    if (i %% plots_per_page == 0 || i == length(effects)) {
      file_num <- ceiling(i / plots_per_page)
      grid_plot <- do.call(grid.arrange, c(plot_list[((file_num - 1) * plots_per_page + 1):i], ncol = 2, nrow = 2))
      ggsave(filename = paste0("results/", dependent_variable, "_plots_page_", file_num, ".png"), plot = grid_plot, dpi = 300, width = 14, height = 10)
      plot_list <- list()
    }
  }
}
