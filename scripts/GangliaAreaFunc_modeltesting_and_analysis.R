library(dplyr)
data <- read.csv("data/HuBMAP_ganglia_data.csv")
data <- data %>% filter(!CaseID=="")
colnames(data) <- gsub("\\.\\.","\\.",colnames(data))
colnames(data) <- gsub("\\.\\.","\\.",colnames(data))

dependent_variable <- "Ganglia.Area.µm.2."
summary_table <- data %>%
  group_by(CaseID, Disease.Status) %>%
  summarize(
    n = n(),
    avg_nNos.cvg=mean(!!sym(dependent_variable)),
    unique_slide_age_count = n_distinct(Slide.Age)
  ) %>%
  ungroup()
data$Slide.Age <- as.factor(data$Slide.Age)

#Data are transformed:
var(data[[dependent_variable]])
# Create the log-transformed variable dynamically
log_var_name <- paste0("log_", dependent_variable)
data[[log_var_name]] <- log(data[[dependent_variable]] + 1)

library(lme4)
library(lmerTest) # Provides p-values for fixed effects

fixed_effect <- data$Disease.Status

# Create models to test using all possible combinations of random effects, no interactions for now,
# keeping in mind the structure of the data which is
# that each CaseID only has one level of Slide.Age, T1D.Duration, Tissue.Region, and Disease.State
# Name models logically based on the random effects included in the model.

# These are the random effects
random_effects <- c("Tissue.Region", "Donor.Age", "T1D.Duration", "Slide.Age", "Ganglia.Type", "CaseID")

# Convert character variables to factors
data$Tissue.Region <- as.factor(data$Tissue.Region)
data$Ganglia.Type <- as.factor(data$Ganglia.Type)
data$CaseID <- as.factor(data$CaseID)

# Step 1: Create the initial binned variable
data$T1D_Duration_Binned <- gsub("^0$","Onset",data$T1D.Duration)
data$T1D_Duration_Binned[is.na(data$T1D_Duration_Binned)] <- "None"
data$T1D_Duration_Binned[data$T1D.Duration > 0 & data$T1D.Duration < 11] <- "Short"
data$T1D_Duration_Binned[data$T1D.Duration >= 11] <- "Long"
data$T1D_Duration_Binned <- as.factor(data$T1D_Duration_Binned)

# Keep Donor.Age and T1D.Duration as numeric if they are continuous
# random_effects list now distinguishes between factors and continuous variables
random_effects <- c("Tissue.Region", "Donor.Age", "T1D_Duration_Binned", "Slide.Age", "Ganglia.Type", "CaseID")

# Remove any factors with only one level
random_effects <- random_effects[sapply(random_effects, function(effect) {
  if (is.factor(data[[effect]])) {
    length(levels(data[[effect]])) > 1
  } else {
    TRUE  # Keep numeric random effects
  }
})]

# Model creation loop
models <- list()

# Generate models with combinations of random effects
for (i in 1:length(random_effects)) {
  for (j in i:length(random_effects)) {
    if (i != j) {
      model_name <- paste0(random_effects[i], "_", random_effects[j])
      formula <- as.formula(paste(paste0("log_",dependent_variable)," ~ Disease.Status + (1 |", random_effects[i], ") + (1 |", random_effects[j], ")"))
      models[[model_name]] <- lmer(formula, data = data)
    } else {
      model_name <- random_effects[i]
      formula <- as.formula(paste(paste0("log_",dependent_variable)," ~ Disease.Status + (1 |", random_effects[i], ")"))
      models[[model_name]] <- lmer(formula, data = data)
    }
  }
}

# Add baseline models
model_formula <- as.formula(paste(log_var_name, "~ Disease.Status"))
models$NoRandom <- lm(model_formula, data = data)

model_formula <- as.formula(paste(log_var_name, "~ Disease.Status + (1 | CaseID) + (1 | Slide.Age)"))
models$CaseID_SlideAge <- lmer(model_formula, data = data)

# Extract AIC and BIC values for each model
aic_values <- sapply(models, AIC)
bic_values <- sapply(models, BIC)

# Rank the models based on AIC and BIC
aic_rank <- rank(aic_values, na.last = TRUE)
bic_rank <- rank(bic_values, na.last = TRUE)

# Extract logLik values for each model
logLik_values <- sapply(models, logLik)
logLik_rank <- rank(-logLik_values, na.last = TRUE) # Use negative to rank from highest to lowest

# Extract p-values for the fixed effects from the ANOVA output
p_values <- sapply(models, function(model) {
  if (inherits(model, "lm")) {
    tidy(anova(model))$p.value[1]
  } else {
    tidy(anova(model))$p.value[1]
  }
})

# List of random effects and the fixed effect for plotting
effects <- c(random_effects, "Disease.Status")

# Initialize a list to store the models
single_effect_models <- list()

# Loop over each random effect and create a linear model
for (effect in random_effects) {
  # Dynamically create the formula
  formula <- as.formula(paste(log_var_name, "~", effect))
  
  # Fit the linear model and store it in the list
  single_effect_models[[effect]] <- lm(formula, data = data)
}
# Collect all models including previously defined ones
all_models <- c(models, single_effect_models)

# Extract AIC and BIC values for each model
aic_values <- sapply(models, function(model) {
  if (inherits(model, "lm") || inherits(model, "lmerMod")) {
    AIC(model)
  } else {
    NA
  }
})

bic_values <- sapply(models, function(model) {
  if (inherits(model, "lm") || inherits(model, "lmerMod")) {
    BIC(model)
  } else {
    NA
  }
})

# Rank the models based on AIC and BIC
aic_rank <- rank(aic_values, na.last = TRUE)
bic_rank <- rank(bic_values, na.last = TRUE)

# Extract logLik values for each model
logLik_values <- sapply(models, function(model) {
  if (inherits(model, "lm") || inherits(model, "lmerMod")) {
    logLik(model)
  } else {
    NA
  }
})

# Rank the models based on logLik (higher is better)
logLik_rank <- rank(-logLik_values, na.last = TRUE)

# Extract p-values for the fixed effects from the ANOVA output
p_values <- sapply(all_models, function(model) {
  if (inherits(model, "lm")) {
    tidy(anova(model))$p.value[1]
  } else {
    tidy(anova(model))$p.value[1]
  }
})

# Create the updated results table
results_table <- data.frame(
  Model = names(all_models),
  Fixed_Effect = c(rep("Disease.Status", length(models)), random_effects),
  Random_Effects = c(sapply(names(models), function(name) {
    if (name == "NoRandom") {
      return("None")
    } else {
      return(gsub("_", ", ", name))
    }
  }), rep("None", length(single_effect_models))),
  Model_Fit_Rank = c(aic_rank, rep(NA, length(single_effect_models))),  # Rank only for Disease.Status models
  Model_Likelihood_Rank = c(logLik_rank, rep(NA, length(single_effect_models))),  # Rank only for Disease.Status models
  Significance_of_Fixed_Effect = p_values
)
results_table$Model_Fit_Rank <- as.integer(results_table$Model_Fit_Rank)
results_table$Model_Likelihood_Rank <- as.integer(results_table$Model_Likelihood_Rank)
results_table$Significance_of_Fixed_Effect <- round(results_table$Significance_of_Fixed_Effect,4)

# Write the results table to a CSV file
results_table <- results_table %>% arrange(Model_Fit_Rank)
write.csv(results_table, file = paste0("results/",dependent_variable,"_model_results.csv"), row.names = FALSE)

# List of random effects and the fixed effect for plotting
effects <- c(random_effects, "Disease.Status")

library(gridExtra)
library(ggplot2)
library(gridExtra)

# Define how many plots per page/file
plots_per_page <- 4

# Generate and store plots in a list
plot_list <- list()

# Loop over each effect and generate the appropriate plot
for (i in seq_along(effects)) {
  effect <- effects[i]
  
  # Determine if the effect is continuous or categorical
  if (is.numeric(data[[effect]])) {
    # Create a line plot for continuous variables
    plot <- ggplot(data, aes_string(x = effect, y = log_var_name)) +
      geom_point() +
      geom_smooth(method = "lm") +
      labs(title = paste("Distribution of", dependent_variable, "by", effect), 
           x = effect, 
           y = paste("log(", log_var_name, " coverage)")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  } else {
    # Create a boxplot for categorical variables
    plot <- ggplot(data, aes_string(x = effect, y = log_var_name, fill = effect)) +
      geom_boxplot() +
      labs(title = paste("Distribution of", dependent_variable, "by", effect), 
           x = effect, 
           y = paste("log(", log_var_name, " coverage)")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels for readability
      scale_fill_discrete(name = effect)
  }
  
  plot_list[[i]] <- plot
  
  # Every `plots_per_page` plots, or at the last plot, arrange and save
  if (i %% plots_per_page == 0 || i == length(effects)) {
    # Determine the file number
    file_num <- ceiling(i / plots_per_page)
    
    # Arrange the plots in a grid
    grid_plot <- do.call(grid.arrange, c(plot_list[((file_num-1)*plots_per_page + 1):i], ncol=2, nrow=2))
    
    # Save the grid of plots
    ggsave(filename = paste0("results/", dependent_variable, "_plots_page_", file_num, ".png"), plot = grid_plot, dpi = 300, width = 14, height = 10)
    
    # Clear the plot list for the next set of plots
    plot_list <- list()
  }
}
