#Conclusion:
#  Best Model: The CaseID model is the most valid model, offering the best fit with the lowest AIC/BIC and no evidence that adding Slide.Age improves the model.
#Fixed Effects: The fixed effects (Disease.Status) are not significant in the CaseID model, indicating no strong evidence that disease status affects log_Ganglia_Area after accounting for CaseID.
#No Need to Consider the SlideAge Model: Given its worse fit and potential overfitting, the SlideAge model should not be preferred over the CaseID model, despite the significance of Other-Diabetes in that model.
#In summary, the CaseID model is preferred, the fixed effect of Disease.Status is not significant in this model, and the intercept reflects the baseline level of log_Ganglia_Area on the log scale.

library(dplyr)
data <- read.csv("data/HuBMAP_ganglia_data.csv")
data <- data %>% filter(!CaseID=="")
colnames(data) <- gsub("\\.\\.","\\.",colnames(data))
colnames(data) <- gsub("\\.\\.","\\.",colnames(data))
summary_table <- data %>%
  group_by(CaseID, Disease.Status) %>%
  summarize(
    n = n(),
    avg_ganglia=mean(Ganglia.Area.µm.2.),
    unique_slide_age_count = n_distinct(Slide.Age)
  ) %>%
  ungroup()
data$Slide.Age <- as.factor(data$Slide.Age)

#Due to very high variance in ganglia area calcualted from the raw data, data are transformed:
var(data$Ganglia.Area.µm.2.)
data$log_Ganglia_Area <- log(data$Ganglia.Area.µm.2. + 1) # Adding 1 to avoid log(0)

library(lme4)
library(lmerTest) # Provides p-values for fixed effects
#Model 1: No random effects
NoRandom <- lm(log_Ganglia_Area ~ Disease.Status,data=data)
#Model 2: Random intercept for CaseID only.
CaseID <- lmer(log_Ganglia_Area ~ Disease.Status + (1 | CaseID), data = data)
#Model 3: Random intercept for Slide.Age only.
SlideAge <- lmer(log_Ganglia_Area ~ Disease.Status + (1 | Slide.Age), data = data)
#Model 4: Random intercept for both CaseID and Slide.Age.
CaseID_SlideAge <- lmer(log_Ganglia_Area ~ Disease.Status + (1 | CaseID) + (1 | Slide.Age), data = data)

# Calculate AIC and BIC for each model
aic_values <- AIC(NoRandom, CaseID, SlideAge, CaseID_SlideAge)
bic_values <- BIC(NoRandom, CaseID, SlideAge, CaseID_SlideAge)

# Rank the models based on AIC and BIC
aic_rank <- rank(aic_values$AIC)
bic_rank <- rank(bic_values$BIC)

print(CaseID)
print(SlideAge)
print(CaseID_SlideAge)

data.frame(AIC(CaseID,SlideAge,CaseID_SlideAge)) %>% arrange(AIC)
data.frame(BIC(CaseID,SlideAge,CaseID_SlideAge)) %>% arrange(BIC)

anova(CaseID, CaseID_SlideAge) #insignificant, CaseID has lower AIC/BIC
anova(SlideAge, CaseID_SlideAge) #significant, CaseID_SlideAge has lower logLik (significant) CaseID_SlideAge has lower AIC, SlideAge has lower BIC

summary(NoRandom) # disease status is not significant, but it is 0.09 which is worth looking at. This is the diff between Other-Diabetes and Aab+ (Intercept)
summary(CaseID) # disease status is not significant
summary(SlideAge) # disease status is significant; This is the diff between Other-Diabetes and Aab+ (Intercept)
summary(CaseID_SlideAge) #boundary singular fit error

# Extract p-values for Disease.Status from the ANOVA output
tidy_anova <- function(model) {
  tidy(anova(model)) %>%
    kable() %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
}

#Run models to assess significance of random effects (as fixed effects) and add to results table
# Run the additional models
SlideAge_lm <- lm(log_Ganglia_Area ~ Slide.Age, data = data)
CaseID_lm <- lm(log_Ganglia_Area ~ CaseID, data = data)
# Extract logLik values for the original models (NA for the linear models)
logLik_values <- c(
  logLik(NoRandom), 
  logLik(CaseID), 
  logLik(SlideAge), 
  logLik(CaseID_SlideAge),
  NA,  # NA for SlideAge_lm
  NA   # NA for CaseID_lm
)

# Rank the models based on logLik (higher is better)
logLik_rank <- rank(-logLik_values, na.last = TRUE)

# Extract AIC values for the original models (NA for SlideAge_lm and CaseID_lm)
aic_values <- c(
  AIC(NoRandom), 
  AIC(CaseID), 
  AIC(SlideAge), 
  AIC(CaseID_SlideAge),
  NA,  # NA for SlideAge_lm
  NA   # NA for CaseID_lm
)

# Rank the models based on AIC (lower is better)
aic_rank <- rank(aic_values, na.last = TRUE)

# Extract p-values for the fixed effects from the ANOVA output
p_values <- c(
  tidy(anova(NoRandom))$p.value[1],
  tidy(anova(CaseID))$p.value[1],
  tidy(anova(SlideAge))$p.value[1],
  tidy(anova(CaseID_SlideAge))$p.value[1],
  tidy(anova(SlideAge_lm))$p.value[1],  # p-value for SlideAge_lm
  tidy(anova(CaseID_lm))$p.value[1]     # p-value for CaseID_lm
)
# Create the updated results table
results_table <- data.frame(
  Model = c("NoRandom", "CaseID", "SlideAge", "CaseID_SlideAge", "SlideAge_lm", "CaseID_lm"),
  Fixed_Effect = c(rep("Disease.Status", 4), "Slide.Age", "CaseID"),
  Random_Effects = c("None", "CaseID", "Slide.Age", "CaseID, Slide.Age", "None", "None"),
  Model_Fit_Rank = c(aic_rank[1:4], NA, NA),  # NA for the two added models
  Model_Likelihood_Rank = c(logLik_rank[1:4], NA, NA),  # NA for the two added models
  Significance_of_Fixed_Effect = p_values
)

#write results table
write.csv(results_table,file="results/model_results.csv",row.names = FALSE)

##Plots
# Create a boxplot of log_Ganglia_Area by CaseID
boxplot_plot <- ggplot(data, aes(x = CaseID, y = log_Ganglia_Area, fill = CaseID)) +
  geom_boxplot() +
  labs(title = "Distribution of log_Ganglia_Area by CaseID", x = "CaseID", y = "log(Ganglia Area)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels for readability
  scale_fill_discrete(name = "CaseID")
ggsave(plot = boxplot_plot, "results/area_by_CaseID_boxplot.png",dpi = 300,width = 10,height = 8)

# Create a boxplot of log_Ganglia_Area by CaseID
boxplot_slidAge <- ggplot(data, aes(x = as.factor(Slide.Age), y = log_Ganglia_Area, fill = as.factor(Slide.Age))) +
  geom_boxplot() +
  labs(title = "Distribution of log_Ganglia_Area by Slide Age", x = "Slide Age", y = "log(Ganglia Area)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels for readability
  scale_fill_discrete(name = "Slide Age")
ggsave(plot = boxplot_slidAge, "results/area_by_SlideAge_boxplot.png",dpi = 300,width = 10,height = 8)

# Create a boxplot of log_Ganglia_Area by CaseID
boxplot_slidAge <- ggplot(data, aes(x = as.factor(Slide.Age), y = log_Ganglia_Area, fill = as.factor(Slide.Age))) +
  geom_boxplot() +
  labs(title = "Distribution of log_Ganglia_Area by Slide Age", x = "Slide Age", y = "log(Ganglia Area)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels for readability
  scale_fill_discrete(name = "Slide Age")
ggsave(plot = boxplot_slidAge, "results/area_by_SlideAge_boxplot.png",dpi = 300,width = 10,height = 8)

boxplot_disease <- ggplot(data, aes(x = Disease.Status, y = log_Ganglia_Area, fill = Disease.Status)) +
  geom_boxplot() +
  labs(title = "Distribution of log_Ganglia_Area by Disease Status", x = " Disease Status", y = "log(Ganglia Area)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels for readability
  scale_fill_discrete(name = "Slide Age")
ggsave(plot = boxplot_slidAge, "results/area_by_Disease.Status_boxplot.png",dpi = 300,width = 10,height = 8)

