
rm(list=ls())
# Load the required packages
library(dplyr)
library(remstats)
library(remstimate)
library(remify)
library(brms)
library(rstan)
library(bayesplot)
library(coda)
library(purrr)
library(tidyr)
library(ggplot2)
library(GGally)
library(projpred)
library(remulate)
library(devtools)
library(igraph)
library(loo)
library(stringr)
library(reshape2)
library(caret)
library(pROC)
library(glmnet)
library(gmodels)
library(tibble)
library(ggpubr)
library(cmdstanr)
library(posterior)

#Converter function from an array to a data frame
poisson_df1 <- function(events, tie_stats, tie_reh, t0 = 0) {
  # Get unique actors in the event data
  actors <- sort(unique(c(events$sender, events$receiver)))
  
  # Creating risk set
  risk_set <- vector("list", dim(tie_stats)[2])
  for (i in 1:dim(tie_stats)[2]) {
    risk_set[[i]] <- getDyad(tie_reh, i)
  }
  
  risk_set <- do.call(rbind, risk_set)[, 2:3]
  
  M <- nrow(events) # number of events
  poisson_list <- vector("list", M) # initialize list to store data frames
  
  for (m in 1:M) {
    # Calculate time difference between current and previous event
    t_curr <- events$time[m]
    t_prev <- if (m == 1)
      t0
    else
      events$time[m - 1]
    delta_m <- t_curr - t_prev
    
    # Get statistics for current event m
    stats_m <- as.data.frame(tie_stats[m, , ])
    
    # Combine risk set with covariates
    df_m <- cbind(risk_set, stats_m)
    
    # Create binary outcome y
    df_m$y <- ifelse(df_m$actor1 == events$sender[m] &
                       df_m$actor2 == events$receiver[m],
                     1,
                     0)
    
    # Add offset
    df_m$logDelta <- log(delta_m)
    
    # Store data frame for event m
    poisson_list[[m]] <- df_m
  }
  
  # Combine all event data frames in list into one data frame
  df_poisson <- do.call(rbind, poisson_list)
  df_poisson <- df_poisson %>%
    select(-actor1, -actor2, -baseline) %>%
    # turn all variables that include "same_" to factors for memory efficiency
    mutate(across(contains("same_"), as.factor)) %>%
    # turn y to integer for memory efficiency
    mutate(y = as.integer(y))
  
  return(df_poisson)
}

url <- "https://raw.githubusercontent.com/Ali9972Rostami/Bayesian-Penalisation-in-REMs-using-WTCPoliceCalls-dataset/main/data/UUsummerschool.rdata"
load(url(url))

names(WTCPoliceCalls)[names(WTCPoliceCalls) == "number"] <- "time"
names(WTCPoliceCalls)[names(WTCPoliceCalls) == "source"] <- "sender"
names(WTCPoliceCalls)[names(WTCPoliceCalls) == "recipient"] <- "receiver"  # Check spelling

dataset=WTCPoliceCalls

head(dataset)
dim(dataset)

set.seed(45)
# Define the actor attributes FIRST
info <- data.frame(
  name = 1:37,
  time = 0,
  age = sample(20:60, 37, replace = TRUE),
  gender = sample(c(0, 1), 37, replace = TRUE)
)

# Stratifying the age variable 
info$age_group <- cut(
  info$age,
  breaks = c(19, 35, 45, 60),  # note: left open
  labels = c(0, 1, 2),
  right = TRUE
)

# Convert to numeric (labels are factors)
info$age_group <- as.numeric(as.character(info$age_group))


eff1 <- ~ 
  # --- Memory-based effects ---
  inertia(scaling = "std") +
  reciprocity(scaling = "std") +
  isp(scaling = "std") +  
  osp(scaling = "std") +     
  itp(scaling = "std") + 
  otp(scaling = "std") +      
  
  # --- Degree-based effects ---
  indegreeSender(scaling = "std") +
  outdegreeReceiver(scaling = "std") +
  indegreeReceiver(scaling = "std") +
  outdegreeSender(scaling = "std") +
  totaldegreeSender(scaling = "std") +
  totaldegreeReceiver(scaling = "std") +
  
  # --- Time-based effects ---
  recencyReceiveReceiver() +
  recencyReceiveSender() +
  recencySendReceiver() +
  recencySendSender() +
  
  # --- Participation Shifts ---
  psABAB() +
  psABAY() +
  psABBA() +
  psABBY() +
  psABXA() +
  psABXB() +
  psABXY() +
  
  #exogenous effects
  send("age_group") +
  receive("age_group") +
  send("gender") +
  receive("gender") +
  
  #interactions
  inertia(scaling = "std"):reciprocity(scaling = "std") +
  isp(scaling = "std"):osp(scaling = "std") +
  itp(scaling = "std"):otp(scaling = "std")


split_ratio=0.8
n_train <- floor(nrow(dataset) * split_ratio)
train_events = dataset[1:n_train, ]
dim(train_events)
test_events  = dataset[(n_train + 1):nrow(dataset), ]
dim(test_events)

reh_train <- remify::remify(edgelist = train_events, model = "tie", directed = TRUE)
reh_test  <- remify::remify(edgelist = test_events, model = "tie", directed = TRUE)

stats_train <- remstats::remstats(reh = reh_train, tie_effects = eff1, attr_actors = info)
stats_test  <- remstats::remstats(reh = reh_test, tie_effects = eff1, attr_actors = info)

dim(stats_test)


# Fit the relational event model
rem_fit <- remstimate(reh = reh_train, stats = stats_train, method = "MLE", ncores = 1)

sink("outputs/rem_model_output.txt")
print(summary(rem_fit))
sink()


poisson_train <- poisson_df1(events = train_events, tie_stats = stats_train, tie_reh = reh_train)
dim(poisson_train)

poisson_test <- poisson_df1(events = test_events, tie_stats = stats_test, tie_reh = reh_test)
dim(poisson_test)



predictors <- names(poisson_train)[1:(which(colnames(poisson_train) == "y") - 1)]
glm_formula <- as.formula(paste("y ~", paste(
  c(predictors, "offset(logDelta)"), collapse = " + ")))

glm_fit <- glm(
  formula = glm_formula,
  data = poisson_train,
  family = poisson(link = "log"))

sink("outputs/glm_model_output.txt")
print(summary(glm_fit))
sink()

y_preds_glm <- predict(glm_fit, newdata = poisson_test, type = "link")
y_true_glm <- poisson_test$y

png("figures/ROC_Curve_glm.png", width = 1200, height = 800)
roc(y_true_glm, y_preds_glm, plot=TRUE, legacy.axes = TRUE, percent = TRUE,
    xlab = "False Positive Percentage", ylab = "True Positive Percentage", col = "#377eb8", lwd = 3, print.auc = TRUE)
legend("bottomright", legend = c("glm Model"))
dev.off()

dyads_per_event <- nrow(poisson_test) / nrow(test_events)  

poisson_test$event_id <- rep(1:nrow(test_events), each = dyads_per_event)
length(poisson_test$event_id)

poisson_test$hazard_glm <- predict(glm_fit, newdata = poisson_test, type = "response")  # exp(Xβ)
poisson_test$true_event <- poisson_test$y

head(poisson_test)


brm_fit <- brm(
  formula = y ~  inertia + reciprocity + isp + osp + itp + otp + indegreeSender + outdegreeReceiver + indegreeReceiver +
    outdegreeSender + totaldegreeSender + totaldegreeReceiver + recencyReceiveReceiver + recencyReceiveSender +
    recencySendReceiver + recencySendSender + psABAB + psABAY + psABBA + psABBY + psABXA + psABXB +
    psABXY + send_age_group + receive_age_group + send_gender + receive_gender + inertia:reciprocity +
    isp:osp + itp:otp + offset(logDelta),
  data = poisson_train,
  family = poisson(link = "log"),
  prior_strong <- set_prior(horseshoe(
    df         = 1,
    scale_global = 0.01,   # small ⇒ strong global shrinkage
    df_global  = 1,
    scale_slab = 0.1,
    df_slab    = 4,
    par_ratio  = 0.05,
    autoscale  = TRUE
  ), class = "b"),
  chains = 4,
  iter = 8000,
  warmup = 2000,
  cores = 4,       # parallel processing
  seed = 123,
  backend = "cmdstan"
)

sink("outputs/brm_model_output.txt")
print(summary(brm_fit), digits = 5)
sink()


# Get fitted values (posterior means of λ) from brm model
hazard_brm_post <- fitted(brm_fit, newdata = poisson_test, scale = "response")
str(hazard_brm_post)
head(hazard_brm_post)

if (length(dim(hazard_brm_post)) == 3) {
  hazard_brm <- colMeans(hazard_brm_post[, , 1])
} else if (length(dim(hazard_brm_post)) == 2) {
  hazard_brm <- colMeans(hazard_brm_post)
} else {
  stop("Unexpected dimensions in fitted() output")
}

poisson_test$hazard_brm <- hazard_brm

head(poisson_test)


preds_brm <- posterior_predict(brm_fit, newdata = poisson_test)
y_prob_brm <- colMeans(preds_brm)

#brm model- ROC curve
png("figures/ROC_Curve_brm.png", width = 1200, height = 800)
roc(poisson_test$y, y_prob_brm, plot=TRUE, legacy.axes = TRUE, percent = TRUE,
    xlab = "False Positive Percentage", ylab = "True Positive Percentage", col = "#377eb8", lwd = 3, print.auc = TRUE)
legend("bottomright", legend = c("brm Model"))
dev.off()


########## ggplot for the rem and glm models ##############

# ---- REM estimates & CIs ----
rem_est <- data.frame(
  term = names(rem_fit$coef),
  estimate = rem_fit$coef,
  se = sqrt(diag(rem_fit$vcov)),
  model = "REM Model"
) %>%
  filter(term != "baseline")


rem_est <- rem_est %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se
  )

glm_se <- sqrt(diag(vcov(glm_fit)))  # Standard errors
z <- 1.96  # For 95% CI

glm_est <- data.frame(
  term = names(coef(glm_fit)),
  estimate = coef(glm_fit),
  se = glm_se,
  lower = coef(glm_fit) - z * glm_se,
  upper = coef(glm_fit) + z * glm_se,
  model = "GLM Model"
)

# Drop intercept
glm_est <- glm_est %>% filter(term != "(Intercept)")

# ---- Combine both ----
all_estimates <- bind_rows(rem_est, glm_est)

p1 <- ggplot(all_estimates, aes(x = term, y = estimate, color = model, shape = model)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.1, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_shape_manual(values = c("REM Model" = 16, "GLM Model" = 17)) +
  scale_color_manual(values = c("REM Model" = "red", "GLM Model" = "deepskyblue3")) +
  labs(
    x = "Effect",
    y = "Estimate",
    color = "Model",
    shape = "Model"
  ) +
  theme(
    strip.text = element_text(size = 8, face = "bold"),  
    axis.title = element_text(size = 8, face = "bold"),  
    axis.text = element_text(size = 8, face = "bold"),  
    legend.title = element_text(size = 8, face = "bold"),  
    legend.text = element_text(size = 8, face = "bold"))


# Save as PNG
ggsave("figures/ggplot_rem_glm.png", plot = p1, width = 8, height = 6, dpi = 300)

########## ggplot for the rem, glm and brm models ###########

extract_glm_data <- function(model, model_name) {
  coefs <- coef(summary(model))
  tibble(
    Variable = rownames(coefs),
    Est = coefs[, "Estimate"],
    LB = Est - 1.96 * coefs[, "Std. Error"],
    UB = Est + 1.96 * coefs[, "Std. Error"],
    Model = model_name
  ) %>%
    filter(Variable != "(Intercept)")   # remove intercept here
}

extract_rem_data <- function(model, model_name) {
  tibble(
    Variable = names(model$coef),
    Est = model$coef,
    LB = model$coef - 1.96 * model$se,
    UB = model$coef + 1.96 * model$se,
    Model = model_name
  ) %>%
    filter(Variable != "baseline")      # remove baseline here
}

extract_brm_data <- function(model, model_name) {
  posterior_samples <- as.data.frame(model)
  
  # keep only coefficients starting with "b_" but exclude intercept
  posterior_effects <- posterior_samples %>%
    select(starts_with("b_")) %>%
    select(-b_Intercept)   # remove intercept safely
  
  # clean names
  names(posterior_effects) <- gsub("^b_", "", names(posterior_effects))
  
  posterior_effects %>%
    pivot_longer(cols = everything(), names_to = "Variable", values_to = "Estimate") %>%
    group_by(Variable) %>%
    summarise(
      Est = mean(Estimate),
      LB = quantile(Estimate, 0.025),
      UB = quantile(Estimate, 0.975),
      .groups = "drop"
    ) %>%
    mutate(Model = model_name)
}


df_glm <- extract_glm_data(glm_fit, "glm MLE")
df_mle <- extract_rem_data(rem_fit, "remstimate MLE")
df_brm_ <- extract_brm_data(brm_fit, "brm with Horseshoe Prior")

plot_df <- bind_rows(df_glm, df_mle, df_brm_)
plot_df$Model <- factor(plot_df$Model, levels = c("glm MLE", "remstimate MLE", "brm with Horseshoe Prior"))

pd <- position_dodge(width = 0.6)

p2 <- ggplot(plot_df, aes(x = Est, y = Variable, color = Model, shape = Model, linetype = Model)) +
  geom_errorbar(aes(xmin = LB, xmax = UB), position = pd, linewidth = 1) +
  geom_point(position = pd, size = 2.5) +
  geom_vline(xintercept = 0, color = "grey40", linetype = "dashed") +
  labs(x = "Coefficient Estimate", y = "Effect") +
  theme(
    strip.text = element_text(size = 8, face = "bold"),  
    axis.title = element_text(size = 8, face = "bold"),  
    axis.text = element_text(size = 8, face = "bold"),  
    legend.title = element_text(size = 8, face = "bold"),  
    legend.text = element_text(size = 8, face = "bold"))


# Save as PNG
ggsave("figures/ggplot_rem_glm_brm.png", plot = p2, width = 8, height = 6, dpi = 300)

########### Diagnostic Tests ##################

capture.output(summary(brm_fit)$fixed, file = "outputs/fixed_effects_summary.txt")

capture.output(neff_ratio(brm_fit), file = "outputs/neff_ratio.txt")


# Function to extract posterior summaries from your brm model
extract_posterior_data <- function(model, model_name) {
  posterior_samples <- as.data.frame(model)
  
  # Clean coefficient names (remove "b_" prefix)
  cleaned_names <- gsub("^b_", "", names(posterior_samples))
  names(posterior_samples) <- cleaned_names
  
  # Select predictors from your formula (exclude intercept)
  df <- posterior_samples %>%
    select(inertia, reciprocity, isp, osp, itp, otp, 
           indegreeSender, outdegreeReceiver, indegreeReceiver, outdegreeSender,
           totaldegreeSender, totaldegreeReceiver, recencyReceiveReceiver, recencyReceiveSender,
           recencySendReceiver, recencySendSender,
           psABAB, psABAY, psABBA, psABBY, psABXA, psABXB, psABXY,
           send_age_group, receive_age_group, send_gender, receive_gender,
           inertia:reciprocity, isp:osp, itp:otp) %>%
    pivot_longer(cols = everything(), names_to = "Variable", values_to = "Estimate") %>%
    group_by(Variable) %>%
    summarise(
      Est = mean(Estimate),
      LB = quantile(Estimate, 0.025),
      UB = quantile(Estimate, 0.975),
      .groups = "drop"
    ) %>%
    mutate(Model = model_name)
  
  return(df)
}

# Extract summaries for your brm_model
df_brm <- extract_posterior_data(brm_fit, "brm_model")

# Plot with ggplot2
p <- ggplot(df_brm, aes(x = Est, y = Variable)) +
  geom_errorbar(aes(xmin = LB, xmax = UB), width = 0.2, linewidth = 1, color = "steelblue") +
  geom_point(size = 2.5, color = "steelblue") +
  geom_vline(xintercept = 0, color = "grey40", linetype = "dashed") +
  labs(x = "Posterior Estimates with 95% Credible Intervals", y = "Predictor") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

# Save as PNG
ggsave("figures/ggplot_brm_model.png", plot = p, width = 8, height = 6, dpi = 300)



png("figures/mcmc_areas_plot.png", width = 1200, height = 800)
mcmc_areas(brm_fit, pars = c("b_inertia", "b_reciprocity", "b_isp", "b_osp", "b_itp", "b_otp", 
                             "b_indegreeSender", "b_outdegreeReceiver", "b_indegreeReceiver", "b_outdegreeSender",
                             "b_totaldegreeSender", "b_totaldegreeReceiver", "b_recencyReceiveReceiver", "b_recencyReceiveSender",
                             "b_recencySendReceiver", "b_recencySendSender",
                             "b_psABAB", "b_psABAY", "b_psABBA", "b_psABBY", "b_psABXA", "b_psABXB", "b_psABXY",
                             "b_send_age_group", "b_receive_age_group", "b_send_gender", "b_receive_gender", 
                             "b_inertia:reciprocity", "b_isp:osp", "b_itp:otp"), prob = 0.95, prob_outer = 1)
dev.off()


png("figures/pp_check_plot.png", width = 1200, height = 800)
pp_check(brm_fit)
dev.off()

# Extract as draws array
posterior1 <- as_draws_array(brm_fit)

# Get parameter names
param_names <- variables(posterior1)

fixed_params <- grep("^b_", param_names, value = TRUE)

trace_plot <- mcmc_trace(posterior1, pars = fixed_params) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", color = NA))
ggsave("figures/trace_plot.png", trace_plot, width = 12, height = 8, bg = "white")


#Density overlay plot
dens_plot <- mcmc_dens_overlay(posterior1, pars = fixed_params) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", color = NA))
ggsave("figures/density_overlay.png", dens_plot, width = 12, height = 8, bg = "white")


#Interval plot (similar to mcmc_areas)
area_plot <- mcmc_areas(posterior1, pars = fixed_params, prob = 0.95) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", color = NA))
ggsave("figures/posterior_intervals.png", area_plot, width = 12, height = 8, bg = "white")


########## The CI for 30%, 50%, 90% and 95% and their variable selection
# Step 1: Extract posterior draws
draws <- as_draws_df(brm_fit)
fixed_effects_names <- grep("^b_", colnames(draws), value = TRUE)

# Step 2: Function to get CI tibble for a given level
get_ci <- function(level) {
  lower <- (1 - level) / 2
  upper <- 1 - lower
  summarise_draws(draws[, fixed_effects_names],
                  ~quantile2(.x, probs = c(lower, upper))) %>%
    rename_with(~paste0("CI", level*100, "_", c("lower", "upper")), -variable)
}

# Step 3: Compute CIs for 30%, 50%, 90%, and 95%
ci30 <- get_ci(0.30)
ci50 <- get_ci(0.50)
ci90 <- get_ci(0.90)
ci95 <- get_ci(0.95)

# Step 4: Combine all CIs
combined_cis <- ci95 %>%
  inner_join(ci90, by = "variable") %>%
  inner_join(ci50, by = "variable") %>%
  inner_join(ci30, by = "variable")

sink("outputs/CIs_outputs.txt")
print(combined_cis, n=Inf)
sink()


# Step 5: Variable selection function
select_vars <- function(ci_df, level) {
  lower_col <- paste0("CI", level, "_lower")
  upper_col <- paste0("CI", level, "_upper")
  ci_df %>%
    filter(!(!!sym(lower_col) < 0 & !!sym(upper_col) > 0)) %>%
    pull(variable)
}

# Step 6: Run variable selection for each CI width
selected_vars_95 <- select_vars(combined_cis, 95)
selected_vars_90 <- select_vars(combined_cis, 90)
selected_vars_50 <- select_vars(combined_cis, 50)
selected_vars_30 <- select_vars(combined_cis, 30)

# Step 7: Output results
Selected_Variables= list(
  CI95_selected = selected_vars_95,
  CI90_selected = selected_vars_90,
  CI50_selected = selected_vars_50,
  CI30_selected = selected_vars_30
)

sink("outputs/Selected_Variables_outputs.txt")
print(Selected_Variables)
sink()



