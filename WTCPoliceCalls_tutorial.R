
#################### The coding of figure 1 in section 2 ###############################

rm(list = ls())

# Load required libraries
library(ggplot2)
library(reshape2)

# Common beta grid
beta <- seq(-3, 3, length.out = 500)

# Define prior scenarios
scenarios <- data.frame(
  prior_sd = c(1, 0.3),
  label = c("Standard Ridge Prior (N(0,1))",
            "Stronger Ridge Prior (N(0,0.3²))")
)

# Initialize empty data frame to store results
df_all <- data.frame()

# Loop through scenarios to compute Prior, Likelihood, Posterior
for (i in 1:nrow(scenarios)) {
  prior_sd <- scenarios$prior_sd[i]
  label <- scenarios$label[i]
  
  # Prior
  prior <- dnorm(beta, mean = 0, sd = prior_sd)
  
  # Likelihood
  like_mean <- 1.2
  like_sd <- 0.5
  likelihood <- dnorm(beta, mean = like_mean, sd = like_sd)
  
  # Posterior (normal-normal conjugacy)
  prior_var <- prior_sd^2
  like_var <- like_sd^2
  post_var <- 1 / (1/prior_var + 1/like_var)
  post_mean <- post_var * (0/prior_var + like_mean/like_var)
  posterior <- dnorm(beta, mean = post_mean, sd = sqrt(post_var))
  
  # Combine into a data frame
  df <- data.frame(
    beta = beta,
    Prior = prior,
    Likelihood = likelihood / max(likelihood),  # scaled for visibility
    Posterior = posterior,
    Scenario = label
  )
  
  df_all <- rbind(df_all, df)
}

# Reshape for ggplot
df_long <- melt(df_all, id.vars = c("beta", "Scenario"))

# Plot with updated Likelihood color and line type
ggplot(df_long, aes(x = beta, y = value, color = variable, linetype = variable)) +
  geom_line(size = 2) +  # Bold lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_vline(xintercept = 1.2, linetype = "dashed", color = "red", alpha = 0.5) +
  facet_grid(. ~ Scenario) +
  labs(
    x = expression(beta),
    y = "Density",
    color = "Distribution",
    linetype = "Distribution"
  ) +
  theme(
    panel.background = element_rect(fill = "gray90", color = NA),  # light gray panel
    plot.background = element_rect(fill = "gray90", color = NA),   # gray plot background
    strip.background = element_rect(fill = "gray80", color = NA),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  scale_color_manual(values = c("blue", "red", "darkgreen")) +  # Likelihood is now red
  scale_linetype_manual(values = c("dashed", "dashed", "solid"))  # Likelihood dashed like Prior


##################### The coding of section 3 and the relevant figures ##################################

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


### ggplot for the rem and glm models ####

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

## ggplot for the rem, glm and brm models ##

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

## Diagnostic Tests ##

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


#### The CI for 30%, 50%, 90% and 95% and their variable selection
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

###### Shrinkem Coding ######

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
library(shrinkem)

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


y_preds_glm <- predict(glm_fit, newdata = poisson_test, type = "link")
y_true_glm <- poisson_test$y

roc(y_true_glm, y_preds_glm, plot=TRUE, legacy.axes = TRUE, percent = TRUE,
    xlab = "False Positive Percentage", ylab = "True Positive Percentage", col = "#377eb8", lwd = 3, print.auc = TRUE)
legend("bottomright", legend = c("glm Model"))


theta_hat <- coef(rem_fit)  # MLE estimates
sigma_hat <- rem_fit$vcov  # Variance-covariance matrix


#Apply shrinkage using shrinkem
#Horseshoe shrinkage prior
shrinkem_horseshoe_strong <- shrinkem(
  x       = theta_hat,
  Sigma   = sigma_hat,
  type    = "horseshoe",
  df1     = 1,      # default
  df2     = 1,      # default
  scale2  = 0.01,   # smaller = stronger shrinkage
  iterations = 50000,
  burnin     = 5000,
  cred.level = 0.95
)

summary(shrinkem_horseshoe_strong)


beta_hs_shrinkem <- as.numeric(shrinkem_horseshoe_strong$estimates$shrunk.mean)
names(beta_hs_shrinkem) <- rownames(shrinkem_horseshoe_strong$estimates)  # keep names aligned  

X_test_shrinkem <- model.matrix(glm_formula, poisson_test)

# Rename the intercept column in X_test to match beta_hs
colnames(X_test_shrinkem)[colnames(X_test_shrinkem) == "(Intercept)"] <- "baseline"

# Keep only columns that match beta names
X_test_shrinkem <- X_test_shrinkem[, names(beta_hs_shrinkem)]

# Predictions (linear predictor)
y_preds_hs_shrinkem <- as.vector(X_test_shrinkem %*% beta_hs_shrinkem)
y_true_hs_shrinkem  <- poisson_test$y

# ROC curve
roc(y_true_hs_shrinkem, y_preds_hs_shrinkem,
    plot = TRUE, legacy.axes = TRUE, percent = TRUE,
    xlab = "False Positive Percentage",
    ylab = "True Positive Percentage",
    col = "#e41a1c", lwd = 3, print.auc = TRUE)
legend("bottomright", legend = c("Horseshoe Shrinkem Model"))



# -----------------------------
# 1. Prepare shrinkem data
# -----------------------------
shrinkem_df_horseshoe <- data.frame(
  parameter = c("inertia", "reciprocity", "isp", "osp", "itp", "otp", 
                "indegreeSender", "outdegreeReceiver", "indegreeReceiver", "outdegreeSender", 
                "totaldegreeSender", "totaldegreeReceiver", "recencyReceiveReceiver", 
                "recencyReceiveSender", "recencySendReceiver", "recencySendSender", 
                "psABAB", "psABAY", "psABBA", "psABBY", "psABXA", "psABXB", "psABXY", 
                "send_age_group", "receive_age_group", "send_gender", "receive_gender", 
                "inertia:reciprocity", "isp:osp", "itp:otp"),
  estimate = c(0.014, 0.027, -0.015, 0.039, 0.047, -0.008, 0.010, -0.150, 0.250,
               0.174, 0.118, 0.168, 0.092, 3.560, 6.998, 0.357, 0.920, -0.058, 2.229, -0.102,
               -0.832, 0.143, -0.250, 0.015, 0.038, 0.094, 0.004, -0.004, -0.002, 0.005),
  lower = c(-0.041, -0.022, -0.123, -0.078, -0.033, -0.106, -0.269, -0.436, -0.069,
            -0.042, -0.111, -0.196, -0.307, 0.414, 3.981, -0.216, -0.050, -0.688, 0.218, -1.061,
            -2.090, -0.364, -1.157, -0.093, -0.071, -0.074, -0.181, -0.008, -0.040, -0.029),
  upper = c(0.079, 0.092, 0.070, 0.188, 0.150, 0.084, 0.241, 0.032, 0.593,
            0.429, 0.533, 0.698, 0.896, 6.190, 10.410, 2.094, 2.068, 0.297, 3.819, 0.360,
            0.081, 0.931, 0.134, 0.143, 0.195, 0.381, 0.199, -0.001, 0.031, 0.045),
  model = "shrinkem"
)

# -----------------------------
# 2. Prepare brms data
# -----------------------------
brms_df_horseshoe <- data.frame(
  parameter = shrinkem_df_horseshoe$parameter,  # same order
  estimate = c(0.01156, 0.02297, -0.00250, 0.02870, 0.02961, -0.00874, 
               0.02474, -0.09175, 0.26778, 0.19632, 0.09776, 0.11495, 0.12550, 
               3.40619, 6.31846, 0.36240, 1.04022, -0.02676, 2.62210, -0.04472, 
               -0.51752, 0.24589, -0.16856, 0.01174, 0.02590, 0.06295, 0.00249, 
               -0.00378, -0.00328, 0.00809),
  lower = c(-0.03652, -0.02525, -0.08773, -0.06052, -0.03242, -0.08508,
            -0.16927, -0.36697, -0.03425, -0.02264, -0.08841, -0.15608, -0.18093,
            0.27546, 3.54761, -0.15743, -0.01139, -0.45725, 0.99588, -0.67696, 
            -1.70740, -0.12136, -1.10615, -0.07846, -0.06424, -0.06857, -0.14747,
            -0.00706, -0.03842, -0.01906),
  upper = c(0.06954, 0.08168, 0.06888, 0.15962, 0.12582, 0.06572, 0.23273,
            0.04872, 0.55236, 0.44144, 0.44188, 0.62038, 1.26090, 5.66337, 9.52706, 
            1.86922, 2.19072, 0.25424, 4.11948, 0.31031, 0.11376, 1.04441, 0.13555,
            0.12645, 0.17939, 0.31856, 0.16164, -0.00056, 0.01944, 0.04394),
  model = "brms"
)

# -----------------------------
# 3. Combine data
# -----------------------------
plot_df_horseshoe <- bind_rows(shrinkem_df_horseshoe, brms_df_horseshoe)

# Optional: order parameters by shrinkem estimate
plot_df_horseshoe$parameter <- factor(plot_df_horseshoe$parameter, levels = shrinkem_df_horseshoe$parameter)

# -----------------------------
# 4. Plot with ggplot2
# -----------------------------
ggplot(plot_df_horseshoe, aes(x = estimate, y = parameter, color = model)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 position = position_dodge(width = 0.6), height = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Coefficient Estimate", y = NULL) +
  theme(
    strip.text = element_text(size = 8, face = "bold"),  
    axis.title = element_text(size = 8, face = "bold"),  
    axis.text = element_text(size = 8, face = "bold"),  
    legend.title = element_text(size = 8, face = "bold"),  
    legend.text = element_text(size = 8, face = "bold"))




# Strong ridge shrinkage prior
shrinkem_ridge_strong <- shrinkem(
  x       = theta_hat,
  Sigma   = sigma_hat,
  type    = "ridge",
  scale2  = 10,   # match brms prior variance
  iterations = 80000,
  burnin     = 20000,
  cred.level = 0.95
)

# Print results
summary(shrinkem_ridge_strong)


beta_ridge_shrinkem <- as.numeric(shrinkem_ridge_strong$estimates$shrunk.mean)
names(beta_ridge_shrinkem) <- rownames(shrinkem_ridge_strong$estimates)  # keep names aligned  

X_test_shrinkem_ridge <- model.matrix(glm_formula, poisson_test)

# Rename the intercept column in X_test to match beta_hs
colnames(X_test_shrinkem_ridge)[colnames(X_test_shrinkem_ridge) == "(Intercept)"] <- "baseline"

# Keep only columns that match beta names
X_test_shrinkem_ridge <- X_test_shrinkem_ridge[, names(beta_ridge_shrinkem)]

# Predictions (linear predictor)
y_preds_ridge_shrinkem <- as.vector(X_test_shrinkem_ridge %*% beta_ridge_shrinkem)
y_true_ridge_shrinkem  <- poisson_test$y

# ROC curve
roc(y_true_ridge_shrinkem, y_preds_ridge_shrinkem,
    plot = TRUE, legacy.axes = TRUE, percent = TRUE,
    xlab = "False Positive Percentage",
    ylab = "True Positive Percentage",
    col = "#e41a1c", lwd = 3, print.auc = TRUE)
legend("bottomright", legend = c("Horseshoe Shrinkem Model"))


## --- shrinkem results (already in your summary) ---
shrinkem_df <- data.frame(
  term = c("inertia","reciprocity","isp","osp","itp","otp","indegreeSender",
           "outdegreeReceiver","indegreeReceiver","outdegreeSender",
           "totaldegreeSender","totaldegreeReceiver","recencyReceiveReceiver",
           "recencyReceiveSender","recencySendReceiver","recencySendSender",
           "psABAB","psABAY","psABBA","psABBY","psABXA","psABXB","psABXY",
           "send_age_group","receive_age_group","send_gender","receive_gender",
           "inertia:reciprocity","isp:osp","itp:otp"),
  estimate = c(-0.006,0.056,-0.084,-0.051,0.136,0.064,-0.093,-0.274,0.298,0.151,
               0.268,0.247,0.185,2.344,4.427,0.914,-0.353,-2.231,1.565,-1.910,
               -1.949,-1.352,-2.434,0.019,0.074,0.241,0.069,-0.005,-0.002,0.009),
  lower = c(-0.089,-0.025,-0.384,-0.473,-0.182,-0.281,-0.998,-0.577,-0.468,-0.253,
            -0.924,-0.729,-0.985,0.237,2.216,-0.560,-1.765,-3.632,-0.029,-3.319,
            -3.360,-2.640,-3.627,-0.151,-0.105,-0.073,-0.260,-0.009,-0.064,-0.054),
  upper = c(0.078,0.137,0.218,0.371,0.453,0.409,0.803,0.024,1.066,0.555,
            1.456,1.230,1.357,4.513,6.751,2.403,0.994,-0.859,3.028,-0.536,
            -0.597,-0.090,-1.238,0.190,0.252,0.559,0.400,-0.001,0.060,0.070)
) %>%
  mutate(model = "shrinkem")

## --- brms results (point estimates only) ---
brms_df <- data.frame(
  term = c("inertia","reciprocity","isp","osp","itp","otp","indegreeSender",
           "outdegreeReceiver","indegreeReceiver","outdegreeSender",
           "totaldegreeSender","totaldegreeReceiver","recencyReceiveReceiver",
           "recencyReceiveSender","recencySendReceiver","recencySendSender",
           "psABAB","psABAY","psABBA","psABBY","psABXA","psABXB","psABXY",
           "send_age_group","receive_age_group","send_gender","receive_gender",
           "inertia:reciprocity","isp:osp","itp:otp"),
  estimate = c(0.11552,0.15360,0.03382,0.06285,0.04074,0.01183,0.14909,
               0.13278,0.22696,0.19318,0.18335,0.19465,0.09850,0.21984,
               0.19760,0.06463,-0.02443,-0.06076,0.66839,-0.11090,-0.03194,
               -0.04818,-0.24779,0.01646,0.03669,0.02142,-0.03782,-0.00997,
               -0.03710,0.03222),
  lower = c(0.07516,0.12152,-0.06184,-0.02517,-0.16058,-0.21726,0.07274,
            0.05400,0.14697,0.10740,0.08490,0.10953,-0.02362,0.12543,-0.38526,
            -0.03915,-0.12813,-0.19583,0.56957,-0.19829,-0.16812,-0.17867,
            -0.35902,-0.06223,-0.05126,-0.08214,-0.14960,-0.01205,-0.08040,
            -0.00416),
  upper = c(0.15204,0.18369,0.37152,0.13991,0.10935,0.09053,0.23869,0.22644,
            0.30706,0.27993,0.25792,0.28187,0.19633,0.31371,0.32944,0.17821,
            0.07155,0.56048,0.77892,-0.02306,0.53767,0.52555,-0.16004,0.10871,
            0.11763,0.12278,0.06132,-0.00773,0.00136,0.07541)
) %>%
  mutate(model = "brms")

## --- combine ---
combined_df <- bind_rows(shrinkem_df, brms_df)

## --- plot ---
ggplot(combined_df, aes(x = estimate, y = term, color = model)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    position = position_dodge(width = 0.6), height = 0.3,
    na.rm = TRUE   # ignore brms rows with NA intervals
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Coefficient Estimate", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(size = 8, face = "bold"),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8, face = "bold")
  )




shrinkem_lasso_strong <- shrinkem(
  x       = theta_hat,
  Sigma   = sigma_hat,
  type    = "lasso",
  scale2  = 0.01,   # match brms prior variance
  iterations = 50000,
  burnin     = 5000,
  cred.level = 0.95
)

# Print results
summary(shrinkem_lasso_strong)


beta_lasso_shrinkem <- as.numeric(shrinkem_lasso_strong$estimates$shrunk.mean)
names(beta_lasso_shrinkem) <- rownames(shrinkem_lasso_strong$estimates)  # keep names aligned  

X_test_shrinkem_lasso <- model.matrix(glm_formula, poisson_test)

# Rename the intercept column in X_test to match beta_hs
colnames(X_test_shrinkem_lasso)[colnames(X_test_shrinkem_lasso) == "(Intercept)"] <- "baseline"

# Keep only columns that match beta names
X_test_shrinkem_lasso <- X_test_shrinkem_lasso[, names(beta_lasso_shrinkem)]

# Predictions (linear predictor)
y_preds_lasso_shrinkem <- as.vector(X_test_shrinkem_lasso %*% beta_lasso_shrinkem)
y_true_lasso_shrinkem  <- poisson_test$y

# ROC curve
roc(y_true_lasso_shrinkem, y_preds_lasso_shrinkem,
    plot = TRUE, legacy.axes = TRUE, percent = TRUE,
    xlab = "False Positive Percentage",
    ylab = "True Positive Percentage",
    col = "#e41a1c", lwd = 3, print.auc = TRUE)
legend("bottomright", legend = c("Horseshoe Shrinkem Model"))


