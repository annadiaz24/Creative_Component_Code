library(tidyverse)
library(scales)
library(emmeans)

# Load and prep the data
df <- read_csv("matrix_dataset.csv", show_col_types = FALSE) %>%
  # Filter out any infinite or NA values in our target columns
  filter(is.finite(eig_error_relative), is.finite(cond_input)) %>%
  rename(dim = size) %>%
  mutate(
    dim_fct = as.factor(dim),
    scale_fct = as.factor(scale),
    scale_label = paste0("Scale (s) = ", scale)
  )

# =======================================================================
# Simulation Results: Algorithm Error
# =======================================================================

p1 <- ggplot(df, aes(x = log10(cond_input), 
                     y = log10(eig_error_relative), 
                     color = dim_fct, shape = scale_fct)) +
  # position_jitterdodge adds randomness inside the dodged columns!
  geom_point(position = position_jitterdodge(jitter.width = 0.15, 
                                             dodge.width = 0.6), 
             size = 2.5, alpha = 0.6) +
  scale_color_viridis_d(option = "viridis") +
  # Format axes to 10^x
  scale_x_continuous(
    breaks = c(2, 4, 6, 8, 12),
    labels = function(x) parse(text = paste0("10^", x))
  ) +
  scale_y_continuous(labels = math_format(10^.x)) +
  labs(
    title = "Simulation Results",
    x = expression("Condition Number (" * kappa * ")"),
    y = expression(e[lambda*",rel"]),
    color = "Dimension (n)", 
    shape = "Scale (s)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

ggsave("R_error_main_plot.png", plot = p1, width = 10, height = 6, dpi = 300,
       bg = "white")
print(p1)

# =======================================================================
# Variability Analysis 
# =======================================================================

# Calculate the standard deviation for each group
df_std <- df %>%
  group_by(cond_input, dim, scale, dim_fct, scale_fct) %>%
  summarise(std_log_error = sd(log10(eig_error_relative)), .groups = 'drop')

p2 <- ggplot(df_std, aes(x = log10(cond_input), y = std_log_error, 
                         color = dim_fct, shape = scale_fct)) +
  # position_dodge neatly separates the shapes side-by-side
  geom_point(position = position_dodge(width = 0.5), size = 4, alpha = 0.8) +
  scale_color_viridis_d(option = "viridis") + 
  # Format the X-axis to show 10^x
 scale_x_continuous(     
   breaks = c(2, 4, 6, 8, 12),     
   labels = function(x) parse(text = paste0("10^", x))   
   ) +
  scale_y_continuous() +
  labs(
    title = "Within-Group Variability Across Condition Number",
    x = expression("Condition Number (" * kappa * ")"),
    y = expression(sigma),
    color = "Dimension (n)", 
    shape = "Scale (s)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

ggsave("R_error_std_plot.png", plot = p2, width = 10, height = 6, dpi = 300,
       bg = "white")
print(p2)

# =======================================================================
# Average Error: Mean of Relative Error 
# =======================================================================

# Calculate the averages for each group
df_mean <- df %>%
  group_by(cond_input, dim, scale, dim_fct, scale_fct) %>%
  summarise(mean_log_error = mean(log10(eig_error_relative)), .groups = 'drop')

p_mean <- ggplot(df_mean, aes(x = log10(cond_input), 
                              y = mean_log_error, 
                              color = dim_fct, shape = scale_fct)) +
  # Add the points on top, dodged to align perfectly
  geom_point(position = position_jitterdodge(jitter.width = 0.15, 
                                             dodge.width = 0.6), 
             size = 4.5, alpha = 0.9) +
  # Enforce the Viridis color scheme
  scale_color_viridis_d(option = "viridis") + 
  # Format the X-axis to show 10^x
 scale_x_continuous(     
   breaks = c(2, 4, 6, 8, 12),     
   labels = function(x) parse(text = paste0("10^", x))
   ) +
  scale_y_continuous(labels = math_format(10^.x)) +
  labs(
    title = "Group Means Across Condition Number",
    x = expression("Condition Number (" * kappa * ")"),
    y = expression("Mean" ~ e[lambda*",rel"]),
    color = "Dimension (n)", 
    shape = "Scale (s)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

ggsave("R_error_mean_plot.png", plot = p_mean, width = 10, height = 6, 
       dpi = 300, bg = "white")
print(p_mean)


# =====================================================================
# Grand Mean 
# =====================================================================

# We group by condition number and dimension, completely ignoring the scales.
df_grand_mean <- df_mean %>%
  group_by(cond_input, dim_fct) %>%
  summarise(grand_mean_error = mean(mean_log_error), .groups = 'drop')

p_grand <- ggplot(df_grand_mean, aes(x = log10(cond_input), 
                                     y = grand_mean_error, 
                                     color = dim_fct)) +

  geom_point(position = position_dodge(width = 0.5), size = 4.5, alpha = 0.9) +
  geom_line(aes(group = dim_fct), position = position_dodge(width = 0.5), 
            linewidth = 1) +
  scale_color_viridis_d(option = "viridis") + 
  scale_x_continuous(     
   breaks = c(2, 4, 6, 8, 12),     
   labels = function(x) parse(text = paste0("10^", x))
  ) +
  scale_y_continuous(labels = math_format(10^.x)) +
  labs(
    title = "Effect of Averaging Across Scaling Factor",
    x = expression("Condition Number (" * kappa * ")"),
    y = expression("Grand Mean "~ e[lambda*",rel"]),
    color = "Dimension (n)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

ggsave("R_flawed_grand_mean_plot.png", plot = p_grand, width = 10, height = 6, 
       dpi = 300, bg = "white")
print(p_grand)

# =======================================================================
# Continuous Linear model: Mean Relative Error 
# =======================================================================

#Fitting ANOVA model for average error
continuous_model <- lm(log10(eig_error_relative) ~ log10(cond_input) 
                       * log10(dim)
                       * log10(scale), data = df)

# upper and lower tell the stepwise specific models list. 
print("=== Running Stepwise Selection: Mean Error Model ===")
final_model <- step(continuous_model, 
                    scope = list(
                      lower = "~ log10(cond_input) + log10(dim) + log10(scale)",
                      upper = "~ log10(cond_input) * log10(dim) * log10(scale)"
                    ) )

print("=== Final Optimized Model ===")
summary(final_model)

#check residuals
plot(final_model)

#Fit the final model to the plot of relative error
pred_grid_cont <- expand.grid(
  cond_input = unique(df$cond_input),
  dim = unique(df$dim),
  scale = unique(df$scale)
) %>%
  mutate(
    dim_fct = as.factor(dim),
    scale_fct = as.factor(scale)
  )

pred_grid_cont$predicted_mean <- predict(final_model, newdata = pred_grid_cont)

p_mean_fit <- ggplot() +
  geom_point(data = df, 
             aes(x = log10(cond_input), 
                 y = log10(eig_error_relative), 
                 color = as.factor(dim), shape = as.factor(scale)), 
             position = position_jitterdodge(jitter.width = 0.15, 
                                             dodge.width = 0.6), 
             size = 2.5, alpha = 0.5) +
  
  geom_line(data = pred_grid_cont, 
            aes(x = log10(cond_input), y = predicted_mean, 
                group = interaction(dim_fct, scale_fct), 
                color = dim_fct,
                linetype = scale_fct), 
            linewidth = 1.2) +
  
  scale_color_viridis_d(option = "viridis") + 
  scale_x_continuous(     
   breaks = c(2, 4, 6, 8, 12),     
   labels = function(x) parse(text = paste0("10^", x))
  ) +
  scale_y_continuous(labels = math_format(10^.x)) +
  labs(
    title = "Fitted Linear Model For Means",
    x = expression("Condition Number (" * kappa * ")"),
    y = expression(e[lambda*",rel"]),
    color = "Dimension (n)", 
    shape = "Scale (s)",      
    linetype = "Scale (s)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

ggsave("R_continuous_model_fit.png", plot = p_mean_fit, width = 10, height = 6, 
       dpi = 300, bg = "white")
print(p_mean_fit)


# =======================================================================
# Continuous Linear Model: Variability
# =======================================================================

# Stepwise Selection for Variability Analysis
full_std_model <- lm(std_log_error ~ log10(cond_input) * log10(dim) 
                     * log10(scale), data = df_std)

print("=== Running Stepwise Selection: Variability Model ===")
best_std_model <- step(full_std_model, 
                       scope = list(
                         lower = "~ log10(cond_input) + log10(dim) + log10(scale)",
                         upper = "~ log10(cond_input) * log10(dim) * log10(scale)"
                       ))

print("=== Final Optimized Variability Model ===")
summary(best_std_model)

# generate predictions using the selected BEST model
df_std$predicted_std <- predict(best_std_model, newdata = df_std)

plot(best_std_model)

#### Final Variability model 
# Create a grid of all combinations you want to predict
pred_grid <- expand.grid(
  cond_input = unique(df_std$cond_input),
  dim = unique(df_std$dim),
  scale = unique(df_std$scale)
) %>%
  mutate(
    # match the variable names in your 'best_std_model'
    log10_cond = log10(cond_input),
    dim_fct = as.factor(dim),
    scale_fct = as.factor(scale)
  )

pred_grid$predicted_std <- predict(best_std_model, newdata = pred_grid)

p_std_fit <- ggplot() +
  # Layer 1: The raw data (from df_std)
  geom_point(data = df_std, 
             aes(x = log10(cond_input), y = std_log_error, 
                 color = dim_fct, shape = scale_fct), 
             position = position_jitterdodge(jitter.width = 0.15, 
                                             dodge.width = 0.6), 
             size = 4, alpha = 0.5) +
  
  # Layer 2: The model lines (from pred_grid)
  geom_line(data = pred_grid, 
            aes(x = log10(cond_input), y = predicted_std, 
                group = interaction(dim_fct, scale_fct), 
                color = dim_fct, linetype = scale_fct), 
            linewidth = 1.2) +
  
  scale_color_viridis_d(option = "viridis") + 
  scale_x_continuous(     
   breaks = c(2, 4, 6, 8, 12),
   labels = function(x) parse(text = paste0("10^", x))   
   ) +
  labs(
    title = "Fitted Linear Model for Variability",
    x = expression("Condition Number (" * kappa * ")"),
    y = expression(sigma),
    color = "Dimension (n)", 
    shape = "Scale (s)",
    linetype = "Scale (s)"
  ) +
  theme_minimal(base_size = 14)

ggsave("R_continuous_var_model_fit.png", plot = p_std_fit, width = 10, 
       height = 6,
       dpi = 300, bg = "white")
print(p_std_fit)

# =====================================================================
# Indicator Model
# =====================================================================

# Creating a TRUE/FALSE column, treating it as 1 or 0 (dummy variable)
df$is_failure <- ifelse(round(log10(df$cond_input)) == 12, 1, 0)

# The scale behaves differently in condition number 10^12, thus we will
# include the interaction in the model

indicator_model <- lm(log10(eig_error_relative) ~ log10(dim) + 
                        log10(cond_input) * log10(scale) + is_failure + 
                        is_failure:log10(scale), data = df)

print("=== Running Stepwise Selection on Indicator Model ===")
best_indicator_model <- step(indicator_model, direction = "both", 
                             trace = 1)

print("=== Final Optimized Indicator Model ===")
summary(best_indicator_model)

### Graph Indicator model
pred_grid_mean <- expand.grid(
  cond_input = unique(df$cond_input),
  dim = unique(df$dim),
  scale = unique(df$scale)
) %>%
  mutate(
    is_failure = ifelse(round(log10(cond_input)) == 12, 1, 0),
    dim_fct = as.factor(dim),
    scale_fct = as.factor(scale)
  )

pred_grid_mean$predicted_mean <- predict(best_indicator_model, pred_grid_mean)

p_indicator <- ggplot() +
  # THE DOTS: Raw data from 'df'
  geom_point(data = df, 
             aes(x = log10(cond_input), y = log10(eig_error_relative), 
                 color = as.factor(dim), shape = as.factor(scale)), 
             position = position_jitterdodge(jitter.width = 0.15, 
                                             dodge.width = 0.6), 
             size = 2.5, alpha = 0.5) +
  
  # THE MODEL LINES: Structural break lines from 'pred_grid_mean'
  # Mapping both shape and linetype to scale_fct with the same name merges them
  geom_line(data = pred_grid_mean, 
            aes(x = log10(cond_input), y = (predicted_mean), 
                group = interaction(dim_fct, scale_fct), 
                color = dim_fct, linetype = scale_fct), 
            linewidth = 1.2) +
  
  scale_color_viridis_d(option = "viridis") + 
  scale_x_continuous(
   breaks = c(2, 4, 6, 8, 12),
   labels = function(x) parse(text = paste0("10^", x))
   ) +
  scale_y_continuous(labels = scales::label_math(10^.x)) +
  labs(
    title = "Fitted Linear Model for Structural Break",
    x = expression("Condition Number (" * kappa * ")"),
    y = expression(e[lambda*",rel"]),
    color = "Dimension (n)", 
    shape = "Scale (s)",      
    linetype = "Scale (s)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

ggsave("Indicator_Model_Plot.png", plot = p_indicator, width = 10, height = 5, 
       bg = "white")
print(p_indicator)

plot(best_indicator_model)

# =====================================================================
# Running F-Test to compare Indicator vs. Non-Indicator Model
# =====================================================================

# Define the reduced model (Continuous only, including log(dim))
reduced_model <- lm(log10(eig_error_relative) ~ log10(cond_input) * 
                      log10(scale) + log10(dim), data = df)

# Define the full model (Indicator + Continuous, including log(dim))
full_indicator_model <- lm(log10(eig_error_relative) ~ log10(dim) + 
                        log10(cond_input) * log10(scale) + is_failure + 
                        is_failure:log10(scale), data = df)

# Compare models using an ANOVA F-test
print("=== F-Test: Continuous Model vs. Structural Break (Indicator) Model ===")
anova(reduced_model, full_indicator_model)

# =======================================================================
# Weighted Regression model: Mean Relative Error (Weighted by Dimension)
# =======================================================================

# Weighting by matrix dimension gives larger matrices more influence,
# since larger matrices may produce more numerically stable estimates.
weighted_continuous_model <- lm(log10(eig_error_relative) ~ log10(cond_input) 
                                * log10(dim) 
                                * log10(scale), 
                                data = df, 
                                weights = dim)

print("=== Running Stepwise Selection: Weighted Mean Error Model ===")
final_weighted_model <- step(weighted_continuous_model, 
                             scope = list(
                               lower = "~ log10(cond_input) + 
                               log10(dim) + log10(scale)",
                               upper = "~ log10(cond_input) * 
                               log10(dim) * log10(scale)"
                             ))

print("=== Final Optimized Weighted Model ===")
summary(final_weighted_model)

# Check residuals of the weighted model
plot(final_weighted_model)
