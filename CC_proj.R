library(tidyverse)
library(scales)
library(emmeans)

# =======================================================================
# Data Ingestion & Preparation
# =======================================================================

df <- read_csv("matrix_dataset.csv", show_col_types = FALSE) %>%
  
  # Numerical failures (e.g., overflow in eigenvalue solvers) produce Inf or -Inf.
  # These would break subsequent log-transformations and linear modelling,
  # so we remove them to keep only physically meaningful, finite computations.

  filter(is.finite(eig_error_relative), is.finite(cond_input)) %>%
  rename(dim = size) %>%
  
  # Many continuous variables are duplicated as factors so that ggplot can map
  # them to discrete colour/shape scales, enabling clear visual separation
  # between distinct experimental conditions rather than a continuous gradient.
  
  mutate(
    dim_fct = as.factor(dim),
    scale_fct = as.factor(scale),
    scale_label = paste0("Scale (s) = ", scale)
  )

# =======================================================================
# Simulation Results: Algorithm Error
# =======================================================================

p1 <- ggplot(df, aes(
                     x = log10(cond_input), 
                     y = log10(eig_error_relative), 
                     color = dim_fct, shape = scale_fct)
             ) +
  
  # Thousands of points are generated at the same discrete (cond, dim) combos.
  # Without jitter and dodge, points would overlap completely, hiding the
  # density and distribution of errors. Jitterdodge spreads them so the
  # reader can assess both central tendency and spread at each condition.
  
  geom_point(position = position_jitterdodge(jitter.width = 0.15, 
                    dodge.width = 0.6), size = 2.5, alpha = 0.6) + scale_color_viridis_d(option = "viridis") +
  
  # Data are plotted on log‑transformed axes to linearise exponential trends.
  # Axis breaks are labelled in power‑of‑10 notation (10^2, 10^4, …) because
  # readers intuitively understand matrix condition numbers in these orders.
  
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

# Standard deviation of the log‑error captures the “volatility” of the
# algorithm’s output. High σ means the order‑of‑magnitude accuracy swings
# wildly across replicates – a warning sign of numerical instability.

df_std <- df %>%
  group_by(cond_input, dim, scale, dim_fct, scale_fct) %>%
  summarise(std_log_error = sd(log10(eig_error_relative)), .groups = 'drop')

p2 <- ggplot(df_std, aes(x = log10(cond_input), y = std_log_error, 
                         color = dim_fct, shape = scale_fct)) +
  
  # Now each point is a summary per group, so overplotting is minimal.
  # Using dodge (without jitter) prevents shape occlusion while preserving
  # the exact position, making it easier to compare variability across groups.
  
  geom_point(position = position_dodge(width = 0.5), size = 4, alpha = 0.8) +
  scale_color_viridis_d(option = "viridis") + 
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

# Averaging on the log scale gives the geometric mean of the relative error.
# This prevents a single catastrophic outlier
# from dominating the central tendency of an otherwise well‑behaved group.

df_mean <- df %>%
  group_by(cond_input, dim, scale, dim_fct, scale_fct) %>%
  summarise(mean_log_error = mean(log10(eig_error_relative)), .groups = 'drop')

p_mean <- ggplot(df_mean, aes(x = log10(cond_input), 
                              y = mean_log_error, 
                              color = dim_fct, shape = scale_fct)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, 
                                             dodge.width = 0.6), 
             size = 4.5, alpha = 0.9) +
  scale_color_viridis_d(option = "viridis") + 
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
  # Minor gridlines are removed because the log‑scale axes already guide the eye;
  # extra lines only add visual noise without aiding interpretation.
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

ggsave("R_error_mean_plot.png", plot = p_mean, width = 10, height = 6, 
       dpi = 300, bg = "white")
print(p_mean)


# =====================================================================
# Grand Mean 
# =====================================================================

# Ignoring the scale factor and averaging over it reveals whether dimension
# alone interacts with condition number. However, because the scale factor
# can alter the error behaviour, this “grand mean” can be misleading
# (hence the cautious filename).

df_grand_mean <- df_mean %>%
  group_by(cond_input, dim_fct) %>%
  summarise(grand_mean_error = mean(mean_log_error), .groups = 'drop')

p_grand <- ggplot(df_grand_mean, aes(x = log10(cond_input), 
                                     y = grand_mean_error, 
                                     color = dim_fct)) +

  geom_point(position = position_dodge(width = 0.5), size = 4.5, alpha = 0.9) +
  
  # Connecting lines are used only for the same dimension, helping the eye
  # follow the trend as condition number increases – a visual guide to the
  # interaction between dimension and condition number.
  
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

# We include all two‑ and three‑way interactions because theory suggests
# that the effects of condition number, dimension, and scale may not be
# additive – e.g., ill‑conditioning may be amplified at large dimensions.
continuous_model <- lm(log10(eig_error_relative) ~ log10(cond_input) 
                       * log10(dim)
                       * log10(scale), data = df)

# We use stepwise selection based on AIC to prune out unnecessary interactions. 
# This prevents overfitting and leaves us with the most parsimonious model 
# that still explains the variance effectively.
print("=== Running Stepwise Selection: Mean Error Model ===")
final_model <- step(continuous_model, 
                    scope = list(
                      lower = "~ log10(cond_input) + log10(dim) + log10(scale)",
                      upper = "~ log10(cond_input) * log10(dim) * log10(scale)"
                    ))

print("=== Final Optimized Model ===")
summary(final_model)
plot(final_model)

# A synthetic grid of all unique experimental conditions is created so we can
# draw smooth, continuous prediction lines rather than jagged segments that
# merely connect the observed discrete points.

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
  
  # Model‑predicted lines overlay the raw data so the reader can judge
  # how well the fitted interactions capture the average behavior.
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

# A separate model for the standard deviation (heteroscedasticity) allows
# us to identify conditions where the algorithm’s output is not only
# inaccurate but also unpredictable.
full_std_model <- lm(std_log_error ~ log10(cond_input) * log10(dim) 
                     * log10(scale), data = df_std)

print("=== Running Stepwise Selection: Variability Model ===")
best_std_model <- step(full_std_model, 
                       scope = list(
                         lower = "~ log10(cond_input) + log10(dim) + 
                         log10(scale)",
                         upper = "~ log10(cond_input) * log10(dim) *
                         log10(scale)"
                       ))

print("=== Final Optimized Variability Model ===")
summary(best_std_model)

df_std$predicted_std <- predict(best_std_model, newdata = df_std)

plot(best_std_model)

pred_grid <- expand.grid(
  cond_input = unique(df_std$cond_input),
  dim = unique(df_std$dim),
  scale = unique(df_std$scale)
) %>%
  mutate(
    log10_cond = log10(cond_input),
    dim_fct = as.factor(dim),
    scale_fct = as.factor(scale)
  )

pred_grid$predicted_std <- predict(best_std_model, newdata = pred_grid)

p_std_fit <- ggplot() +
  geom_point(data = df_std, 
             aes(x = log10(cond_input), y = std_log_error, 
                 color = dim_fct, shape = scale_fct), 
             position = position_jitterdodge(jitter.width = 0.15, 
                                             dodge.width = 0.6), 
             size = 4, alpha = 0.5) +

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

# Floating‑point precision in matrix algorithms often collapses abruptly
# beyond a critical condition number. By introducing a binary failure indicator
# at κ = 10^12, the model can capture a structural break rather
# than forcing a smooth polynomial through two distinct regimes.

df$is_failure <- ifelse(round(log10(df$cond_input)) == 12, 1, 0)

# The interaction `is_failure:log10(scale)` is included because visual
# exploration suggests that different scaling factors react in opposite ways
# once the solver effectively fails – an effect that a main‑effect‑only
# indicator would miss.

indicator_model <- lm(log10(eig_error_relative) ~ log10(dim) + 
                        log10(cond_input) * log10(scale) + is_failure + 
                        is_failure:log10(scale), data = df)

print("=== Running Stepwise Selection on Indicator Model ===")
best_indicator_model <- step(indicator_model, direction = "both", 
                             trace = 1)

print("=== Final Optimized Indicator Model ===")
summary(best_indicator_model)

pred_grid_mean <- expand.grid(
  cond_input = unique(df$cond_input),
  dim = unique(df$dim),
  scale = unique(df$scale)
) %>%
  mutate(
    # The synthetic grid must map the exact same failure logic, otherwise 
    # the predictions will miss the structural break.
    is_failure = ifelse(round(log10(cond_input)) == 12, 1, 0),
    dim_fct = as.factor(dim),
    scale_fct = as.factor(scale)
  )

pred_grid_mean$predicted_mean <- predict(best_indicator_model, pred_grid_mean)

p_indicator <- ggplot() +
  geom_point(data = df, 
             aes(x = log10(cond_input), y = log10(eig_error_relative), 
                 color = as.factor(dim), shape = as.factor(scale)), 
             position = position_jitterdodge(jitter.width = 0.15, 
                                             dodge.width = 0.6), 
             size = 2.5, alpha = 0.5) +
  
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
# F-Tests to compare CMEM vs. CMEM+ and SBLM vs. CMEM+
# =====================================================================

# A formal F‑test on nested models quantifies whether the added parameters
# (failure indicator + its interaction with scale) significantly improve the
# model fit, justifying the extra complexity of the structural‑break approach.

# CMEM: Continuous Mean Error Model (AIC-selected model)
# Terms: dim, cond, scale, cond:scale
model_cmem <- lm(log10(eig_error_relative) ~ log10(dim) + 
                 log10(cond_input) * log10(scale), 
                 data = df)

# SBLM: Structural Break Linear Model (AIC-selected model)
# Terms: dim, cond, scale, failure, failure:scale 
# (Added log10(dim) back in to ensure perfect nesting)
model_sblm <- lm(log10(eig_error_relative) ~ log10(dim) + 
                   log10(cond_input) + log10(scale) + 
                   is_failure + is_failure:log10(scale), 
                 data = df)

# CMEM+: Minimal Supermodel (The Union)
# Terms: dim, cond, scale, cond:scale, failure, failure:scale
# This contains exactly the terms from CMEM and SBLM combined, and NOTHING else.
model_cmem_plus <- lm(log10(eig_error_relative) ~ log10(dim) + 
                        log10(cond_input) * log10(scale) + 
                        is_failure + is_failure:log10(scale), 
                      data = df)


# --- F-TEST COMPARISONS ---

# Test 1: CMEM vs. CMEM+
# Question: Do we need the structural break terms (is_failure and is_failure:scale)?
print("=== Test 1: Base Model (CMEM) vs. Minimal Supermodel (CMEM+) ===")
anova(model_cmem, model_cmem_plus)

# Test 2: SBLM vs. CMEM+
# Question: Do we need the continuous interaction term (cond_input:scale)?
print("=== Test 2: Selected Break Model (SBLM) vs. Minimal Supermodel (CMEM+) ===")
anova(model_sblm, model_cmem_plus)

# =======================================================================
# Weighted Regression model: Mean Relative Error (Weighted by Dimension)
# =======================================================================

# Larger matrices provide a greater number of internal floating‑point operations,
# often averaging out random noise and yielding more theoretically “stable”
# errors. Weighting the regression by dimension pulls the fitted line toward
# these more reliable observations.

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

plot(final_weighted_model)
